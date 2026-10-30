# Recoupling symbols: N, F (6j), R, A, B, and the 1j / 2j / 3j symbols.
#
# Layouts and conventions follow TensorKitSectors.jl, so that a Kiwi irrep can be
# used as a sector there (see the KiwiTensorKitSectorsExt extension):
#
#   fusion tensor   C^{ab}_c[i, j, k, μ] = ⟨a i; b j | c k; μ⟩     (see `fusion_tensor`)
#   F-symbol        F^{abc}_d[e, f][μ, ν, κ, λ]
#                     = Σ  C^{af}_d[i, n, 1, λ] C^{bc}_f[j, k, n, κ] C^{ab}_e[i, j, m, μ] C^{ec}_d[m, k, 1, ν]
#                   (the coefficient of the vertex pair (κ, λ) in the F-move of (μ, ν))
#   R-symbol        R^{ab}_c[μ, ν] = Σ C^{ba}_c[j, i, 1, ν] C^{ab}_c[i, j, 1, μ]
#
# The index 1 of the coupled irrep is its highest-weight vector; by Schur's lemma any
# basis vector gives the same F and R.  Kiwi evaluates these contractions directly
# on highest-weight vectors, without forming the full Clebsch–Gordan tensors.
#
# 1j and 2j symbols
# -----------------
# The 2j symbol of `a` is the unitary matrix
#     Z_a[i, j] = √dim(a) · C^{a ā}_1[i, j, 1, 1],
# i.e. the normalised invariant of V_a ⊗ V_ā scaled to be orthogonal.  Kiwi's
# convention fixes Z_a[hw(a), lw(ā)] > 0 (the multiplicity convention of
# `intertwiner` applied to a ⊗ ā → 1).  By Schur's lemma Z_ā = ϰ_a Z_aᵀ with
# ϰ_a = ±1, the 1j symbol (`onej`).  For self-conjugate a, ϰ_a is the
# Frobenius–Schur indicator (+1 real, −1 pseudoreal); for a ≠ ā, ϰ_a = ϰ_ā is a
# convention-dependent sign.  ϰ_a equals TensorKit's Frobenius–Schur phase
# sign(F^{a ā a}_a[1, 1]), which ties the line-bending (A, B) symbols to Z.
#
# 3j symbol
# ---------
#     (a b c)[i, j, k, μ] = Σ_m C^{ab}_{c̄}[i, j, m, μ] Z_c[k, m] / √dim(c),
# the invariant of V_a ⊗ V_b ⊗ V_c obtained by closing a ⊗ b → c̄ with the 2j
# symbol; Σ |(a b c)|² = 1 for each μ.  For su(2) it is related to Wigner's 3j symbol
# by (j₁ j₂ j₃)_Wigner = (−1)^{j₁−j₂−j₃} (j₁ j₂ j₃)_Kiwi.

const _NSYMBOL_CACHE = Dict{Tuple{LieAlgebra, Vector{Int}, Vector{Int}}, Dict{Vector{Int}, Int}}()
const _NSYMBOL_LOCK = ReentrantLock()

function _fusion_rules(a::Irrep, b::Irrep)
    key = (a.algebra, a.dynkin_labels, b.dynkin_labels)
    found = lock(() -> get(_NSYMBOL_CACHE, key, nothing), _NSYMBOL_LOCK)
    found === nothing || return found
    tp = tensor_product(a, b)
    d = Dict{Vector{Int}, Int}(ir.dynkin_labels => Int(m) for (ir, m) in tp.components)
    lock(() -> (_NSYMBOL_CACHE[key] = d), _NSYMBOL_LOCK)
    return d
end

"""
    nsymbol(a, b, c) -> Int

Multiplicity of `c` in `a ⊗ b`.
"""
nsymbol(a::Irrep, b::Irrep, c::Irrep) = get(_fusion_rules(a, b), c.dynkin_labels, 0)

_unit(g::LieAlgebra) = Irrep(g, zeros(Int, lie_rank(g)))

# ---------------------------------------------------------------------------
# Vectors of V_a ⊗ V_b ⊗ V_c of a fixed weight, as blocks keyed by weight indices
# ---------------------------------------------------------------------------

struct _TripleVectors{T}
    keys::Vector{NTuple{3, Int}}
    off::Dict{NTuple{3, Int}, Int}
    total::Int
end

function _triple_layout(Ma, Mb, Mc, keys)
    ks = sort!(collect(keys))
    off = Dict{NTuple{3, Int}, Int}()
    t = 0
    for k in ks
        off[k] = t
        t += Ma.ws.mult[k[1]] * Mb.ws.mult[k[2]] * Mc.ws.mult[k[3]]
    end
    return ks, off, t
end

# L_{μν} = (φ^{ab→e}_μ ⊗ 1)(hw vector ν of e ⊗ c → d), for all μ, ν, as blocks
function _left_vectors(A::Intertwiner{T}, B::Intertwiner{T}) where {T}
    Ma, Mb, Mc = A.Ma, A.Mb, B.Mb
    LB = B.layouts[1]
    XB = B.phi[1]
    out = Dict{NTuple{3, Int}, Array{T, 5}}()      # (ka,kb,kc) => ma × mb × mc × N_A × N_B
    for p in eachindex(LB.ka)
        ke, kc = LB.ka[p], LB.kb[p]
        me, mc = B.Ma.ws.mult[ke], Mc.ws.mult[kc]
        ΦA = _phi(A, ke)                             # total_ab × me × N_A
        LA = _layout(A, ke)
        for ν in 1:B.N
            Bblk = reshape(view(XB, _blockrange(LB, p, me, mc), 1, ν), me, mc)
            for μ in 1:A.N
                P = view(ΦA, :, :, μ) * Bblk           # total_ab × mc
                for q in eachindex(LA.ka)
                    ka, kb = LA.ka[q], LA.kb[q]
                    ma, mb = Ma.ws.mult[ka], Mb.ws.mult[kb]
                    blk = get!(() -> zeros(T, ma, mb, mc, A.N, B.N), out, (ka, kb, kc))
                    blk[:, :, :, μ, ν] .+= reshape(P[_blockrange(LA, q, ma, mb), :], ma, mb, mc)
                end
            end
        end
    end
    return out
end

# R_{κλ} = (1 ⊗ φ^{bc→f}_κ)(hw vector λ of a ⊗ f → d), for all κ, λ, as blocks
function _right_vectors(C::Intertwiner{T}, D::Intertwiner{T}) where {T}
    Ma, Mb, Mc = D.Ma, C.Ma, C.Mb
    LD = D.layouts[1]
    XD = D.phi[1]
    out = Dict{NTuple{3, Int}, Array{T, 5}}()      # (ka,kb,kc) => ma × mb × mc × N_C × N_D
    for p in eachindex(LD.ka)
        ka, kf = LD.ka[p], LD.kb[p]
        ma, mf = Ma.ws.mult[ka], D.Mb.ws.mult[kf]
        ΦC = _phi(C, kf)                             # total_bc × mf × N_C
        LC = _layout(C, kf)
        for λ in 1:D.N
            Dblk = reshape(view(XD, _blockrange(LD, p, ma, mf), 1, λ), ma, mf)
            for κ in 1:C.N
                Q = view(ΦC, :, :, κ) * transpose(Dblk)   # total_bc × ma
                for q in eachindex(LC.ka)
                    kb, kc = LC.ka[q], LC.kb[q]
                    mb, mc = Mb.ws.mult[kb], Mc.ws.mult[kc]
                    blk = get!(() -> zeros(T, ma, mb, mc, C.N, D.N), out, (ka, kb, kc))
                    blk[:, :, :, κ, λ] .+= permutedims(reshape(Q[_blockrange(LC, q, mb, mc), :], mb, mc, ma), (3, 1, 2))
                end
            end
        end
    end
    return out
end

function _metric_apply(M::IrrepModule{T}, k::Int, x::AbstractMatrix) where {T}
    isorthonormal(M) ? x : M.gram[k] * x
end

# Σ over blocks of ⟨R_{κλ}, L_{μν}⟩ (product metric), as an N_A × N_B × N_C × N_D array
function _overlap(Ma, Mb, Mc, L::Dict, R::Dict, NA, NB, NC, ND, ::Type{T}) where {T}
    F = zeros(T, NA, NB, NC, ND)
    for (key, Lb) in L
        Rb = get(R, key, nothing)
        Rb === nothing && continue
        ka, kb, kc = key
        ma, mb, mc = size(Lb, 1), size(Lb, 2), size(Lb, 3)
        Lm = reshape(Lb, ma * mb * mc, NA * NB)
        if !(isorthonormal(Ma) && isorthonormal(Mb) && isorthonormal(Mc))
            Ga = isorthonormal(Ma) ? Matrix{T}(I, ma, ma) : Ma.gram[ka]
            Gb = isorthonormal(Mb) ? Matrix{T}(I, mb, mb) : Mb.gram[kb]
            Gc = isorthonormal(Mc) ? Matrix{T}(I, mc, mc) : Mc.gram[kc]
            Lm = kron(Gc, kron(Gb, Ga)) * Lm
        end
        Rm = reshape(Rb, ma * mb * mc, NC * ND)
        F .+= permutedims(reshape(transpose(Rm) * Lm, NC, ND, NA, NB), (3, 4, 1, 2))
    end
    return F
end

# ---------------------------------------------------------------------------
# Raw (unnormalised) and normalised symbols
#
# Raw gauge.  With an exact scalar type (Rational{BigInt}, Fp{P}) the modules use the
# unnormalised lowering basis and the intertwiners X_μ of each vertex are mutually
# orthogonal but unnormalised, with squared norms ⟨X_μ, X_μ⟩ = `norms2[μ]` (the theta
# value of the vertex).  The raw symbols are the same contractions as the normalised
# ones, evaluated on these vertices: every entry is a rational function of the data,
# with no square roots.  The normalised symbol is recovered by dividing each entry by
# the square root of the product of the norms² of the vertices it involves.  With a
# floating-point type the vertices are orthonormal and all norms² are one.
# ---------------------------------------------------------------------------

_check_normalisable(T) = (T <: Fp) &&
    throw(ArgumentError("normalised symbols are not defined over a prime field; use normalized = false"))

# entrywise s·√(x²/n) as SqrtRational (exact) or x/√n (floats)
_normalise(x, n, ::Type{T}) where {T} =
    _isfloat(T) ? x / sqrt(n) : SqrtRational(sign(x), Rational{BigInt}(x)^2 / Rational{BigInt}(n))

"""
    fsymbol(a, b, c, d, e, f; basis = :lowering, T = Float64, normalized = true)

The F-symbol (6j symbol) `F^{abc}_d[e, f]` in TensorKit's convention, an array of
size `(N^{ab}_e, N^{ec}_d, N^{bc}_f, N^{af}_d)` relating the fusion trees
((a ⊗ b → e) ⊗ c → d) and (a ⊗ (b ⊗ c → f) → d).

* `normalized = true`: unitary F-symbol.  Floating-point `T` gives an `Array{T, 4}`;
  `T = Rational{BigInt}` gives exact entries as [`SqrtRational`](@ref) (any
  multiplicity).
* `normalized = false`: the raw F-symbol in the rational gauge, for any `T`
  (including `Fp{P}`), as a named tuple `(F, norms)` with `F::Array{T, 4}` and
  `norms = (nᵃᵇₑ, nᵉᶜ_d, nᵇᶜ_f, nᵃᶠ_d)` the squared norms of the four vertices.  Then
  `F_unitary[μ, ν, κ, λ] = F[μ, ν, κ, λ] / √(nᵃᵇₑ[μ] nᵉᶜ_d[ν] nᵇᶜ_f[κ] nᵃᶠ_d[λ])`.
"""
function fsymbol(a::Irrep, b::Irrep, c::Irrep, d::Irrep, e::Irrep, f::Irrep;
                 basis::Symbol = :lowering, T::Type = Float64, normalized::Bool = true)
    key = (:F, a.algebra, a.dynkin_labels, b.dynkin_labels, c.dynkin_labels, d.dynkin_labels,
           e.dynkin_labels, f.dynkin_labels, basis, T)
    raw = _memo(() -> _fsymbol_raw(a, b, c, d, e, f, basis, T), key)
    normalized || return raw
    _isfloat(T) && return raw.F
    _check_normalisable(T)
    nA, nB, nC, nD = raw.norms
    return [_normalise(raw.F[μ, ν, κ, λ], nA[μ] * nB[ν] * nC[κ] * nD[λ], T)
            for μ in axes(raw.F, 1), ν in axes(raw.F, 2), κ in axes(raw.F, 3), λ in axes(raw.F, 4)]
end

function _fsymbol_raw(a, b, c, d, e, f, basis, T)
    NA, NB, NC, ND = nsymbol(a, b, e), nsymbol(e, c, d), nsymbol(b, c, f), nsymbol(a, f, d)
    if NA * NB * NC * ND == 0
        return (F = zeros(T, NA, NB, NC, ND), norms = (ones(T, NA), ones(T, NB), ones(T, NC), ones(T, ND)))
    end
    A = intertwiner(a, b, e; basis, T); B = intertwiner(e, c, d; basis, T)
    C = intertwiner(b, c, f; basis, T); D = intertwiner(a, f, d; basis, T)
    L = _left_vectors(A, B)
    R = _right_vectors(C, D)
    F = _overlap(A.Ma, A.Mb, B.Mb, L, R, NA, NB, NC, ND, T)
    return (F = F, norms = (A.norms2, B.norms2, C.norms2, D.norms2))
end

"""
    rsymbol(a, b, c; basis = :lowering, T = Float64, normalized = true)

The R-symbol `R^{ab}_c[μ, ν]` (TensorKit convention): the swap V_a ⊗ V_b → V_b ⊗ V_a
maps the intertwiner μ of `c → a ⊗ b` to `Σ_ν R[μ, ν]` times the intertwiner ν of
`c → b ⊗ a` (the vertex swap).

`normalized = false` returns `(R, norms)` in the rational gauge, for any `T`, with
`norms = (nᵃᵇ_c, nᵇᵃ_c)` and `R_unitary[μ, ν] = R[μ, ν] / √(nᵃᵇ_c[μ] nᵇᵃ_c[ν])`.
With `T = Rational{BigInt}` and `normalized = true` the entries are `SqrtRational`.
"""
function rsymbol(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64,
                 normalized::Bool = true)
    key = (:R, a.algebra, a.dynkin_labels, b.dynkin_labels, c.dynkin_labels, basis, T)
    raw = _memo(() -> _rsymbol_raw(a, b, c, basis, T), key)
    normalized || return raw
    _isfloat(T) && return raw.R
    _check_normalisable(T)
    nA, nB = raw.norms
    return [_normalise(raw.R[μ, ν], nA[μ] * nB[ν], T) for μ in axes(raw.R, 1), ν in axes(raw.R, 2)]
end

function _rsymbol_raw(a, b, c, basis, T)
    N = nsymbol(a, b, c)
    N == 0 && return (R = zeros(T, 0, 0), norms = (T[], T[]))
    A = intertwiner(a, b, c; basis, T)
    B = intertwiner(b, a, c; basis, T)
    LA, LB = A.layouts[1], B.layouts[1]
    XA, XB = A.phi[1], B.phi[1]
    R = zeros(T, N, N)
    for p in eachindex(LA.ka)
        ka, kb = LA.ka[p], LA.kb[p]
        q = get(LB.block, (kb, ka), 0)
        q == 0 && continue
        ma, mb = A.Ma.ws.mult[ka], A.Mb.ws.mult[kb]
        for μ in 1:N, ν in 1:N
            Ablk = reshape(view(XA, _blockrange(LA, p, ma, mb), 1, μ), ma, mb)
            Bblk = reshape(view(XB, _blockrange(LB, q, mb, ma), 1, ν), mb, ma)
            Ga = isorthonormal(A.Ma) ? I : A.Ma.gram[ka]
            Gb = isorthonormal(A.Mb) ? I : A.Mb.gram[kb]
            R[μ, ν] += tr(Ablk * (Gb * Bblk * Ga))
        end
    end
    return (R = R, norms = (A.norms2, B.norms2))
end

# full (block-diagonal) Gram matrix of a module's basis; the identity if orthonormal
function _full_gram(M::IrrepModule{T}) where {T}
    d = dimension(M)
    G = zeros(T, d, d)
    ws = M.ws
    for k in 1:length(ws)
        r = (ws.offset[k] + 1):(ws.offset[k] + ws.mult[k])
        G[r, r] = isorthonormal(M) ? Matrix{T}(I, ws.mult[k], ws.mult[k]) : M.gram[k]
    end
    return G
end

"""
    twoj(a; basis = :lowering, T = Float64, normalized = true)

The 2j symbol `Z_a = √dim(a) · C^{a ā}_1[:, :, 1, 1]`, a real orthogonal matrix
(rows: basis of `a`, columns: basis of `ā`) with `Z_a[hw(a), lw(ā)] > 0`.

`normalized = false` gives the cup in the rational gauge, for any `T`: a named tuple
`(X, norm2, gram_a, gram_abar)` where `X[i, j]` are the coordinates of the
(unnormalised) invariant vector of `V_a ⊗ V_ā` in the module bases, `norm2` its
squared norm and `gram_a`, `gram_abar` the Gram matrices of the bases.  With
`G = RᵀR` (Cholesky), `Z_a = √dim(a) · R_a X R_āᵀ / √norm2`.
"""
function twoj(a::Irrep; basis::Symbol = :lowering, T::Type = Float64, normalized::Bool = true)
    if !normalized
        X = intertwiner(a, conjugate(a), _unit(a.algebra); basis, T)
        C = _fusion_tensor(X)
        return (X = C[:, :, 1, 1], norm2 = X.norms2[1], gram_a = _full_gram(X.Ma), gram_abar = _full_gram(X.Mb))
    end
    _isfloat(T) || throw(ArgumentError("the normalised 2j symbol needs a floating-point type; use normalized = false"))
    C = fusion_tensor(a, conjugate(a), _unit(a.algebra); basis, T)
    return sqrt(T(dimension(a))) * C[:, :, 1, 1]
end

"""
    onej(a; basis = :lowering) -> Int

The 1j symbol ϰ_a = ±1 defined by `twoj(ā) = ϰ_a · transpose(twoj(a))`.  For
self-conjugate `a` this is the Frobenius–Schur indicator; it always equals
TensorKit's Frobenius–Schur phase `sign(F^{a ā a}_a[1, 1])`.  In the lowering basis
it is computed once, in exact arithmetic, and cached per irrep.
"""
function onej(a::Irrep; basis::Symbol = :lowering, T::Type = Float64)
    return _memo(() -> _onej(a, basis), (:onej, a.algebra, a.dynkin_labels, basis))::Int
end

function _onej(a, basis)
    # basis-dependent sign; exact in the lowering basis, floating point for GT
    T = basis == :lowering ? Rational{BigInt} : Float64
    ab = conjugate(a)
    Za = twoj(a; basis, T, normalized = false)
    Zab = twoj(ab; basis, T, normalized = false)
    # ⟨cup(ā), swap(cup(a))⟩ in the metric of V_ā ⊗ V_a
    s = tr(transpose(Zab.X) * Zab.gram_a * transpose(Za.X) * Za.gram_a)
    _isfloat(T) && abs(abs(s) / sqrt(Za.norm2 * Zab.norm2) - 1) > 1e-8 &&
        error("2j symbols of $a and its conjugate are not related by a sign")
    return s > 0 ? 1 : -1
end

"""
    threej(a, b, c; basis = :lowering, T = Float64) -> Array{T, 4}

The 3j symbol `(a b c)[i, j, k, μ] = Σ_m C^{ab}_{c̄}[i, j, m, μ] Z_c[k, m] / √dim(c)`, an
invariant tensor of `V_a ⊗ V_b ⊗ V_c` normalised to `Σ |·|² = 1` for each
multiplicity label μ of `a ⊗ b → c̄`.  For su(2),
`wigner3j(j₁, j₂, j₃) = (−1)^{j₁−j₂−j₃} · threej`.
"""
function threej(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64)
    cb = conjugate(c)
    C = fusion_tensor(a, b, cb; basis, T)
    Z = twoj(c; basis, T)
    da, db, dcb, N = size(C)
    out = zeros(T, da, db, size(Z, 1), N)
    for μ in 1:N
        M = reshape(view(C, :, :, :, μ), da * db, dcb) * transpose(Z)
        out[:, :, :, μ] = reshape(M, da, db, size(Z, 1)) ./ sqrt(T(dimension(c)))
    end
    return out
end

"""
    asymbol(a, b, c; basis = :lowering, T = Float64, normalized = true)

TensorKit's A-symbol `A^{ab}_c = √(dim a dim b / dim c) · ϰ_a F^{ā a b}_b[1, c]`,
reshaped to `(N^{ab}_c, N^{āc}_b)`: it relates the vertex `c → a ⊗ b` to the vertex
`b → ā ⊗ c` obtained by bending the leg `a`.

`normalized = false` returns, for any `T`, `(A, norms, scale2)` with
`A = ϰ_a · F_raw^{ā a b}_b[1, c]` (rational gauge),
`norms = (n_cup, n_unit, nᵃᵇ_c, nᵃ̄ᶜ_b)` (the vertices ā ⊗ a → 1, 1 ⊗ b → b, a ⊗ b → c,
ā ⊗ c → b) and `scale2 = dim a · dim b // dim c`, so that
`A_unitary[μ, ν] = √scale2 · A[μ, ν] / √(n_cup · n_unit · nᵃᵇ_c[μ] · nᵃ̄ᶜ_b[ν])`.
"""
function asymbol(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64,
                 normalized::Bool = true)
    ab = conjugate(a)
    shape = (nsymbol(a, b, c), nsymbol(ab, c, b))
    κ = onej(a; basis)
    scale2 = Rational{BigInt}(dimension(a)) * dimension(b) // dimension(c)
    if !normalized
        raw = fsymbol(ab, a, b, b, _unit(a.algebra), c; basis, T, normalized = false)
        n1, n2, n3, n4 = raw.norms
        return (A = κ .* reshape(raw.F, shape), norms = (n1[1], n2[1], n3, n4), scale2 = scale2)
    end
    F = fsymbol(ab, a, b, b, _unit(a.algebra), c; basis, T)
    _isfloat(T) && return reshape(sqrt(T(scale2)) * κ .* F, shape)
    return reshape([SqrtRational(κ * x.sign, x.square * scale2) for x in F], shape)
end

"""
    bsymbol(a, b, c; basis = :lowering, T = Float64, normalized = true)

TensorKit's B-symbol `B^{ab}_c = √(dim a dim b / dim c) · F^{a b b̄}_a[c, 1]`, reshaped
to `(N^{ab}_c, N^{c b̄}_a)`: it relates the vertex `c → a ⊗ b` to the vertex
`a → c ⊗ b̄` obtained by bending the leg `b`.

`normalized = false` returns, for any `T`, `(B, norms, scale2)` with
`B = F_raw^{a b b̄}_a[c, 1]`, `norms = (nᵃᵇ_c, nᶜᵇ̄_a, n_cup, n_unit)` (the vertices
a ⊗ b → c, c ⊗ b̄ → a, b ⊗ b̄ → 1, a ⊗ 1 → a) and `scale2 = dim a · dim b // dim c`:
`B_unitary[μ, ν] = √scale2 · B[μ, ν] / √(nᵃᵇ_c[μ] · nᶜᵇ̄_a[ν] · n_cup · n_unit)`.
"""
function bsymbol(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64,
                 normalized::Bool = true)
    bb = conjugate(b)
    shape = (nsymbol(a, b, c), nsymbol(c, bb, a))
    scale2 = Rational{BigInt}(dimension(a)) * dimension(b) // dimension(c)
    if !normalized
        raw = fsymbol(a, b, bb, a, c, _unit(a.algebra); basis, T, normalized = false)
        n1, n2, n3, n4 = raw.norms
        return (B = reshape(raw.F, shape), norms = (n1, n2, n3[1], n4[1]), scale2 = scale2)
    end
    F = fsymbol(a, b, bb, a, c, _unit(a.algebra); basis, T)
    _isfloat(T) && return reshape(sqrt(T(scale2)) .* F, shape)
    return reshape([SqrtRational(x.sign, x.square * scale2) for x in F], shape)
end

# ---------------------------------------------------------------------------
# Exact square roots of rationals
# ---------------------------------------------------------------------------

"""
    SqrtRational

An exact real number `s · √q` with `s ∈ {−1, 0, 1}` and `q ≥ 0` rational, the form
taken by multiplicity-free recoupling coefficients in exact arithmetic.
"""
struct SqrtRational <: Real
    sign::Int
    square::Rational{BigInt}
    function SqrtRational(s::Integer, q::Rational)
        q < 0 && throw(DomainError(q, "square must be nonnegative"))
        iszero(q) && return new(0, zero(Rational{BigInt}))
        return new(Int(sign(s)), Rational{BigInt}(q))
    end
end
SqrtRational(s, q) = SqrtRational(Int(sign(s)), Rational{BigInt}(q))
Base.zero(::Type{SqrtRational}) = SqrtRational(0, 0 // 1)
Base.float(x::SqrtRational) = x.sign * sqrt(Float64(x.square))
Base.Float64(x::SqrtRational) = float(x)
Base.BigFloat(x::SqrtRational) = x.sign * sqrt(BigFloat(x.square))
Base.convert(::Type{Float64}, x::SqrtRational) = float(x)
Base.:(==)(x::SqrtRational, y::SqrtRational) = x.sign == y.sign && x.square == y.square
Base.isapprox(x::SqrtRational, y::Real; kw...) = isapprox(float(x), y; kw...)
Base.isapprox(y::Real, x::SqrtRational; kw...) = isapprox(float(x), y; kw...)
Base.:-(x::SqrtRational) = SqrtRational(-x.sign, x.square)
function Base.show(io::IO, x::SqrtRational)
    x.sign == 0 && return print(io, "0")
    print(io, x.sign < 0 ? "-" : "", "√(", numerator(x.square), "/", denominator(x.square), ")")
end
