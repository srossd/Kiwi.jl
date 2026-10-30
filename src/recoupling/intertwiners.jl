# Intertwiners V_c → V_a ⊗ V_b (Clebsch–Gordan coefficients).
#
# Highest-weight vectors
# ----------------------
# A copy of V_c in V_a ⊗ V_b is fixed by its highest-weight vector X, a vector of
# weight c killed by every Eᵢ.  Write X = Σ X_{μ,ν} (μ + ν = c) over pairs of weight
# spaces.  Two facts make the computation cheap:
#
# 1. (Parthasarathy–Ranga Rao–Varadarajan.)  The component of X along v_a ⊗ V_b
#    is v_a ⊗ u with u ∈ V_b(c − a) satisfying Eᵢ^{⟨a,αᵢ^∨⟩+1} u = 0 for all i, and
#    X ↦ u is a bijection onto that subspace.  So the multiplicity space is a null
#    space inside one weight space of V_b.
# 2. The condition Eᵢ X = 0 at weight (μ + αᵢ, ν) reads
#        Eᵢ X_{μ,ν} + X_{μ+αᵢ, ν−αᵢ} Eᵢᵀ = 0      (matrix form, X_{μ,ν}: V_a(μ) × V_b(ν)).
#    For μ ≠ a the stacked map (Eᵢ)ᵢ on V_a(μ) is injective, so X_{μ,ν} follows from
#    the blocks above it with a (cached) left inverse.  Going down the weights of V_a
#    reconstructs X from u.
#
# Multiplicity convention
# -----------------------
# Let K be the space of highest-weight vectors (dim N) and f₁, f₂, … the basis of
# V_b(c − a).  The intertwiner basis X₁, …, X_N is the Gram–Schmidt
# orthonormalisation of the projections Π_K(v_a ⊗ f_j), j = 1, 2, … in order,
# skipping those dependent on earlier ones.  In particular ⟨X_k, v_a ⊗ f_{j_k}⟩ > 0,
# and in the multiplicity-free case the first nonzero coefficient of X along
# v_a ⊗ V_b(c − a) is positive.  For su(2) this is the Condon–Shortley convention
# ⟨j₁ j₁; j₂ (J − j₁) | J J⟩ > 0.
#
# Lowering
# --------
# The image φ(V_c(λ)) of a lower weight space follows from φ(Fᵢ v) = Δ(Fᵢ) φ(v):
# with S = [F₁ F₂ …] : ⊕ᵢ V_c(λ + αᵢ) → V_c(λ) (surjective for λ ≠ c),
# φ_λ = [Δ(Fᵢ) φ_{λ+αᵢ}]ᵢ S⁺ for any right inverse S⁺.

"""
    PairLayout

Layout of the weight-μ subspace of `V_a ⊗ V_b`: blocks `(ka, kb)` of weight-space
indices with `w_a + w_b = μ`, stored contiguously (each block column-major,
`mult_a × mult_b`), ordered by the depth of `ka`.
"""
struct PairLayout
    ka::Vector{Int}
    kb::Vector{Int}
    off::Vector{Int}
    total::Int
    block::Dict{Tuple{Int, Int}, Int}
end

function PairLayout(wa::WeightStructure, wb::WeightStructure, μ::AbstractVector{Int})
    pairs = Tuple{Int, Int}[]
    if length(wa) <= length(wb)
        for (ka, w) in enumerate(wa.weights)
            kb = get(wb.index, μ .- w, 0)
            kb == 0 || push!(pairs, (ka, kb))
        end
    else
        for (kb, w) in enumerate(wb.weights)
            ka = get(wa.index, μ .- w, 0)
            ka == 0 || push!(pairs, (ka, kb))
        end
        sort!(pairs)
    end
    off = Int[]; t = 0
    block = Dict{Tuple{Int, Int}, Int}()
    for (p, (ka, kb)) in enumerate(pairs)
        push!(off, t)
        block[(ka, kb)] = p
        t += wa.mult[ka] * wb.mult[kb]
    end
    return PairLayout(first.(pairs), last.(pairs), off, t, block)
end

@inline _blockrange(L::PairLayout, p::Int, ma::Int, mb::Int) = (L.off[p] + 1):(L.off[p] + ma * mb)

"""
    _inner(Ma, Mb, L, x, y)

Inner product of two vectors of `V_a ⊗ V_b` with layout `L` (product metric).
"""
function _inner(Ma::IrrepModule{T}, Mb::IrrepModule{T}, L::PairLayout, x::AbstractVector,
                y::AbstractVector) where {T}
    (isorthonormal(Ma) && isorthonormal(Mb)) && return dot(x, y)
    s = zero(T)
    for p in eachindex(L.ka)
        ma, mb = Ma.ws.mult[L.ka[p]], Mb.ws.mult[L.kb[p]]
        r = _blockrange(L, p, ma, mb)
        X = reshape(view(x, r), ma, mb); Y = reshape(view(y, r), ma, mb)
        Ga = isorthonormal(Ma) ? I : Ma.gram[L.ka[p]]
        Gb = isorthonormal(Mb) ? I : Mb.gram[L.kb[p]]
        s += tr(transpose(X) * (Ga * Y * Gb))
    end
    return s
end

# cached left inverse of the stacked raising operators on V(μ), μ ≠ highest weight
function _raise_left_inverse(M::IrrepModule{T}, k::Int) where {T}
    lock(M.lock) do
        cache = get!(() -> Dict{Int, Matrix{T}}(), M.cache, :raise_linv)::Dict{Int, Matrix{T}}
        get!(cache, k) do
            ws = M.ws
            A = reduce(vcat, (M.E[i, k] for i in 1:size(ws.up, 1) if ws.up[i, k] != 0))
            _left_inverse(A)
        end
    end
end

# cached right inverse of S = [F₁ F₂ …] : ⊕ᵢ V(μ + αᵢ) → V(μ)
function _lower_right_inverse(M::IrrepModule{T}, k::Int) where {T}
    lock(M.lock) do
        cache = get!(() -> Dict{Int, Matrix{T}}(), M.cache, :lower_rinv)::Dict{Int, Matrix{T}}
        get!(cache, k) do
            ws = M.ws
            S = reduce(hcat, (M.F[i, k] for i in 1:size(ws.up, 1) if ws.up[i, k] != 0))
            _right_inverse(S)
        end
    end
end

"""
    _highest_weight_vectors(Ma, Mb, c) -> (L, X, norms2)

Basis of the highest-weight vectors of weight `c` in `V_a ⊗ V_b`, in Kiwi's
multiplicity convention: columns of `X` (layout `L`).  Floating-point modules give
orthonormal columns (`norms2` all one); exact modules give mutually orthogonal,
unnormalised columns with squared norms `norms2`.
"""
function _highest_weight_vectors(Ma::IrrepModule{T}, Mb::IrrepModule{T}, c::Vector{Int}) where {T}
    wa, wb = Ma.ws, Mb.ws
    a = wa.irrep.dynkin_labels
    n = length(a)
    L = PairLayout(wa, wb, c)
    kb0 = get(wb.index, c .- a, 0)
    (kb0 == 0 || isempty(L.ka) || L.ka[1] != 1) && return L, zeros(T, L.total, 0), T[]
    mb0 = wb.mult[kb0]
    # 1. PRV conditions on u ∈ V_b(c − a)
    ops = Matrix{T}[]
    for i in 1:n
        op = Matrix{T}(I, mb0, mb0)
        cur = kb0
        alive = true
        for _ in 1:(a[i] + 1)
            nxt = wb.up[i, cur]
            if nxt == 0
                alive = false
                break
            end
            op = Mb.E[i, cur] * op
            cur = nxt
        end
        alive && push!(ops, op)
    end
    U = isempty(ops) ? Matrix{T}(I, mb0, mb0) : _nullspace(reduce(vcat, ops))
    N = size(U, 2)
    N == 0 && return L, zeros(T, L.total, 0), T[]
    # 2. propagate each u down the weights of V_a
    X = zeros(T, L.total, N)
    for k in 1:N
        x = view(X, :, k)
        x[_blockrange(L, 1, 1, mb0)] = U[:, k]
        for p in 2:length(L.ka)
            ka, kb = L.ka[p], L.kb[p]
            ma, mb = wa.mult[ka], wb.mult[kb]
            rhs = Matrix{T}[]
            for i in 1:n
                kau = wa.up[i, ka]
                kau == 0 && continue
                blk = zeros(T, wa.mult[kau], mb)
                kbd = wb.down[i, kb]
                if kbd != 0
                    q = get(L.block, (kau, kbd), 0)
                    if q != 0
                        Y = reshape(view(x, _blockrange(L, q, wa.mult[kau], wb.mult[kbd])),
                                    wa.mult[kau], wb.mult[kbd])
                        blk = -(Y * transpose(Mb.E[i, kbd]))
                    end
                end
                push!(rhs, blk)
            end
            x[_blockrange(L, p, ma, mb)] = vec(_raise_left_inverse(Ma, ka) * reduce(vcat, rhs))
        end
    end
    # 3. multiplicity convention: Gram–Schmidt of the projections of v_a ⊗ f_j
    G = [_inner(Ma, Mb, L, view(X, :, k), view(X, :, l)) for k in 1:N, l in 1:N]
    Gb = isorthonormal(Mb) ? Matrix{T}(I, mb0, mb0) : Mb.gram[kb0]
    t = transpose(U) * Gb                       # t[k, j] = ⟨X_k, v_a ⊗ f_j⟩
    coeffs = _solve(G, t)                       # column j: Π_K(v_a ⊗ f_j) in the X basis
    basis = Vector{Vector{T}}()
    norms2 = T[]
    for j in 1:mb0
        y = coeffs[:, j]
        for (z, nz) in zip(basis, norms2)
            y = y - ((transpose(z) * G * y) / nz) * z
        end
        r = transpose(y) * G * y
        _is_negligible(r, transpose(coeffs[:, j]) * G * coeffs[:, j], T) && continue
        if _isfloat(T)
            y = y / sqrt(r)
            r = one(T)
        end
        push!(basis, y)
        push!(norms2, r)
        length(basis) == N && break
    end
    length(basis) == N || error("could not fix the multiplicity basis")
    return L, X * reduce(hcat, basis), norms2
end

"""
    Intertwiner{T}

Basis of `Hom(V_c, V_a ⊗ V_b)` in Kiwi's convention, with its images of the weight
spaces of `V_c` computed on demand.
"""
struct Intertwiner{T}
    Ma::IrrepModule{T}
    Mb::IrrepModule{T}
    Mc::IrrepModule{T}
    N::Int
    norms2::Vector{T}
    layouts::Dict{Int, PairLayout}
    phi::Dict{Int, Array{T, 3}}            # weight index of c ↦ total × mult_c × N
    lock::ReentrantLock
end

Base.show(io::IO, X::Intertwiner{T}) where {T} =
    print(io, "Intertwiner{$T}(", X.Mc.ws.irrep, " → ", X.Ma.ws.irrep, " ⊗ ", X.Mb.ws.irrep,
          ", multiplicity ", X.N, ")")

function Intertwiner(Ma::IrrepModule{T}, Mb::IrrepModule{T}, Mc::IrrepModule{T}) where {T}
    c = Mc.ws.irrep.dynkin_labels
    L, X, norms2 = _highest_weight_vectors(Ma, Mb, c)
    N = size(X, 2)
    phi = Dict{Int, Array{T, 3}}()
    layouts = Dict{Int, PairLayout}(1 => L)
    N > 0 && (phi[1] = reshape(X, L.total, 1, N))
    return Intertwiner{T}(Ma, Mb, Mc, N, norms2, layouts, phi, ReentrantLock())
end

function _layout(X::Intertwiner, kc::Int)
    lock(X.lock) do
        get!(X.layouts, kc) do
            PairLayout(X.Ma.ws, X.Mb.ws, X.Mc.ws.weights[kc])
        end
    end
end

"""
    _apply_lower(Ma, Mb, i, Lsrc, Y, Ldst) -> Z

`Δ(Fᵢ) = Fᵢ ⊗ 1 + 1 ⊗ Fᵢ` applied to the columns of `Y` (layout `Lsrc`, weight
μ + αᵢ), giving columns with layout `Ldst` (weight μ).
"""
function _apply_lower(Ma::IrrepModule{T}, Mb::IrrepModule{T}, i::Int, Lsrc::PairLayout,
                      Y::AbstractMatrix{T}, Ldst::PairLayout) where {T}
    wa, wb = Ma.ws, Mb.ws
    m = size(Y, 2)
    Z = zeros(T, Ldst.total, m)
    for p in eachindex(Ldst.ka)
        ka, kb = Ldst.ka[p], Ldst.kb[p]
        ma, mb = wa.mult[ka], wb.mult[kb]
        r = _blockrange(Ldst, p, ma, mb)
        kau = wa.up[i, ka]
        q1 = kau == 0 ? 0 : get(Lsrc.block, (kau, kb), 0)
        kbu = wb.up[i, kb]
        q2 = kbu == 0 ? 0 : get(Lsrc.block, (ka, kbu), 0)
        (q1 == 0 && q2 == 0) && continue
        if q1 != 0
            r1 = _blockrange(Lsrc, q1, wa.mult[kau], mb)
            Fa = Ma.F[i, ka]
            for col in 1:m
                Zb = reshape(view(Z, r, col), ma, mb)
                mul!(Zb, Fa, reshape(view(Y, r1, col), wa.mult[kau], mb), one(T), one(T))
            end
        end
        if q2 != 0
            r2 = _blockrange(Lsrc, q2, ma, wb.mult[kbu])
            FbT = transpose(Mb.F[i, kb])
            for col in 1:m
                Zb = reshape(view(Z, r, col), ma, mb)
                mul!(Zb, reshape(view(Y, r2, col), ma, wb.mult[kbu]), FbT, one(T), one(T))
            end
        end
    end
    return Z
end

"""
    _phi(X::Intertwiner, kc) -> Array (total × mult_c(kc) × N)

Images of the basis of the weight space `kc` of `V_c` under the intertwiners.
"""
function _phi(X::Intertwiner{T}, kc::Int) where {T}
    lock(X.lock) do
        haskey(X.phi, kc) && return X.phi[kc]
        wc = X.Mc.ws
        Ld = _layout(X, kc)
        ups = [i for i in 1:size(wc.up, 1) if wc.up[i, kc] != 0]
        Rinv = _lower_right_inverse(X.Mc, kc)
        out = zeros(T, Ld.total, wc.mult[kc], X.N)
        for k in 1:X.N
            pieces = Matrix{T}[]
            for i in ups
                ku = wc.up[i, kc]
                Yu = _phi(X, ku)
                push!(pieces, _apply_lower(X.Ma, X.Mb, i, _layout(X, ku), view(Yu, :, :, k), Ld))
            end
            out[:, :, k] = reduce(hcat, pieces) * Rinv
        end
        X.phi[kc] = out
        return out
    end
end

# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

const _INTERTWINER_CACHE = Dict{Any, Any}()
const _INTERTWINER_LOCK = ReentrantLock()

const _MEMO = Dict{Any, Any}()
const _MEMO_LOCK = ReentrantLock()

function _clear_intertwiner_cache!()
    lock(() -> empty!(_INTERTWINER_CACHE), _INTERTWINER_LOCK)
    lock(() -> empty!(_MEMO), _MEMO_LOCK)
end

# memoised results (fusion tensors, recoupling symbols); computed outside the lock
function _memo(f, key)
    found = lock(() -> get(_MEMO, key, nothing), _MEMO_LOCK)
    found === nothing || return found
    val = f()
    lock(() -> get!(_MEMO, key, val), _MEMO_LOCK)
end

"""
    intertwiner(a, b, c; basis = :lowering, T = Float64) -> Intertwiner

Basis of `Hom(V_c, V_a ⊗ V_b)` in Kiwi's convention (cached).
"""
function intertwiner(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64)
    key = (a.algebra, a.dynkin_labels, b.dynkin_labels, c.dynkin_labels, basis, T)
    found = lock(() -> get(_INTERTWINER_CACHE, key, nothing), _INTERTWINER_LOCK)
    found === nothing || return found::Intertwiner{T}
    X = Intertwiner(irrep_module(a; basis, T), irrep_module(b; basis, T), irrep_module(c; basis, T))
    lock(_INTERTWINER_LOCK) do
        get!(_INTERTWINER_CACHE, key, X)
    end::Intertwiner{T}
end

# ---------------------------------------------------------------------------
# Clebsch–Gordan tensor
# ---------------------------------------------------------------------------

"""
    fusion_tensor(a, b, c; basis = :lowering, T = Float64, sparse = false)

Clebsch–Gordan coefficients `C[i, j, k, μ] = ⟨a i; b j | c k; μ⟩` in orthonormal bases
(TensorKit's `fusiontensor` layout): the intertwiner `μ` maps the basis vector `k` of
`V_c` to `Σᵢⱼ C[i, j, k, μ] eᵢ ⊗ eⱼ`.  For each `μ` the columns are orthonormal.
Basis vectors are ordered as in [`basis_weights`](@ref); the highest weight of each
irrep comes first.

With `sparse = true` the same coefficients are returned as a `SparseMatrixCSC` of size
`(dim a · dim b) × (dim c · N)`, with row `i + dim(a)·(j − 1)` and column
`k + dim(c)·(μ − 1)`: the isometry `V_c ⊗ ℂᴺ → V_a ⊗ V_b`.  Only entries with
`|C| > 1e-14` are stored.
"""
function fusion_tensor(a::Irrep, b::Irrep, c::Irrep; basis::Symbol = :lowering, T::Type = Float64,
                       sparse::Bool = false)
    _isfloat(T) || error("fusion_tensor needs a floating-point type (orthonormal basis)")
    X = intertwiner(a, b, c; basis, T)
    if sparse
        return _memo((:fusion_tensor_sparse, X)) do
            _fusion_tensor_sparse(X)
        end::SparseMatrixCSC{T, Int}
    end
    n = dimension(X.Ma) * dimension(X.Mb) * dimension(X.Mc) * X.N
    n > 10^7 && return _fusion_tensor(X)          # too large to keep
    return _memo((:fusion_tensor, X)) do
        _fusion_tensor(X)
    end::Array{T, 4}
end

# visit every coefficient: f(i, j, k, μ, value)
function _foreach_coefficient(f, X::Intertwiner{T}) where {T}
    wa, wb, wc = X.Ma.ws, X.Mb.ws, X.Mc.ws
    for kc in 1:length(wc)
        Y = _phi(X, kc)
        L = _layout(X, kc)
        for p in eachindex(L.ka)
            ka, kb = L.ka[p], L.kb[p]
            ma, mb = wa.mult[ka], wb.mult[kb]
            r = _blockrange(L, p, ma, mb)
            for μ in 1:X.N, m in 1:wc.mult[kc]
                blk = reshape(view(Y, r, m, μ), ma, mb)
                for j in 1:mb, i in 1:ma
                    f(wa.offset[ka] + i, wb.offset[kb] + j, wc.offset[kc] + m, μ, blk[i, j])
                end
            end
        end
    end
end

function _fusion_tensor(X::Intertwiner{T}) where {T}
    C = zeros(T, dimension(X.Ma), dimension(X.Mb), dimension(X.Mc), X.N)
    X.N == 0 && return C
    _foreach_coefficient(X) do i, j, k, μ, v
        C[i, j, k, μ] = v
    end
    return C
end

function _fusion_tensor_sparse(X::Intertwiner{T}) where {T}
    da, db, dc = dimension(X.Ma), dimension(X.Mb), dimension(X.Mc)
    rows = Int[]; cols = Int[]; vals = T[]
    X.N == 0 || _foreach_coefficient(X) do i, j, k, μ, v
        if abs(v) > 1e-14
            push!(rows, i + da * (j - 1)); push!(cols, k + dc * (μ - 1)); push!(vals, v)
        end
    end
    return sparse(rows, cols, vals, da * db, dc * X.N)
end
