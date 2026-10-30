# Explicit irreducible modules: weight-space bases and the matrices of the Chevalley
# generators Eᵢ, Fᵢ (simple roots only) between neighbouring weight spaces.
#
# Conventions (the "Kiwi basis", `basis = :lowering`)
# ---------------------------------------------------
# * Weights of V(λ) are ordered by depth (the height of λ − μ), and within a depth
#   by decreasing Dynkin labels (lexicographically).  The highest weight comes first.
# * The contravariant (Shapovalov) form ⟨·,·⟩ is normalised by ⟨v_λ, v_λ⟩ = 1 and
#   satisfies ⟨Fᵢx, y⟩ = ⟨x, Eᵢy⟩; it is positive definite (compact real form).
# * The candidates spanning V(λ)_μ are the vectors Fᵢ b, for i = 1, …, n in turn and
#   b running over the basis of V(λ)_{μ+αᵢ} in order.  The basis of V(λ)_μ is the
#   first independent subfamily of the candidates, taken greedily in this order.
# * Exact arithmetic keeps these vectors (Gram matrix = Shapovalov form).  Over Q the
#   choice is decided exactly; over a prime field Fp{P} it is decided modulo the
#   reference prime 2⁶¹ − 1 and reused for every other P, so all primes reduce the same
#   rational basis.  The
#   orthonormal basis is their Gram–Schmidt orthonormalisation in the same order,
#   so that each new vector has positive overlap with its candidate.  With that
#   normalisation Fᵢ has nonnegative entries along the candidate it came from, and
#   for su(2) it reproduces the Condon–Shortley basis |j m⟩ with F = J₋.
#
# In an orthonormal basis Eᵢ = Fᵢᵀ.  In the exact basis Eᵢ = G⁻¹ Fᵢᵀ G.

"""
    WeightStructure

Weights of an irrep with multiplicities, in Kiwi's canonical order, and the
neighbour tables `up[i, k]` / `down[i, k]` giving the index of `wₖ ± αᵢ` (or 0).
"""
struct WeightStructure
    irrep::Irrep
    weights::Vector{Vector{Int}}
    depth::Vector{Int}
    mult::Vector{Int}
    offset::Vector{Int}          # global basis index of the first vector of each weight space, minus 1
    index::Dict{Vector{Int}, Int}
    up::Matrix{Int}
    down::Matrix{Int}
end

Base.length(ws::WeightStructure) = length(ws.weights)

function WeightStructure(rep::Irrep)
    g = rep.algebra
    ad = algebra_data(g)
    n = ad.n
    λ = rep.dynkin_labels
    W, M = _full_character(g, λ)
    mult = Dict{Vector{Int}, Int}(zip(W, M))
    # depth by lowering from the top (every weight is reached by simple lowerings)
    depth = Dict{Vector{Int}, Int}(λ => 0)
    layer = [λ]
    d = 0
    while !isempty(layer)
        next = Vector{Vector{Int}}()
        d += 1
        for w in layer, i in 1:n
            v = w .- view(ad.α, :, i)
            if haskey(mult, v) && !haskey(depth, v)
                depth[v] = d
                push!(next, v)
            end
        end
        layer = next
    end
    length(depth) == length(mult) || error("weight system of $rep is not connected by lowering")
    weights = collect(keys(depth))
    sort!(weights, by = w -> (depth[w], map(-, w)))
    nw = length(weights)
    index = Dict(w => k for (k, w) in enumerate(weights))
    up = zeros(Int, n, nw); down = zeros(Int, n, nw)
    for (k, w) in enumerate(weights), i in 1:n
        up[i, k] = get(index, w .+ view(ad.α, :, i), 0)
        down[i, k] = get(index, w .- view(ad.α, :, i), 0)
    end
    m = [mult[w] for w in weights]
    off = cumsum(vcat(0, m))[1:nw]
    return WeightStructure(rep, weights, [depth[w] for w in weights], m, off, index, up, down)
end

"""
    IrrepModule{T}

An irreducible module with explicit matrices of the simple-root generators.

`F[i, k]` maps the weight space `k + αᵢ` (index `ws.up[i, k]`) to the weight space
`k`; `E[i, k]` maps weight space `k` to `k + αᵢ`.  Missing neighbours give 0×0
placeholders.  For an orthonormal basis `gram` is empty; otherwise `gram[k]` is the
Gram matrix of the basis of weight space `k`.
"""
struct IrrepModule{T}
    ws::WeightStructure
    basis::Symbol
    F::Matrix{Matrix{T}}
    E::Matrix{Matrix{T}}
    gram::Vector{Matrix{T}}
    provenance::Vector{Vector{Tuple{Int, Int}}}   # (i, b): basis vector ∝ Fᵢ b (lowering basis)
    cache::Dict{Symbol, Any}
    lock::ReentrantLock
end

isorthonormal(M::IrrepModule) = isempty(M.gram)
Base.eltype(::IrrepModule{T}) where {T} = T
dimension(M::IrrepModule) = sum(M.ws.mult)

function Base.show(io::IO, M::IrrepModule{T}) where {T}
    print(io, "IrrepModule{$T}(", M.ws.irrep, ", basis = :", M.basis,
          isorthonormal(M) ? ", orthonormal" : ", Shapovalov", ")")
end

# ---------------------------------------------------------------------------
# Construction in the lowering basis
# ---------------------------------------------------------------------------

"""
    _lowering_module(rep, T; orthonormal, selection) -> IrrepModule{T}

Build V(rep) weight space by weight space.  For candidate vectors Fᵢu, Fⱼv
(u ∈ V_{μ+αᵢ}, v ∈ V_{μ+αⱼ}) the contravariant form is

    ⟨Fᵢu, Fⱼv⟩ = ⟨u, EᵢFⱼv⟩ = ⟨u, FⱼEᵢv⟩ + δᵢⱼ ⟨μ+αᵢ, αᵢ^∨⟩ ⟨u, v⟩,

which only involves weight spaces built earlier.  `selection[k]`, if given, fixes
which candidates form the basis of weight space `k` (it must be the greedy choice;
it is passed in when it was decided in exact arithmetic).
"""
function _lowering_module(rep::Irrep, ::Type{T}; orthonormal::Bool = _isfloat(T),
                          selection = nothing) where {T}
    orthonormal && !_isfloat(T) && error("an orthonormal basis needs a floating-point type")
    ws = WeightStructure(rep)
    n = size(ws.up, 1)
    nw = length(ws)
    empty = zeros(T, 0, 0)
    F = fill(empty, n, nw)
    E = fill(empty, n, nw)
    gram = orthonormal ? Matrix{T}[] : Vector{Matrix{T}}(undef, nw)
    orthonormal || (gram[1] = ones(T, 1, 1))
    prov = Vector{Vector{Tuple{Int, Int}}}(undef, nw)
    prov[1] = Tuple{Int, Int}[]
    for k in 2:nw
        w = ws.weights[k]
        ups = [i for i in 1:n if ws.up[i, k] != 0]
        cands = Tuple{Int, Int}[]
        start = Dict{Int, Int}()
        for i in ups
            start[i] = length(cands)
            for u in 1:ws.mult[ws.up[i, k]]
                push!(cands, (i, u))
            end
        end
        nc = length(cands)
        # K[i, j] = matrix of Eᵢ Fⱼ : V_{μ+αⱼ} → V_{μ+αᵢ}
        K = Dict{Tuple{Int, Int}, Matrix{T}}()
        G = zeros(T, nc, nc)
        for i in ups, j in ups
            ki, kj = ws.up[i, k], ws.up[j, k]
            Kij = zeros(T, ws.mult[ki], ws.mult[kj])
            top = ws.up[j, ki]                       # μ + αᵢ + αⱼ
            if top != 0
                Kij .+= F[j, ki] * E[i, kj]
            end
            if i == j
                h = ws.weights[ki][i]                 # ⟨μ + αᵢ, αᵢ^∨⟩
                for u in 1:ws.mult[ki]
                    Kij[u, u] += h
                end
            end
            K[(i, j)] = Kij
            Gij = orthonormal ? Kij : gram[ki] * Kij
            r0, c0 = start[i], start[j]
            G[(r0 + 1):(r0 + size(Gij, 1)), (c0 + 1):(c0 + size(Gij, 2))] = Gij
        end
        sel = selection === nothing ? _greedy_select(G) : selection[k]
        length(sel) == ws.mult[k] ||
            throw(_SelectionError(rep, w, length(sel), ws.mult[k]))
        GS = G[sel, sel]
        if orthonormal
            R = cholesky(Symmetric(GS)).U
            Y = transpose(R) \ G[sel, :]               # all candidates in the new orthonormal basis
            for i in ups
                cols = (start[i] + 1):(start[i] + ws.mult[ws.up[i, k]])
                F[i, k] = Y[:, cols]
                E[i, k] = Matrix(transpose(F[i, k]))
            end
        else
            X = _solve(GS, G[sel, :])
            gram[k] = GS
            for i in ups
                cols = (start[i] + 1):(start[i] + ws.mult[ws.up[i, k]])
                F[i, k] = X[:, cols]
                Ei = zeros(T, ws.mult[ws.up[i, k]], ws.mult[k])
                for (m, s) in enumerate(sel)
                    (j, v) = cands[s]
                    Ei[:, m] = K[(i, j)][:, v]
                end
                E[i, k] = Ei
            end
        end
        prov[k] = cands[sel]
    end
    return IrrepModule{T}(ws, :lowering, F, E, gram, prov, Dict{Symbol, Any}(), ReentrantLock())
end

struct _SelectionError <: Exception
    rep::Irrep
    weight::Vector{Int}
    found::Int
    expected::Int
end
Base.showerror(io::IO, e::_SelectionError) =
    print(io, "basis selection for $(e.rep) at weight $(e.weight) found $(e.found) independent ",
          "vectors, expected $(e.expected)")

"""
    _exact_selection(rep) -> Vector{Vector{Int}}

The greedy basis choice of every weight space, decided in exact arithmetic modulo a
large prime (a false dependence would show up as a wrong multiplicity, which is
checked against Freudenthal's formula).
"""
function _exact_selection(rep::Irrep)
    M = irrep_module(rep; T = Fp{_P61})
    ws = M.ws
    out = Vector{Vector{Int}}(undef, length(ws))
    out[1] = Int[]
    for k in 2:length(ws)
        # recover candidate positions from provenance
        n = size(ws.up, 1)
        pos = Dict{Tuple{Int, Int}, Int}()
        c = 0
        for i in 1:n
            ws.up[i, k] == 0 && continue
            for u in 1:ws.mult[ws.up[i, k]]
                c += 1
                pos[(i, u)] = c
            end
        end
        out[k] = [pos[p] for p in M.provenance[k]]
    end
    return out
end

# ---------------------------------------------------------------------------
# Cache and public constructors
# ---------------------------------------------------------------------------

const _MODULE_CACHE = Dict{Tuple{LieAlgebra, Vector{Int}, Symbol, DataType}, Any}()
const _MODULE_LOCK = ReentrantLock()

"""
    irrep_module(rep::Irrep; basis = :lowering, T = Float64) -> IrrepModule

Explicit module for `rep` with matrices of the simple-root generators.

* `basis = :lowering` — Kiwi's basis (see the conventions at the top of
  `recoupling/modules.jl`), for every simple Lie algebra.
* `basis = :gelfand_tsetlin` — the Gelfand–Tsetlin basis (su(N) only).

`T` is the scalar type: `Float64` or `BigFloat` give an orthonormal basis;
`Rational{BigInt}` or `Fp{P}` give the exact, unnormalised lowering basis (with its
Gram matrices).
"""
function irrep_module(rep::Irrep; basis::Symbol = :lowering, T::Type = Float64)
    key = (rep.algebra, rep.dynkin_labels, basis, T)
    lock(_MODULE_LOCK) do
        get!(_MODULE_CACHE, key) do
            _build_module(rep, basis, T)
        end
    end::IrrepModule{T}
end

function _build_module(rep::Irrep, basis::Symbol, ::Type{T}) where {T}
    if basis == :lowering
        if _isfloat(T)
            try
                return _lowering_module(rep, T)
            catch err
                err isa _SelectionError || rethrow()
                return _lowering_module(rep, T; selection = _exact_selection(rep))
            end
        end
        # Over other primes, reuse the selection decided modulo the reference prime 2⁶¹ − 1,
        # so that every prime computes residues of the same canonical (rational) basis; a
        # selection that is singular modulo P then fails in the exact solve.
        if T <: Fp && T !== Fp{_P61}
            return _lowering_module(rep, T; orthonormal = false, selection = _exact_selection(rep))
        end
        return _lowering_module(rep, T; orthonormal = false)
    elseif basis == :gelfand_tsetlin
        return _gt_module(rep, T)
    else
        throw(ArgumentError("unknown basis $basis (use :lowering or :gelfand_tsetlin)"))
    end
end

"""
    clear_recoupling_cache!()

Empty the caches of modules, intertwiners and recoupling symbols.
"""
function clear_recoupling_cache!()
    lock(_MODULE_LOCK) do
        empty!(_MODULE_CACHE)
    end
    _clear_intertwiner_cache!()
    return nothing
end

# ---------------------------------------------------------------------------
# Generators as full matrices
# ---------------------------------------------------------------------------

"""
    generators(rep::Irrep; basis = :lowering, T = Float64) -> (E, F, H)

Matrices of the Chevalley generators `Eᵢ`, `Fᵢ`, `Hᵢ` (i = 1, …, rank) in the basis of
[`irrep_module`](@ref), as vectors of sparse matrices.  They satisfy
`[Eᵢ, Fⱼ] = δᵢⱼ Hᵢ`, `[Hᵢ, Eⱼ] = Cⱼᵢ Eⱼ` and the Serre relations; in an orthonormal
basis `Eᵢ = Fᵢᵀ`.  Basis vectors are ordered by weight (see [`basis_weights`](@ref)).
"""
function generators(rep::Irrep; basis::Symbol = :lowering, T::Type = Float64)
    M = irrep_module(rep; basis, T)
    ws = M.ws
    n = size(ws.up, 1)
    d = dimension(M)
    Es = SparseMatrixCSC{T, Int}[]
    Fs = SparseMatrixCSC{T, Int}[]
    Hs = SparseMatrixCSC{T, Int}[]
    for i in 1:n
        Ei = spzeros(T, d, d); Fi = spzeros(T, d, d); Hi = spzeros(T, d, d)
        for k in 1:length(ws)
            r = (ws.offset[k] + 1):(ws.offset[k] + ws.mult[k])
            for a in r
                Hi[a, a] = ws.weights[k][i]
            end
            ku = ws.up[i, k]
            ku == 0 && continue
            ru = (ws.offset[ku] + 1):(ws.offset[ku] + ws.mult[ku])
            Fi[r, ru] = M.F[i, k]
            Ei[ru, r] = M.E[i, k]
        end
        push!(Es, Ei); push!(Fs, Fi); push!(Hs, Hi)
    end
    return Es, Fs, Hs
end

"""
    basis_weights(rep::Irrep; basis = :lowering) -> Vector{Vector{Int}}

The weight (Dynkin labels) of each basis vector of [`irrep_module`](@ref), in order.
"""
function basis_weights(rep::Irrep; basis::Symbol = :lowering)
    ws = basis == :lowering ? WeightStructure(rep) : irrep_module(rep; basis).ws
    out = Vector{Vector{Int}}()
    for (k, w) in enumerate(ws.weights), _ in 1:ws.mult[k]
        push!(out, w)
    end
    return out
end
