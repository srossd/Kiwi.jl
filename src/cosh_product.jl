"""
Products of the form ∏_w (e^{w/2} + e^{-w/2}) over the weights of a self-conjugate
representation, decomposed into irreps.

For a real representation r: g → so(N) the product over one weight from each
±pair is the character of the spinor module of so(N) restricted to g (up to the
powers of 2 contributed by the zero weights).  For the adjoint representation it is
Kostant's 2^{rk g} · V_ρ.
"""

"""
    frobenius_schur_indicator(rep::Irrep) -> Int

Return `1` if `rep` is real (orthogonal), `-1` if it is pseudo-real (symplectic)
and `0` if it is not self-conjugate.

For a self-conjugate irrep with highest weight λ the indicator is
`(-1)^⟨λ, 2ρ^∨⟩`, where `2ρ^∨` is the sum of the positive coroots.
"""
function frobenius_schur_indicator(rep::Irrep)
    is_self_conjugate(rep) || return 0
    ad = algebra_data(rep.algebra)
    s = 0
    for a in ad.coroots
        s += dot(a, rep.dynkin_labels)
    end
    return iseven(s) ? 1 : -1
end

"""
    is_real(rep::Irrep) -> Bool

Whether `rep` is a real (orthogonal) representation, i.e. its Frobenius–Schur
indicator is `+1`.
"""
is_real(rep::Irrep) = frobenius_schur_indicator(rep) == 1

"""
    cosh_product(r; half=true, levi=:auto, verbose=false) -> Rep

Decompose the character

    χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2})

into irreps, where `Φ⁺(r)` contains one weight from each pair `±w` of non-zero
weights of the self-conjugate representation `r` (with multiplicity) and `m₀` is
the multiplicity of the zero weight.  Each zero weight contributes a factor
`e^0 + e^0 = 2`.  For a real `r` of dimension `N`, `χ/2^{⌈m₀/2⌉}` is the
character of the spinor module of `so(N)` (the Dirac spinor, i.e. both
chiralities, when `N` is even) restricted to `g`; for the adjoint
representation `χ = 2^{rk g} · χ_ρ` (Kostant).

With `half=false` the product runs over *all* weights of `r`,
`χ = ∏_w (e^{w/2} + e^{-w/2})`, which equals the character of the exterior
algebra `Λ•(r)`.

`r` may be an `Irrep` or a `Rep`; it must be self-conjugate (as a multiset of
weights), and for `half=true` the result must be integral (which holds for any
real representation).

# Algorithm
Write `A_λ = Σ_{w∈W} sign(w) e^{wλ}`.  The multiplicity of the irrep `ν` in `χ` is
the coefficient of `A_{ν+ρ}` in `χ · A_ρ`, so it suffices to compute that product
in the ring of alternants.  All weights are kept in doubled coordinates, so that
the half-weights `w/2` are integral.

1. The non-zero weights of `r` are grouped into Weyl orbits; the product over the
   ±pairs in an orbit (or in an orbit together with its negative) is Weyl
   invariant.
2. A parabolic subgroup `W_J ⊂ W` (a Levi subalgebra) is chosen, and every orbit
   block is split further into `W_J`-invariant pieces, graded by the level
   `⟨w, ξ⟩` along the deleted nodes.  `A_ρ` is rewritten as a signed sum of
   `W_J`-alternants over the cosets `W_J \\ W`.
3. Each piece is expanded as a product of binomials, keeping only `J`-dominant
   weights (partial sums that can no longer reach the chamber are pruned), and is
   multiplied in with the Brauer–Klimyk rule `A^J_z · F = Σ_y F(y) A^J_{z+y}`
   followed by Racah–Speiser reflection by `W_J`.  States that can no longer
   become dominant for the full algebra are discarded.
4. The coefficients of the strictly dominant labels are the multiplicities.

`J = ` all nodes (`levi = :none`) uses the full Weyl group throughout; smaller `J`
gives smaller pieces at the price of more cosets.  With `levi = :auto` the full
algebra and the maximal Levi subalgebras obtained by deleting an end node of the
Dynkin diagram are compared using the sizes of their expanded pieces.  A vector
of node indices may also be passed.  Arithmetic is exact: `Int`, `Int128`, or a
multi-modular representation reconstructed with the Chinese remainder theorem
when multiplicities exceed 128 bits.  The Racah–Speiser steps use all available
threads (`julia -t auto`).

# Example
```julia
cosh_product(Irrep(SU(2), [4]))            # 2 × [3]
cosh_product(adjoint_irrep(G_series(2)))   # 4 × [1,1]
cosh_product(Irrep(SO(16), [0,0,0,0,0,0,0,1]); verbose = true)   # 135 irreps (E8/Spin(16))
```
"""
cosh_product(r::Irrep; kwargs...) = cosh_product(Rep(r); kwargs...)

function cosh_product(r::Rep; half::Bool = true, verbose::Bool = false, levi = :auto)
    g = r.algebra
    ad = algebra_data(g)
    n = ad.n

    # Dominant weights of r with multiplicities
    domm = Dict{Vector{Int}, Int}()
    for (irrep, mult) in r.components
        mult == 0 && continue
        doms, mults = dominant_character(g, irrep.dynkin_labels)
        for (μ, m) in zip(doms, mults)
            domm[μ] = get(domm, μ, 0) + Int(mult) * m
        end
    end
    m0 = get(domm, zeros(Int, n), 0)

    # Group the non-zero weights into blocks of ±pairs closed under W
    blocks = Tuple{Vector{Vector{Int}}, Int}[]   # (pair representatives, copies)
    done = Set{Vector{Int}}()
    for μ in sort!(collect(keys(domm)))
        (all(iszero, μ) || μ in done) && continue
        m = domm[μ]
        μbar = -μ
        _reflect_to_dominant!(μbar, ad)
        get(domm, μbar, 0) == m ||
            error("cosh_product requires a self-conjugate representation")
        orb, _ = _orbit(ad, μ)
        if μbar == μ
            pairs = [w for w in orb if _lex_positive(w)]
        else
            pairs = orb
            push!(done, μbar)
        end
        push!(done, μ)
        push!(blocks, (pairs, half ? m : 2m))
    end

    npairs = sum((length(p) * c for (p, c) in blocks); init = 0)
    E = npairs + m0                              # dim χ = 2^E
    T = _exact_int_type(E + 2)

    # Choose the Levi subalgebra whose Weyl group symmetry is used in the
    # intermediate products
    cache = Dict{Vector{Vector{Int}}, Any}()
    if levi === :auto
        Jl, pieces, cache = _choose_levi(ad, blocks, T)
    else
        Jl = (levi === :none || levi === :full) ? collect(1:n) : sort!(collect(Int, levi))
        pieces = _levi_pieces(ad, blocks, Jl)
    end
    verbose && println("cosh_product: $(length(blocks)) W-orbits, $npairs binomials, dim = 2^$E, ",
                       "arithmetic $T, Levi nodes $Jl, $(length(pieces)) W_J-invariant pieces, ",
                       "$(Threads.nthreads()) threads")

    out = _levi_alternant_product(ad, pieces, Jl, T(big(2)^m0), verbose; cache = cache)

    # Divide by A_ρ: z = 2(ν + ρ)
    res = Dict{Vector{Int}, BigInt}()
    for (z, c) in out
        all(iseven, z) || error("cosh_product: the product is not an integral character " *
                                "(is the representation pseudo-real?)")
        cb = _to_bigint(c)
        cb < 0 && error("cosh_product: negative multiplicity (internal error)")
        res[(z .÷ 2) .- 1] = cb
    end
    return _dict_to_rep(g, res)
end

# One representative of each ±pair: first non-zero coordinate positive
function _lex_positive(w::Vector{Int})
    for x in w
        x != 0 && return x > 0
    end
    return false
end

# ---------------------------------------------------------------------------
# Packed weight keys
# ---------------------------------------------------------------------------

"""
Bit-packing of integer vectors with coordinates in `lo[i]:hi[i]` into a single
unsigned integer of type `U`.  Adding a packed difference vector is a single
(wrapping) integer addition as long as every coordinate stays in range.
"""
struct _Packer{U<:Unsigned}
    lo::Vector{Int}
    shifts::Vector{Int}
    masks::Vector{U}
end

function _packer_type(lo::Vector{Int}, hi::Vector{Int})
    bits = sum(i -> max(1, ndigits(hi[i] - lo[i], base = 2)), eachindex(lo))
    bits <= 64 && return UInt64
    bits <= 128 && return UInt128
    return nothing
end

function _Packer{U}(lo::Vector{Int}, hi::Vector{Int}) where {U}
    shifts = Int[]
    masks = U[]
    sh = 0
    for i in eachindex(lo)
        w = max(1, ndigits(hi[i] - lo[i], base = 2))
        push!(shifts, sh)
        push!(masks, (one(U) << w) - one(U))
        sh += w
    end
    _Packer{U}(copy(lo), shifts, masks)
end

@inline function _pack(p::_Packer{U}, x::AbstractVector{Int}) where {U}
    k = zero(U)
    @inbounds for i in eachindex(p.lo)
        k |= U(x[i] - p.lo[i]) << p.shifts[i]
    end
    return k
end

@inline function _field(p::_Packer{U}, k::U, i::Int) where {U}
    @inbounds return Int((k >> p.shifts[i]) & p.masks[i])
end

@inline function _unpack!(x::AbstractVector{Int}, p::_Packer{U}, k::U) where {U}
    @inbounds for i in eachindex(p.lo)
        x[i] = _field(p, k, i) + p.lo[i]
    end
    return x
end

# Packed representation of a difference vector (wrapping arithmetic)
function _pack_delta(p::_Packer{U}, v::AbstractVector{Int}) where {U}
    d = big(0)
    for i in eachindex(v)
        d += big(v[i]) << p.shifts[i]
    end
    return U(mod(d, big(2)^(8 * sizeof(U))))
end

# ---------------------------------------------------------------------------
# Parabolic (Levi) machinery
# ---------------------------------------------------------------------------

"""
    _reflect_J!(x, ad, Jl) -> sign

Reflect `x` into the dominant chamber of the parabolic subgroup `W_J` generated
by the simple reflections in `Jl`.
"""
@inline function _reflect_J!(x::AbstractVector{Int}, ad::AlgebraData, Jl::Vector{Int})
    n = ad.n
    A = ad.α
    s = 1
    @inbounds while true
        i = 0
        for j in Jl
            if x[j] < 0
                i = j
                break
            end
        end
        i == 0 && return s
        c = x[i]
        for j in 1:n
            x[j] -= c * A[j, i]
        end
        s = -s
    end
end

@inline function _rs_reflect_J!(x::AbstractVector{Int}, ad::AlgebraData, Jl::Vector{Int})
    s = _reflect_J!(x, ad, Jl)
    @inbounds for j in Jl
        x[j] == 0 && return 0
    end
    return s
end

# W_J-orbit of a J-dominant weight (each element once)
function _orbit_J(ad::AlgebraData, μ::Vector{Int}, Jl::Vector{Int})
    length(Jl) == ad.n && return first(_orbit(ad, μ))
    n = ad.n
    A = ad.α
    out = [copy(μ)]
    layer = [copy(μ)]
    while !isempty(layer)
        next = Vector{Vector{Int}}()
        nextset = Set{Vector{Int}}()
        for ν in layer, i in Jl
            c = ν[i]
            c > 0 || continue
            η = ν .- c .* view(A, :, i)
            if !(η in nextset)
                push!(nextset, η)
                push!(next, η)
            end
        end
        append!(out, next)
        layer = next
    end
    return out
end

"""
    _coset_alternants(ad, Jl) -> Vector{Tuple{Vector{Int}, Int}}

The Weyl alternant `A_ρ = Σ_{w ∈ W} sign(w) e^{wρ}` written as a sum of
`W_J`-alternants: the `J`-dominant elements `x` of the orbit `W(2ρ)` (doubled
coordinates) with their signs, one per coset in `W_J \\ W`.
"""
function _coset_alternants(ad::AlgebraData, Jl::Vector{Int})
    n = ad.n
    # Cosets W_J w are explored through right multiplication w ↦ w sᵢ (the Cayley
    # graph of W projects onto the coset graph).  Each coset is represented by the
    # J-dominant point x = w(2ρ) together with the matrix of w in the Dynkin basis.
    start = fill(2, n)
    seen = Dict{Vector{Int}, Int}(start => 1)
    queue = [(start, Matrix{Int}(I, n, n), 1)]
    k = 1
    while k <= length(queue)
        x, Mw, sx = queue[k]
        for i in 1:n
            # w sᵢ (2ρ) = w(2ρ - 2αᵢ) = x - 2 w(αᵢ)
            wαi = Mw * view(ad.α, :, i)
            y = x .- 2 .* wαi
            Mi = copy(Mw)
            Mi .-= wαi * transpose(Matrix{Int}(I, n, n)[i, :])     # Mw · sᵢ
            sy = -sx
            # move to the J-dominant representative of the coset (left W_J action)
            while true
                j = findfirst(jj -> y[jj] < 0, Jl)
                j === nothing && break
                jj = Jl[j]
                c = y[jj]
                y .-= c .* view(ad.α, :, jj)
                # left multiplication by s_jj: M ↦ M - α_jj (row jj of M)
                Mi .-= view(ad.α, :, jj) * transpose(Mi[jj, :])
                sy = -sy
            end
            if !haskey(seen, y)
                seen[y] = sy
                push!(queue, (y, Mi, sy))
            end
        end
        k += 1
    end
    return collect(seen)
end

"""
    _levi_pieces(ad, blocks, Jl)

Split the `W`-invariant blocks of binomials into `W_J`-invariant pieces.  The
weights of a block (`±` the pair representatives) are grouped into `W_J`-orbits
and graded by the level `⟨w, ξ⟩`, `ξ = Σ_{i ∉ J} ωᵢ^∨`, which is `W_J`-invariant.
Orbits of positive level form a piece on their own (one weight of each pair);
level-zero orbits are paired with their negatives (or halved when closed under
negation).
"""
function _levi_pieces(ad::AlgebraData, blocks, Jl::Vector{Int})
    n = ad.n
    length(Jl) == n && return [(p, c) for (p, c) in blocks]
    notJ = setdiff(1:n, Jl)
    Cinv = inv(Rational{Int}.(transpose(ad.C)))        # Dynkin → simple-root coordinates
    ℓ = vec(sum(Cinv[notJ, :], dims = 1))
    level(w) = dot(ℓ, w)
    pieces = Tuple{Vector{Vector{Int}}, Int}[]
    for (pairs, copies) in blocks
        allw = vcat(pairs, [-w for w in pairs])
        seen = Set{Vector{Int}}()
        for w in allw
            w in seen && continue
            x = copy(w)
            _reflect_J!(x, ad, Jl)
            orb = _orbit_J(ad, x, Jl)
            union!(seen, orb)
            lv = level(x)
            if lv > 0
                push!(pieces, (orb, copies))
            elseif lv == 0
                xm = -x
                _reflect_J!(xm, ad, Jl)
                if xm == x
                    push!(pieces, ([v for v in orb if _lex_positive(v)], copies))
                elseif x > xm
                    push!(pieces, (orb, copies))
                end
            end
        end
    end
    return pieces
end

# Order of the stabiliser in W_J of a J-dominant weight μ
function _stabilizer_order_J(ad::AlgebraData, μ::Vector{Int}, Jl::Vector{Int})
    ν = ones(Int, ad.n)
    for j in Jl
        ν[j] = μ[j]
    end
    return _stabilizer_order(ad, ν)
end

# Number of weights obtained by expanding J-dominant weights into W_J-orbits
function _orbit_count_J(ad::AlgebraData, D::Vector{Vector{Int}}, Jl::Vector{Int})
    μJ = ones(Int, ad.n)
    μJ[Jl] .= 0
    WJ = _stabilizer_order(ad, μJ)
    return sum((div(WJ, _stabilizer_order_J(ad, μ, Jl)) for μ in D); init = big(0))
end

"""
    _choose_levi(ad, blocks, T) -> (Jl, pieces, cache)

Choose the Levi subalgebra whose Weyl group symmetry is used in the intermediate
products.  Candidates are the full algebra and the maximal Levi subalgebras
obtained by deleting an end node of the Dynkin diagram (these have the fewest
cosets).  Each candidate's `W_J`-invariant pieces are expanded (with a budget on
the number of states) and the cost is estimated as

    (number of cosets W_J\\W) × (total number of weights of all pieces),

a proxy for the Racah–Speiser work.  If every candidate exceeds the budget, the
maximal Levi with the smallest largest piece is used.  The expansions of the
chosen candidate are returned in `cache` for reuse.
"""
function _choose_levi(ad::AlgebraData, blocks, ::Type{T}; budget::Int = 2_000_000) where {T}
    n = ad.n
    full = collect(1:n)
    maxblock = maximum((length(p) for (p, _) in blocks); init = 0)
    if maxblock <= 16 || n == 1
        return full, _levi_pieces(ad, blocks, full), Dict{Vector{Vector{Int}}, Any}()
    end
    degree(i) = count(j -> j != i && ad.C[i, j] != 0, 1:n)
    candidates = vcat([full], [setdiff(full, k) for k in 1:n if degree(k) <= 1])
    best = nothing
    for Jl in candidates
        pieces = _levi_pieces(ad, blocks, Jl)
        μJ = ones(Int, n)
        μJ[Jl] .= 0
        ncos = div(big(ad.weyl_order), _stabilizer_order(ad, μJ))
        cache = Dict{Vector{Vector{Int}}, Any}()
        total = big(0)
        for (vecs, copies) in sort(pieces, by = p -> -length(p[1]))
            haskey(cache, vecs) && continue
            res = _binomial_dominant(vecs, n, T; coords = Jl, budget = budget)
            if res === nothing
                total = big(-1)
                break
            end
            cache[vecs] = res
            total += copies * _orbit_count_J(ad, res[1], Jl)
        end
        total < 0 && continue
        cost = ncos * total
        if best === nothing || cost < best[1]
            best = (cost, Jl, pieces, cache)
        end
    end
    if best === nothing
        # Everything is large: use the maximal Levi with the smallest pieces
        options = [(maximum(length(p) for (p, _) in _levi_pieces(ad, blocks, setdiff(full, k))), k)
                   for k in 1:n]
        k = minimum(options)[2]
        Jl = setdiff(full, k)
        return Jl, _levi_pieces(ad, blocks, Jl), Dict{Vector{Vector{Int}}, Any}()
    end
    return best[2], best[3], best[4]
end

"""
    _levi_alternant_product(ad, pieces, Jl, c0, verbose)

Compute `c0 · A_ρ · ∏ pieces` in the ring of `W_J`-alternants and return the
coefficients of the `g`-strictly-dominant labels (doubled coordinates), which are
the multiplicities of the irreps in the product.
"""
function _levi_alternant_product(ad::AlgebraData, pieces, Jl::Vector{Int}, c0::T,
                                 verbose::Bool; cache = Dict{Vector{Vector{Int}}, Any}()) where {T}
    n = ad.n
    notJ = setdiff(1:n, Jl)

    # Expand every piece on the J-dominant chamber, then into full W_J-orbits
    tasks = [haskey(cache, vecs) ? cache[vecs] : Threads.@spawn(_binomial_dominant(vecs, n, T; coords = Jl))
             for (vecs, _) in pieces]
    expanded = Tuple{Vector{Vector{Int}}, Vector{T}, Int}[]
    for (i, (vecs, copies)) in enumerate(pieces)
        D, DM = tasks[i] isa Task ? fetch(tasks[i]) : tasks[i]
        W = Vector{Vector{Int}}()
        M = T[]
        for (μ, m) in zip(D, DM)
            orb = _orbit_J(ad, μ, Jl)
            append!(W, orb)
            append!(M, fill(m, length(orb)))
        end
        verbose && println("  piece of $(length(vecs)) binomials: $(length(D)) J-dominant / ",
                           "$(length(W)) weights (× $copies)")
        push!(expanded, (W, M, copies))
    end
    sort!(expanded, by = e -> -length(e[1]))
    steps = Tuple{Vector{Vector{Int}}, Vector{T}}[]
    for (W, M, copies) in expanded, _ in 1:copies
        push!(steps, (W, M))
    end
    S = length(steps)

    # rem[:, s] = largest possible increase of the non-J coordinates after step s
    maxc = [isempty(W) ? zeros(Int, n) : [maximum(y[i] for y in W) for i in 1:n] for (W, _) in steps]
    rem = zeros(Int, n, S + 1)
    for s in S:-1:1
        rem[:, s] = rem[:, s + 1] .+ maxc[s]
    end

    # Starting point: A_ρ as a sum of W_J-alternants
    start = _coset_alternants(ad, Jl)
    verbose && println("  $(length(start)) cosets of W_J in W")

    # Coordinate bounds for packing
    nrm(v) = sqrt(dot(v, ad.G * v) / ad.Gden)
    N = nrm(fill(2, n)) + sum((maximum(nrm, W) for (W, _) in steps if !isempty(W)); init = 0.0)
    lo = ones(Int, n)
    hi = [ceil(Int, 2N / sqrt(ad.norms[i] / ad.Gden)) + 2 for i in 1:n]
    for i in notJ
        lo[i] = 1 - rem[i, 1]
        hi[i] = maximum(x[i] for (x, _) in start) + rem[i, 1]
    end

    U = _packer_type(lo, hi)
    if U === nothing
        X = Dict{Vector{Int}, T}()
        for (x, s) in start
            all(i -> x[i] + rem[i, 1] >= 1, notJ) || continue
            X[x] = s > 0 ? c0 : -c0
        end
        for (si, (W, M)) in enumerate(steps)
            X = _levi_step_generic(ad, Jl, notJ, view(rem, :, si + 1), X, W, M)
            verbose && println("  alternant terms: $(length(X))")
        end
        return X
    end
    pk = _Packer{U}(lo, hi)
    X = Dict{U, T}()
    for (x, s) in start
        all(i -> x[i] + rem[i, 1] >= 1, notJ) || continue
        X[_pack(pk, x)] = s > 0 ? c0 : -c0
    end
    for (si, (W, M)) in enumerate(steps)
        X = _levi_step(ad, Jl, notJ, rem[:, si + 1], pk, X, W, M)
        verbose && println("  alternant terms: $(length(X))")
    end
    out = Dict{Vector{Int}, T}()
    for (k, c) in X
        out[_unpack!(zeros(Int, n), pk, k)] = c
    end
    return out
end

function _levi_step(ad::AlgebraData, Jl::Vector{Int}, notJ::Vector{Int}, remv::Vector{Int},
                    pk::_Packer{U}, X::Dict{U, T}, W::Vector{Vector{Int}}, M::Vector{T}) where {U, T}
    n = ad.n
    entries = collect(X)
    isempty(entries) && return X
    nchunks = min(length(entries), 4 * Threads.nthreads())
    results = Vector{Dict{U, T}}(undef, nchunks)
    Threads.@threads for ci in 1:nchunks
        Yc = Dict{U, T}()
        z = zeros(Int, n)
        buf = zeros(Int, n)
        for e in ci:nchunks:length(entries)
            key, c = entries[e]
            _unpack!(z, pk, key)
            @inbounds for k in eachindex(W)
                y = W[k]
                for j in 1:n
                    buf[j] = z[j] + y[j]
                end
                s = _rs_reflect_J!(buf, ad, Jl)
                s == 0 && continue
                ok = true
                for i in notJ
                    if buf[i] + remv[i] < 1
                        ok = false
                        break
                    end
                end
                ok || continue
                key2 = _pack(pk, buf)
                v = c * M[k]
                idx = Base.ht_keyindex(Yc, key2)
                if idx > 0
                    Yc.vals[idx] = s > 0 ? Yc.vals[idx] + v : Yc.vals[idx] - v
                else
                    Yc[key2] = s > 0 ? v : -v
                end
            end
        end
        results[ci] = Yc
    end
    Y = results[1]
    for ci in 2:nchunks
        for (k, v) in results[ci]
            idx = Base.ht_keyindex(Y, k)
            if idx > 0
                Y.vals[idx] += v
            else
                Y[k] = v
            end
        end
    end
    filter!(p -> !iszero(p.second), Y)
    return Y
end

function _levi_step_generic(ad::AlgebraData, Jl, notJ, remv, X::Dict{Vector{Int}, T},
                            W::Vector{Vector{Int}}, M::Vector{T}) where {T}
    Y = Dict{Vector{Int}, T}()
    for (z, c) in X, k in eachindex(W)
        buf = z .+ W[k]
        s = _rs_reflect_J!(buf, ad, Jl)
        s == 0 && continue
        all(i -> buf[i] + remv[i] >= 1, notJ) || continue
        v = c * M[k]
        Y[buf] = s > 0 ? get(Y, buf, zero(T)) + v : get(Y, buf, zero(T)) - v
    end
    filter!(p -> !iszero(p.second), Y)
    return Y
end

# ---------------------------------------------------------------------------
# Products of binomials, restricted to the dominant chamber
# ---------------------------------------------------------------------------

"""
    _binomial_dominant(vecs, n, T) -> (doms, mults)

Dominant part of ∏_{v ∈ vecs} (e^{v/2} + e^{-v/2}) in doubled coordinates.

The product is expanded one binomial at a time; a partial sum `x` is discarded
as soon as it can no longer reach the dominant chamber, i.e. when
`f(x) + Σ_{remaining v} |f(v)| < 0` for one of the functionals `f` = coordinate
`i` or the sum of coordinates (a separating hyperplane must have a normal in the
dual cone).  Binomials with large coordinates are processed first, which makes
the pruning effective early.
"""
function _binomial_dominant(vecs::Vector{Vector{Int}}, n::Int, ::Type{T};
                            coords::Vector{Int} = collect(1:n), budget::Int = typemax(Int)) where {T}
    vecs = sort(vecs, by = v -> -sum(i -> abs(v[i]), coords; init = 0))
    P = length(vecs)
    F = n + 1
    # functionals: coordinate i (only i ∈ coords are used) and the sum over coords
    proj(v, k) = k <= n ? v[k] : sum(i -> v[i], coords; init = 0)
    R = zeros(Int, F, P + 1)
    for t in P:-1:1, k in 1:F
        R[k, t] = R[k, t + 1] + abs(proj(vecs[t], k))
    end
    S = [R[i, 1] for i in 1:n]
    U = _packer_type(-S, S)
    if U === nothing
        return _binomial_dominant_generic(vecs, n, T, R, coords, budget)
    end
    return _binomial_dominant_packed(_Packer{U}(-S, S), vecs, n, T, R, coords, budget)
end

function _binomial_dominant_packed(pk::_Packer{U}, vecs, n, ::Type{T}, R, coords, budget) where {U, T}
    P = length(vecs)
    S = [-pk.lo[i] for i in 1:n]
    Sc = sum(S[i] for i in coords; init = 0)
    cur = Dict{U, T}(_pack(pk, zeros(Int, n)) => one(T))
    for t in 1:P
        dv = _pack_delta(pk, vecs[t])
        thr = [S[i] - R[i, t + 1] for i in coords]  # field_i ≥ thr_i
        thrsum = Sc - R[n + 1, t + 1]               # Σ_{i ∈ coords} field_i ≥ thrsum
        nxt = Dict{U, T}()
        sizehint!(nxt, length(cur) + length(cur) >> 1)
        for (k, c) in cur
            for k2 in (k + dv, k - dv)
                ok = true
                tot = 0
                @inbounds for (ii, i) in enumerate(coords)
                    f = _field(pk, k2, i)
                    if f < thr[ii]
                        ok = false
                        break
                    end
                    tot += f
                end
                (ok && tot >= thrsum) || continue
                idx = Base.ht_keyindex(nxt, k2)
                if idx > 0
                    @inbounds nxt.vals[idx] += c
                else
                    nxt[k2] = c
                end
            end
        end
        cur = nxt
        length(cur) > budget && return nothing
    end
    doms = Vector{Vector{Int}}()
    mults = T[]
    for (k, c) in cur
        push!(doms, _unpack!(zeros(Int, n), pk, k))
        push!(mults, c)
    end
    return doms, mults
end

function _binomial_dominant_generic(vecs, n, ::Type{T}, R, coords, budget) where {T}
    P = length(vecs)
    cur = Dict{Vector{Int}, T}(zeros(Int, n) => one(T))
    for t in 1:P
        v = vecs[t]
        nxt = Dict{Vector{Int}, T}()
        for (x, c) in cur, s in (1, -1)
            y = x .+ s .* v
            all(i -> y[i] + R[i, t + 1] >= 0, coords) || continue
            sum(i -> y[i], coords; init = 0) + R[n + 1, t + 1] >= 0 || continue
            nxt[y] = get(nxt, y, zero(T)) + c
        end
        cur = nxt
        length(cur) > budget && return nothing
    end
    return collect(keys(cur)), collect(values(cur))
end

