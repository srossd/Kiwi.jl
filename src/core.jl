"""
Fast integer kernels used throughout Kiwi.

Everything in this file works with plain `Vector{Int}` weights written in the
fundamental-weight (Dynkin label) basis, and with an integer-scaled version of
the invariant bilinear form.  All per-algebra data (Cartan matrix, metric,
positive roots, coroots, ...) is computed once and cached in an
[`AlgebraData`](@ref) object, so the hot loops never touch rational arithmetic,
never invert matrices and never allocate `Weight` objects.

The public, `Weight`/`Irrep`-based API is implemented on top of these kernels.
"""

"""
    AlgebraData

Cached, integer-valued data for a simple Lie algebra.

Weights are vectors of Dynkin labels.  The inner product of two weights `x`, `y`
is `dot(x, G*y) // Gden`, with the normalisation used by [`inner_product`](@ref)
(long roots have squared length 2).
"""
struct AlgebraData
    algebra::LieAlgebra
    n::Int
    C::Matrix{Int}                         # Cartan matrix
    α::Matrix{Int}                         # α[:, i] = simple root αᵢ in Dynkin labels (= C[i, :])
    G::Matrix{Int}                         # scaled metric in the fundamental-weight basis
    Gden::Int                              # (x, y) = x'Gy / Gden
    norms::Vector{Int}                     # Gden * (αᵢ, αᵢ)
    pos_roots::Vector{Vector{Int}}         # positive roots, Dynkin labels, sorted by height
    pos_roots_simple::Vector{Vector{Int}}  # same roots, simple-root coordinates
    heights::Vector{Int}
    coroots::Vector{Vector{Int}}           # positive coroots in simple-coroot coordinates
    Gρ::Vector{Int}                        # G * ρ
    ρρ::Int                                # Gden * (ρ, ρ)
    weyl_order::Int
    freud_cache::Dict{UInt64, Vector{Tuple{Vector{Int}, Int, Vector{Int}, Int}}}
    stab_cache::Dict{UInt64, BigInt}
    lock::ReentrantLock
end

const _ALGEBRA_DATA = Dict{LieAlgebra, AlgebraData}()
const _ALGEBRA_DATA_LOCK = ReentrantLock()

"""
    algebra_data(g::LieAlgebra) -> AlgebraData

Return (and cache) the integer data used by the fast kernels.
"""
function algebra_data(g::LieAlgebra)
    lock(_ALGEBRA_DATA_LOCK) do
        get!(() -> _build_algebra_data(g), _ALGEBRA_DATA, g)
    end
end

function _build_algebra_data(g::LieAlgebra)
    C = Matrix{Int}(cartan_matrix(g))
    n = size(C, 1)
    α = Matrix{Int}(transpose(C))

    # Metric (same normalisation as the original rational implementation)
    lengths_sq = simple_root_squared_lengths(g)
    Cq = Rational{Int}.(C)
    target = Cq * Diagonal(lengths_sq) / 2
    Gq = inv(Cq) * target * inv(Matrix(transpose(Cq)))
    Gden = 1
    for x in Gq, y in (denominator(x),)
        Gden = lcm(Gden, y)
    end
    for x in lengths_sq
        Gden = lcm(Gden, denominator(x))
    end
    G = Int.(Gq .* Gden)
    norms = [Int(lengths_sq[i] * Gden) for i in 1:n]

    # Positive roots by the standard string algorithm (simple-root coordinates)
    roots_s = Vector{Vector{Int}}()
    seen = Set{Vector{Int}}()
    for i in 1:n
        e = zeros(Int, n); e[i] = 1
        push!(roots_s, e); push!(seen, e)
    end
    level = copy(roots_s)
    while !isempty(level)
        nextlevel = Vector{Vector{Int}}()
        for β in level
            dl = α * β                     # Dynkin labels of β
            for i in 1:n
                # p = how far down the αᵢ-string through β goes
                p = 0
                γ = copy(β)
                while true
                    γ[i] -= 1
                    (γ[i] >= 0 && γ in seen) || break
                    p += 1
                end
                q = p - dl[i]
                if q > 0
                    δ = copy(β); δ[i] += 1
                    if !(δ in seen)
                        push!(seen, δ); push!(roots_s, δ); push!(nextlevel, δ)
                    end
                end
            end
        end
        level = nextlevel
    end
    sort!(roots_s, by = r -> (sum(r), -reverse(r)))  # by height (stable, simple roots first)
    # keep simple roots in natural order at the front
    heights = [sum(r) for r in roots_s]
    roots_f = [α * r for r in roots_s]
    coroots = Vector{Vector{Int}}()
    for (rs, rf) in zip(roots_s, roots_f)
        nrm = dot(rf, G * rf)              # Gden * (β, β)
        push!(coroots, [div(rs[j] * norms[j], nrm) for j in 1:n])
    end
    ρ = ones(Int, n)
    Gρ = G * ρ
    ρρ = dot(ρ, Gρ)
    AlgebraData(g, n, C, α, G, Gden, norms, roots_f, roots_s, heights, coroots, Gρ, ρρ,
                weyl_group_order(g),
                Dict{UInt64, Vector{Tuple{Vector{Int}, Int, Vector{Int}, Int}}}(),
                Dict{UInt64, BigInt}(), ReentrantLock())
end

# ---------------------------------------------------------------------------
# Basic Weyl group kernels
# ---------------------------------------------------------------------------

"""
    _reflect_to_dominant!(x, ad) -> sign

Reflect the integer weight `x` (Dynkin labels) into the dominant chamber in
place, returning the sign `(-1)^ℓ(w)` of the Weyl group element used.
"""
@inline function _reflect_to_dominant!(x::AbstractVector{Int}, ad::AlgebraData)
    n = ad.n
    A = ad.α
    s = 1
    @inbounds while true
        i = 0
        for j in 1:n
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

"""
    _rs_reflect!(x, ad) -> sign

Racah–Speiser step: reflect `x` (a ρ-shifted weight) to the dominant chamber
and return the sign, or 0 if `x` lies on a wall (i.e. is fixed by a reflection).
"""
@inline function _rs_reflect!(x::AbstractVector{Int}, ad::AlgebraData)
    s = _reflect_to_dominant!(x, ad)
    @inbounds for j in 1:ad.n
        x[j] == 0 && return 0
    end
    return s
end

@inline _is_dominant_vec(x::AbstractVector{Int}) = all(>=(0), x)

"""
    _orbit(ad, μ) -> (weights, parities)

Weyl orbit of the dominant integer weight `μ`.  Each orbit element is generated
exactly once; `parities[k]` is the length (mod 2) of the minimal coset
representative mapping `μ` to `weights[k]`.
"""
function _orbit(ad::AlgebraData, μ::Vector{Int})
    n = ad.n
    A = ad.α
    out = [copy(μ)]
    par = [0]
    layer = [copy(μ)]
    p = 0
    while !isempty(layer)
        p += 1
        next = Vector{Vector{Int}}()
        nextset = Set{Vector{Int}}()
        for ν in layer
            @inbounds for i in 1:n
                c = ν[i]
                c > 0 || continue
                η = similar(ν)
                for j in 1:n
                    η[j] = ν[j] - c * A[j, i]
                end
                if !(η in nextset)
                    push!(nextset, η)
                    push!(next, η)
                end
            end
        end
        append!(out, next)
        append!(par, fill(p & 1, length(next)))
        layer = next
    end
    return out, par
end

"""
    _stabilizer_order(ad, μ) -> BigInt

Order of the stabiliser of the dominant weight `μ` (the parabolic subgroup
generated by the simple reflections with `μᵢ = 0`), computed with Macdonald's
formula |W_J| = ∏_{α ∈ Φ_J⁺} (ht α + 1)/ht α.
"""
function _stabilizer_order(ad::AlgebraData, μ::AbstractVector{Int})
    J = _zero_mask(μ)
    lock(ad.lock) do
    get!(ad.stab_cache, J) do
        num = big(1); den = big(1)
        for (rs, h) in zip(ad.pos_roots_simple, ad.heights)
            if _support_in(rs, J)
                num *= h + 1
                den *= h
            end
        end
        div(num, den)
    end
    end
end

_orbit_size(ad::AlgebraData, μ::AbstractVector{Int}) = div(big(ad.weyl_order), _stabilizer_order(ad, μ))

@inline function _zero_mask(μ::AbstractVector{Int})
    J = UInt64(0)
    @inbounds for i in eachindex(μ)
        μ[i] == 0 && (J |= UInt64(1) << (i - 1))
    end
    return J
end

@inline function _support_in(rs::Vector{Int}, J::UInt64)
    @inbounds for i in eachindex(rs)
        rs[i] != 0 && (J >> (i - 1)) & 1 == 0 && return false
    end
    return true
end

# ---------------------------------------------------------------------------
# Freudenthal's formula on dominant weights (Moody–Patera orbit trick)
# ---------------------------------------------------------------------------

"""
Roots needed by Freudenthal's formula at a dominant weight with zero set `J`,
grouped into orbits of the stabiliser `W_J`.

Returns `(β, coeff, Gβ, (β,β)·Gden)` tuples such that

    2 Σ_{α>0} Σ_{k≥1} (μ+kα, α) m(μ+kα) = Σ_{(β,c)} c Σ_{k≥1} (μ+kβ, β) m(μ+kβ).

`W_J` permutes `Φ⁺ ∖ Φ_J⁺` (orbit weight `2|O|`), and `Φ_J` (orbit weight `|O|`,
using the symmetry of root strings through μ for roots orthogonal to μ).
"""
function _freudenthal_roots(ad::AlgebraData, J::UInt64)
    lock(ad.lock) do
    get!(ad.freud_cache, J) do
        n = ad.n
        A = ad.α
        Jidx = [i for i in 1:n if (J >> (i - 1)) & 1 == 1]
        out = Tuple{Vector{Int}, Int, Vector{Int}, Int}[]
        covered = Set{Vector{Int}}()
        function orbitJ(β)
            orb = Set{Vector{Int}}([β])
            stack = [β]
            while !isempty(stack)
                x = pop!(stack)
                for i in Jidx
                    c = x[i]
                    c == 0 && continue
                    y = x .- c .* A[:, i]
                    if !(y in orb)
                        push!(orb, y); push!(stack, y)
                    end
                end
            end
            orb
        end
        for (β, βs) in zip(ad.pos_roots, ad.pos_roots_simple)
            β in covered && continue
            orb = orbitJ(β)
            union!(covered, orb)
            inJ = _support_in(βs, J)
            c = inJ ? length(orb) : 2 * length(orb)
            Gβ = ad.G * β
            push!(out, (β, c, Gβ, dot(β, Gβ)))
        end
        out
    end
    end
end

"""
    _dominant_character(ad, λ) -> (doms, mults)

Dominant weights of the irrep with highest weight `λ` (sorted by depth below
`λ`) and their multiplicities, computed with Freudenthal's formula using only
dominant weights.
"""
function _dominant_character(ad::AlgebraData, λ::Vector{Int})
    n = ad.n
    # 1. dominant weights: every dominant μ ≤ λ is reachable from λ through
    #    dominant weights by subtracting positive roots (Stembridge).
    doms = [copy(λ)]
    depth = [0]
    index = Dict{Vector{Int}, Int}(doms[1] => 1)
    k = 1
    while k <= length(doms)
        μ = doms[k]
        d = depth[k]
        for (β, h) in zip(ad.pos_roots, ad.heights)
            ν = μ .- β
            _is_dominant_vec(ν) || continue
            haskey(index, ν) && continue
            push!(doms, ν); push!(depth, d + h)
            index[ν] = length(doms)
        end
        k += 1
    end
    perm = sortperm(depth)          # stable: higher weights first
    doms = doms[perm]
    for (i, μ) in enumerate(doms)
        index[μ] = i
    end
    mults = zeros(Int, length(doms))
    mults[1] = 1

    # 2. Freudenthal recursion
    Gλρ = ad.G * (λ .+ 1)
    λρ2 = dot(λ .+ 1, Gλρ)
    buf = zeros(Int, n)
    Gμ = zeros(Int, n)
    @inbounds for i in 2:length(doms)
        μ = doms[i]
        mul!(Gμ, ad.G, μ)
        μμ = dot(μ, Gμ)
        μρ = dot(μ, ad.Gρ)
        lhs = λρ2 - (μμ + 2μρ + ad.ρρ)
        rhs = Int128(0)
        for (β, c, _, ββ) in _freudenthal_roots(ad, _zero_mask(μ))
            μβ = dot(Gμ, β)
            kk = 1
            while true
                for j in 1:n
                    buf[j] = μ[j] + kk * β[j]
                end
                _reflect_to_dominant!(buf, ad)
                idx = get(index, buf, 0)
                (idx == 0 || idx >= i) && break
                m = mults[idx]
                rhs += Int128(c) * Int128(μβ + kk * ββ) * Int128(m)
                kk += 1
            end
        end
        q, r = divrem(rhs, Int128(lhs))
        r == 0 || error("Freudenthal recursion produced a non-integer multiplicity at $μ")
        mults[i] = Int(q)
    end
    return doms, mults
end

const _DOMCHAR_CACHE = Dict{Tuple{LieAlgebra, Vector{Int}}, Tuple{Vector{Vector{Int}}, Vector{Int}}}()
const _DOMCHAR_LOCK = ReentrantLock()

"""
    dominant_character(g::LieAlgebra, λ::Vector{Int}) -> (doms, mults)

Cached dominant weights and multiplicities of the irrep with Dynkin labels `λ`.
"""
function dominant_character(g::LieAlgebra, λ::Vector{Int})
    key = (g, λ)
    r = lock(_DOMCHAR_LOCK) do
        get(_DOMCHAR_CACHE, key, nothing)
    end
    r === nothing || return r
    r = _dominant_character(algebra_data(g), copy(λ))
    lock(_DOMCHAR_LOCK) do
        _DOMCHAR_CACHE[(g, copy(λ))] = r
    end
    return r
end

"""
    _full_character(g, λ) -> (weights, mults)

All weights (with multiplicity) of the irrep with Dynkin labels `λ`.
"""
function _full_character(g::LieAlgebra, λ::Vector{Int})
    ad = algebra_data(g)
    doms, mults = dominant_character(g, λ)
    W = Vector{Vector{Int}}()
    M = Vector{Int}()
    for (μ, m) in zip(doms, mults)
        orb, _ = _orbit(ad, μ)
        append!(W, orb)
        append!(M, fill(m, length(orb)))
    end
    return W, M
end

# ---------------------------------------------------------------------------
# Dimension
# ---------------------------------------------------------------------------

"""
    _dimension(ad, λ) -> BigInt

Weyl dimension formula ∏_{α>0} (λ+ρ, α^∨)/(ρ, α^∨) in exact integer arithmetic.
"""
function _dimension(ad::AlgebraData, λ::AbstractVector{Int})
    num = Int128(1); den = Int128(1)
    big_mode = false
    bnum = big(1); bden = big(1)
    for a in ad.coroots
        p = Int128(0); q = Int128(0)
        @inbounds for j in eachindex(a)
            p += a[j] * (λ[j] + 1)
            q += a[j]
        end
        if !big_mode
            g1 = gcd(p, den); p1 = div(p, g1); den1 = div(den, g1)
            g2 = gcd(q, num); q1 = div(q, g2); num1 = div(num, g2)
            n2, o1 = Base.mul_with_overflow(num1, p1)
            d2, o2 = Base.mul_with_overflow(den1, q1)
            if o1 || o2
                big_mode = true
                bnum = big(num) * p
                bden = big(den) * q
            else
                num, den = n2, d2
            end
        else
            bnum *= p; bden *= q
        end
    end
    return big_mode ? div(bnum, bden) : big(div(num, den))
end

# ---------------------------------------------------------------------------
# Racah–Speiser / Klimyk products
# ---------------------------------------------------------------------------

"""
    _klimyk!(out, ad, base, coeff, weights, mults; scale=1, shift=1)

Accumulate `coeff * Σ_y mults[y] * A_{base + scale*y}` into `out`, where
`A_z` is the Weyl alternant: `base + scale*y` is reflected to the dominant
chamber with its sign, wall weights are dropped and the result is stored keyed
by `z - shift` (i.e. by the highest weight when `base = λ + shift·ρ`).
"""
function _klimyk!(out::Dict{Vector{Int}, T}, ad::AlgebraData, base::Vector{Int}, coeff,
                  weights::Vector{Vector{Int}}, mults::AbstractVector;
                  scale::Int = 1, shift::Int = 1) where {T}
    n = ad.n
    buf = zeros(Int, n)
    @inbounds for k in eachindex(weights)
        y = weights[k]
        for j in 1:n
            buf[j] = base[j] + scale * y[j]
        end
        s = _rs_reflect!(buf, ad)
        s == 0 && continue
        for j in 1:n
            buf[j] -= shift
        end
        v = T(coeff) * T(mults[k])
        idx = Base.ht_keyindex(out, buf)
        if idx > 0
            out.vals[idx] = s > 0 ? out.vals[idx] + v : out.vals[idx] - v
        else
            out[copy(buf)] = s > 0 ? v : -v
        end
    end
    return out
end

"""
    _virtual_times_character(ad, X, weights, mults; scale=1) -> Dict

Multiply the virtual representation `X` (highest weight ⇒ coefficient) by the
character with the given weights (each scaled by `scale`, i.e. an Adams operation
when `scale > 1`) and decompose the result into irreps.
"""
function _virtual_times_character(ad::AlgebraData, X::Dict{Vector{Int}, T},
                                  weights::Vector{Vector{Int}}, mults::AbstractVector;
                                  scale::Int = 1) where {T}
    out = Dict{Vector{Int}, T}()
    for (ν, c) in X
        c == 0 && continue
        _klimyk!(out, ad, ν .+ 1, c, weights, mults; scale = scale)
    end
    filter!(p -> p.second != 0, out)
    return out
end

# Merge a list of (possibly repeated) weights into unique weights with multiplicities
function _merge_weights(W::Vector{Vector{Int}}, M::AbstractVector{T}) where {T}
    d = Dict{Vector{Int}, T}()
    for (w, m) in zip(W, M)
        d[w] = get(d, w, zero(T)) + m
    end
    filter!(p -> p.second != 0, d)
    return collect(keys(d)), collect(values(d))
end

"""
Choose an integer type wide enough for quantities bounded by `bound`.
"""
function _int_type_for(bound::Integer)
    b = big(bound)
    b < big(2)^61 && return Int
    b < big(2)^125 && return Int128
    return BigInt
end
