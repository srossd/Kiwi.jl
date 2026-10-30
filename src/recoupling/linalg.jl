# Arithmetic and small dense linear algebra for the recoupling code.
#
# All constructions in `recoupling/` are written once, generically over the scalar
# type `T`.  Three kinds of `T` are used:
#
#   * `Float64` / `BigFloat`  — orthonormal bases, the default numerical path;
#   * `Rational{BigInt}`      — exact arithmetic in an unnormalised basis;
#   * `Fp{P}`                 — arithmetic modulo a prime `P`, used to decide linear
#                               independence exactly and for the multi-modular path.
#
# The routines below dispatch on whether `T` is a floating-point type (use LAPACK,
# tolerances) or an exact field (fraction-free decisions, nonzero pivots).

# ---------------------------------------------------------------------------
# Prime field
# ---------------------------------------------------------------------------

"""
    Fp{P} <: Number

Element of the prime field ℤ/Pℤ, `P < 2^62`, stored as an `Int64` in `0:P-1`.
"""
struct Fp{P} <: Number
    v::Int64
    Fp{P}(v::Int64, ::Val{:raw}) where {P} = new{P}(v)
end

Fp{P}(x::Integer) where {P} = Fp{P}(Int64(mod(x, P)), Val(:raw))
Fp{P}(x::Rational) where {P} = Fp{P}(numerator(x)) / Fp{P}(denominator(x))
Fp{P}(x::Fp{P}) where {P} = x

Base.zero(::Type{Fp{P}}) where {P} = Fp{P}(Int64(0), Val(:raw))
Base.one(::Type{Fp{P}}) where {P} = Fp{P}(Int64(1), Val(:raw))
Base.zero(::Fp{P}) where {P} = zero(Fp{P})
Base.one(::Fp{P}) where {P} = one(Fp{P})
Base.iszero(a::Fp) = a.v == 0
Base.:(==)(a::Fp{P}, b::Fp{P}) where {P} = a.v == b.v
Base.hash(a::Fp{P}, h::UInt) where {P} = hash(a.v, hash(P, h))

function Base.:+(a::Fp{P}, b::Fp{P}) where {P}
    s = a.v + b.v
    return Fp{P}(s >= P ? s - P : s, Val(:raw))
end
function Base.:-(a::Fp{P}, b::Fp{P}) where {P}
    s = a.v - b.v
    return Fp{P}(s < 0 ? s + P : s, Val(:raw))
end
Base.:-(a::Fp{P}) where {P} = Fp{P}(a.v == 0 ? Int64(0) : P - a.v, Val(:raw))
Base.:*(a::Fp{P}, b::Fp{P}) where {P} = Fp{P}(Int64(mod(Int128(a.v) * b.v, P)), Val(:raw))
function Base.inv(a::Fp{P}) where {P}
    iszero(a) && throw(DivideError())
    r = Int128(1); x = Int128(a.v); e = P - 2
    while e > 0
        isodd(e) && (r = mod(r * x, P))
        x = mod(x * x, P)
        e >>= 1
    end
    return Fp{P}(Int64(r), Val(:raw))
end
Base.:/(a::Fp{P}, b::Fp{P}) where {P} = a * inv(b)
Base.promote_rule(::Type{Fp{P}}, ::Type{<:Integer}) where {P} = Fp{P}
Base.convert(::Type{Fp{P}}, x::Integer) where {P} = Fp{P}(x)
Base.show(io::IO, a::Fp{P}) where {P} = print(io, a.v, " mod ", P)

# Deterministic Miller–Rabin for 64-bit integers
function _isprime64(n::Int64)
    n < 2 && return false
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        n % p == 0 && return n == p
    end
    d = n - 1; s = 0
    while iseven(d); d >>= 1; s += 1; end
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        x = powermod(Int128(a), d, Int128(n))
        (x == 1 || x == n - 1) && continue
        composite = true
        for _ in 1:(s - 1)
            x = mod(x * x, n)
            if x == n - 1
                composite = false
                break
            end
        end
        composite && return false
    end
    return true
end

"""
    _primes_below(bound, k) -> Vector{Int64}

The `k` largest primes below `bound` (used for the multi-modular path).
"""
function _primes_below(bound::Int64, k::Int)
    out = Int64[]
    n = isodd(bound) ? bound - 2 : bound - 1
    while length(out) < k
        _isprime64(n) && push!(out, n)
        n -= 2
    end
    return out
end

const _P61 = Int64(2)^61 - 1          # Mersenne prime, used for exact selection

# ---------------------------------------------------------------------------
# Field-generic helpers
# ---------------------------------------------------------------------------

_isfloat(::Type{T}) where {T} = T <: AbstractFloat

"""
    _is_negligible(x, ref, T)

Exact types: `x == 0`.  Floating types: `|x| ≤ tol · max(|ref|, 1)`.
"""
@inline _is_negligible(x, ref, ::Type{T}) where {T} =
    _isfloat(T) ? abs(x) <= _tol(T) * max(abs(ref), one(T)) : iszero(x)

_tol(::Type{Float64}) = 1e-9
_tol(::Type{BigFloat}) = BigFloat(10)^(-precision(BigFloat) ÷ 5)
_tol(::Type{T}) where {T} = zero(T)

"""
    _greedy_select(G, T) -> Vector{Int}

Indices of a maximal linearly independent prefix-greedy subset of a family of
vectors with Gram matrix `G`: vector `j` is taken iff it is independent of the
vectors taken before it.  Uses symmetric Gaussian elimination (Schur complements).
"""
function _greedy_select(G::AbstractMatrix{T}) where {T}
    n = size(G, 1)
    W = Matrix{T}(G)
    sel = Int[]
    for j in 1:n
        d = W[j, j]
        if _is_negligible(d, G[j, j], T)
            continue
        end
        push!(sel, j)
        # eliminate column j from the later rows/columns
        for k in (j + 1):n
            f = W[k, j] / d
            iszero(f) && continue
            for l in (j + 1):n
                W[k, l] -= f * W[j, l]
            end
        end
    end
    return sel
end

"""
    _solve(A, B) -> X

Solve `A X = B` for square nonsingular `A`.
"""
_solve(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T <: AbstractFloat} = A \ B
function _solve(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    M = hcat(Matrix{T}(A), Matrix{T}(B))
    m = size(M, 2)
    for c in 1:n
        p = findfirst(r -> !iszero(M[r, c]), c:n)
        p === nothing && error("singular matrix in exact solve")
        p += c - 1
        if p != c
            M[c, :], M[p, :] = M[p, :], M[c, :]
        end
        piv = inv(M[c, c])
        for l in c:m
            M[c, l] *= piv
        end
        for r in 1:n
            r == c && continue
            f = M[r, c]
            iszero(f) && continue
            for l in c:m
                M[r, l] -= f * M[c, l]
            end
        end
    end
    return M[:, (n + 1):m]
end

"""
    _rref!(M) -> pivots

Reduced row echelon form in place (exact fields), returning the pivot columns.
"""
function _rref!(M::AbstractMatrix{T}) where {T}
    nr, nc = size(M)
    pivots = Int[]
    r = 1
    for c in 1:nc
        r > nr && break
        p = findfirst(i -> !iszero(M[i, c]), r:nr)
        p === nothing && continue
        p += r - 1
        if p != r
            for l in 1:nc
                M[r, l], M[p, l] = M[p, l], M[r, l]
            end
        end
        piv = inv(M[r, c])
        for l in c:nc
            M[r, l] *= piv
        end
        for i in 1:nr
            i == r && continue
            f = M[i, c]
            iszero(f) && continue
            for l in c:nc
                M[i, l] -= f * M[r, l]
            end
        end
        push!(pivots, c)
        r += 1
    end
    return pivots
end

"""
    _nullspace(A) -> Matrix

Basis of the right null space of `A` (columns).  Exact fields: the standard basis
attached to the reduced row echelon form.  Floats: from the SVD.
"""
function _nullspace(A::AbstractMatrix{T}) where {T <: AbstractFloat}
    m, n = size(A)
    (m == 0 || n == 0) && return Matrix{T}(I, n, n)
    F = svd(Matrix(A); full = true)
    tol = _tol(T) * max(1, maximum(F.S; init = zero(T)))
    r = count(>(tol), F.S)
    return Matrix(F.V[:, (r + 1):n])
end
function _nullspace(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    m == 0 && return Matrix{T}(I, n, n)
    M = Matrix{T}(A)
    piv = _rref!(M)
    free = setdiff(1:n, piv)
    Z = zeros(T, n, length(free))
    for (k, f) in enumerate(free)
        Z[f, k] = one(T)
        for (r, p) in enumerate(piv)
            Z[p, k] = -M[r, f]
        end
    end
    return Z
end

"""
    _left_inverse(A) -> L

A left inverse (`L * A = I`) of a matrix with full column rank.  On a consistent
system `A X = B` every left inverse gives the same `X = L B`.
"""
function _left_inverse(A::AbstractMatrix{T}) where {T <: AbstractFloat}
    size(A, 2) == 0 && return zeros(T, 0, size(A, 1))
    return pinv(Matrix(A))
end
function _left_inverse(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    n == 0 && return zeros(T, 0, m)
    M = Matrix{T}(transpose(A))
    rows = _rref!(M)                      # independent rows of A
    length(rows) == n || error("matrix does not have full column rank")
    L = zeros(T, n, m)
    L[:, rows] = _solve(Matrix{T}(A[rows, :]), Matrix{T}(I, n, n))
    return L
end

"""
    _right_inverse(S) -> R

A right inverse (`S * R = I`) of a matrix with full row rank.
"""
_right_inverse(S::AbstractMatrix{T}) where {T} = Matrix{T}(transpose(_left_inverse(transpose(S))))
