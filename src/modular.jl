"""
Multi-modular integers used for exact computations whose results exceed 128 bits.

A `_MM{K}` value stores an integer modulo `2^128` (plain wrapping `Int128`
arithmetic) together with its residues modulo `K` primes just below `2^31`.
Ring operations are cheap and allocation-free; the exact integer is recovered at
the end with the Chinese remainder theorem, provided its absolute value is less
than half the product of the moduli.
"""

const _MM_PRIMES = Int64[]
const _MM_LOCK = ReentrantLock()

function _isprime_small(p::Int64)
    p < 2 && return false
    p % 2 == 0 && return p == 2
    d = 3
    while d * d <= p
        p % d == 0 && return false
        d += 2
    end
    return true
end

function _mm_primes(K::Int)
    lock(_MM_LOCK) do
        while length(_MM_PRIMES) < K
            p = isempty(_MM_PRIMES) ? Int64(2)^31 - 1 : _MM_PRIMES[end] - 2
            while !_isprime_small(p)
                p -= 2
            end
            push!(_MM_PRIMES, p)
        end
    end
    return _MM_PRIMES
end

struct _MM{K}
    w::Int128
    r::NTuple{K, Int64}
end

@inline _mmp(k) = @inbounds _MM_PRIMES[k]

function _MM{K}(x::Integer) where {K}
    _mm_primes(K)
    w = x isa Union{Int64, Int128} ? Int128(x) :
        Int128(mod(big(x) + big(2)^127, big(2)^128) - big(2)^127)
    return _MM{K}(w, ntuple(k -> Int64(mod(x, _mmp(k))), Val(K)))
end
_MM{K}(x::_MM{K}) where {K} = x

Base.zero(::Type{_MM{K}}) where {K} = _MM{K}(Int128(0), ntuple(_ -> Int64(0), Val(K)))
Base.one(::Type{_MM{K}}) where {K} = _MM{K}(Int128(1), ntuple(_ -> Int64(1), Val(K)))
Base.zero(::_MM{K}) where {K} = zero(_MM{K})

@inline function Base.:+(a::_MM{K}, b::_MM{K}) where {K}
    _MM{K}(a.w + b.w, ntuple(Val(K)) do k
        s = a.r[k] + b.r[k]
        p = _mmp(k)
        s >= p ? s - p : s
    end)
end

@inline function Base.:-(a::_MM{K}, b::_MM{K}) where {K}
    _MM{K}(a.w - b.w, ntuple(Val(K)) do k
        s = a.r[k] - b.r[k]
        s < 0 ? s + _mmp(k) : s
    end)
end

@inline function Base.:-(a::_MM{K}) where {K}
    _MM{K}(-a.w, ntuple(Val(K)) do k
        a.r[k] == 0 ? Int64(0) : _mmp(k) - a.r[k]
    end)
end

@inline function Base.:*(a::_MM{K}, b::_MM{K}) where {K}
    _MM{K}(a.w * b.w, ntuple(k -> (a.r[k] * b.r[k]) % _mmp(k), Val(K)))
end

Base.iszero(a::_MM) = iszero(a.w) && all(iszero, a.r)
Base.:(==)(a::_MM{K}, b::_MM{K}) where {K} = a.w == b.w && a.r == b.r
Base.:(==)(a::_MM{K}, b::Integer) where {K} = a == _MM{K}(b)

"""
    _to_bigint(a::_MM) -> BigInt

Recover the integer represented by `a` (assumed to lie in the symmetric range).
"""
function _to_bigint(a::_MM{K}) where {K}
    M = big(2)^128
    x = big(reinterpret(UInt128, a.w))
    for k in 1:K
        p = _mmp(k)
        t = mod((a.r[k] - mod(x, p)) * invmod(mod(M, p), p), p)
        x += M * t
        M *= p
    end
    return x > M ÷ 2 ? x - M : x
end
_to_bigint(a::Integer) = big(a)

"""
    _exact_int_type(bits) -> Type

An integer-like type that represents integers of absolute value below `2^bits`
exactly: `Int`, `Int128`, or a multi-modular `_MM{K}`.
"""
function _exact_int_type(bits::Integer)
    bits <= 62 && return Int
    bits <= 126 && return Int128
    K = cld(bits + 2 - 127, 30)
    _mm_primes(K)
    return _MM{K}
end
