# Gelfand–Tsetlin basis of su(N) irreps.
#
# A GT pattern has rows m[l] = (m_{1,l}, …, m_{l,l}), l = 1, …, N, with the top row
# m[N] the highest weight as a partition (m_{k,N} = Σ_{j≥k} a_j, m_{N,N} = 0) and the
# betweenness conditions m_{k,l} ≥ m_{k,l-1} ≥ m_{k+1,l}.  With σ_l = Σ_k m_{k,l}
# (σ₀ = 0) the weight has occupation numbers w_l = σ_l − σ_{l-1}, and Dynkin labels
# w_i − w_{i+1}.
#
# The raising operator Eₗ (simple root αₗ = εₗ − εₗ₊₁) raises one entry of row l:
#
#   Eₗ |m⟩ = Σ_k  √c_k(m) |m + δ_{k,l}⟩,
#
#   c_k = − ∏_{k'=1}^{l+1} (m_{k',l+1} − m_{k,l} + k − k')
#           ∏_{k'=1}^{l-1} (m_{k',l-1} − m_{k,l} + k − k' − 1)
#         / ∏_{k'≠k}      (m_{k',l} − m_{k,l} + k − k')(m_{k',l} − m_{k,l} + k − k' − 1),
#
# with nonnegative square roots, and Fₗ = Eₗᵀ (the basis is orthonormal).  This is
# the standard GT basis (e.g. Alex, Kalus, Huckleberry and von Delft, J. Math. Phys.
# 52, 023507 (2011)), the same as in SUNRepresentations.jl.  Within a weight space,
# patterns are ordered by decreasing (m[N-1], m[N-2], …, m[1]) lexicographically.

"""
    gt_patterns(rep::Irrep) -> Vector{Vector{Vector{Int}}}

All Gelfand–Tsetlin patterns of an su(N) irrep, as vectors of rows `m[1], …, m[N]`.
"""
function gt_patterns(rep::Irrep)
    rep.algebra.series == :A || throw(ArgumentError("Gelfand–Tsetlin patterns are defined for su(N)"))
    a = rep.dynkin_labels
    N = length(a) + 1
    top = [sum(a[k:end]; init = 0) for k in 1:N]
    out = Vector{Vector{Vector{Int}}}()
    rows = Vector{Vector{Int}}(undef, N)
    rows[N] = top
    function rec(l)
        if l == 0
            push!(out, [copy(r) for r in rows])
            return
        end
        above = rows[l + 1]
        ranges = [above[k + 1]:above[k] for k in 1:l]
        for t in Iterators.product(ranges...)
            rows[l] = collect(t)
            rec(l - 1)
        end
    end
    rec(N - 1)
    return out
end

function _gt_weight(m::Vector{Vector{Int}})
    N = length(m)
    σ = [l == 0 ? 0 : sum(m[l]) for l in 0:N]
    w = [σ[l + 1] - σ[l] for l in 1:N]
    return [w[i] - w[i + 1] for i in 1:(N - 1)]
end

function _gt_coefficient(m::Vector{Vector{Int}}, k::Int, l::Int)
    num = big(-1)
    mkl = m[l][k]
    for k′ in 1:(l + 1)
        num *= m[l + 1][k′] - mkl + k - k′
    end
    for k′ in 1:(l - 1)
        num *= m[l - 1][k′] - mkl + k - k′ - 1
    end
    den = big(1)
    for k′ in 1:l
        k′ == k && continue
        den *= (m[l][k′] - mkl + k - k′) * (m[l][k′] - mkl + k - k′ - 1)
    end
    return iszero(num) ? big(0) // 1 : num // den
end

function _gt_module(rep::Irrep, ::Type{T}) where {T}
    _isfloat(T) || error("the Gelfand–Tsetlin basis is orthonormal; use a floating-point type")
    ws = WeightStructure(rep)
    n = size(ws.up, 1)
    nw = length(ws)
    pats = gt_patterns(rep)
    byweight = [Vector{Vector{Vector{Int}}}() for _ in 1:nw]
    for m in pats
        push!(byweight[ws.index[_gt_weight(m)]], m)
    end
    key(m) = reduce(vcat, (m[l] for l in (length(m) - 1):-1:1); init = Int[])
    for k in 1:nw
        length(byweight[k]) == ws.mult[k] || error("GT pattern count disagrees with Freudenthal")
        sort!(byweight[k], by = key, rev = true)
    end
    pos = Dict{Vector{Vector{Int}}, Int}()
    for k in 1:nw, (j, m) in enumerate(byweight[k])
        pos[m] = j
    end
    empty = zeros(T, 0, 0)
    E = fill(empty, n, nw)
    F = fill(empty, n, nw)
    for k in 1:nw, l in 1:n
        ku = ws.up[l, k]
        ku == 0 && continue
        El = zeros(T, ws.mult[ku], ws.mult[k])
        for (j, m) in enumerate(byweight[k]), kk in 1:l
            c = _gt_coefficient(m, kk, l)
            iszero(c) && continue
            m′ = [copy(r) for r in m]
            m′[l][kk] += 1
            j′ = get(pos, m′, 0)
            j′ == 0 && continue
            El[j′, j] = sqrt(T(c))
        end
        E[l, k] = El
        F[l, k] = Matrix(transpose(El))
    end
    prov = [Tuple{Int, Int}[] for _ in 1:nw]
    M = IrrepModule{T}(ws, :gelfand_tsetlin, F, E, Matrix{T}[], prov, Dict{Symbol, Any}(), ReentrantLock())
    M.cache[:patterns] = byweight
    return M
end
