using Test
using Kiwi
using LinearAlgebra
using SparseArrays

const RK = Kiwi

# ---------------------------------------------------------------------------
# Closed forms for su(2) (Condon–Shortley), used as references
# ---------------------------------------------------------------------------

_fact(x) = factorial(big(Int(x)))

function racah_cg(j1, m1, j2, m2, J, M)
    m1 + m2 == M || return 0.0
    pre = sqrt((2J + 1) * _fact(J + j1 - j2) * _fact(J - j1 + j2) * _fact(j1 + j2 - J) / _fact(j1 + j2 + J + 1)) *
          sqrt(_fact(J + M) * _fact(J - M) * _fact(j1 - m1) * _fact(j1 + m1) * _fact(j2 - m2) * _fact(j2 + m2))
    s = big(0) // 1
    for k in 0:Int(j1 + j2 + J)
        args = (k, j1 + j2 - J - k, j1 - m1 - k, j2 + m2 - k, J - j2 + m1 + k, J - j1 - m2 + k)
        any(<(0), args) && continue
        s += big(-1)^k // prod(_fact.(args))
    end
    return Float64(pre * s)
end

wigner3j(j1, m1, j2, m2, j3, m3) =
    (-1)^Int(j1 - j2 - m3) / sqrt(2j3 + 1) * racah_cg(j1, m1, j2, m2, j3, -m3)

function wigner6j(a, b, c, d, e, f)   # {a b c; d e f}
    tri(x, y, z) = abs(x - y) <= z <= x + y && isinteger(x + y + z)
    (tri(a, b, c) && tri(a, e, f) && tri(d, b, f) && tri(d, e, c)) || return 0.0
    Δ(x, y, z) = _fact(x + y - z) * _fact(x - y + z) * _fact(-x + y + z) // _fact(x + y + z + 1)
    pre = sqrt(Δ(a, b, c) * Δ(a, e, f) * Δ(d, b, f) * Δ(d, e, c))
    s = big(0) // 1
    for t in 0:(Int(floor(a + b + c + d + e + f)) + 1)
        args = (t - a - b - c, t - a - e - f, t - d - b - f, t - d - e - c,
                a + b + d + e - t, a + c + d + f - t, b + c + e + f - t)
        any(<(0), args) && continue
        s += big(-1)^t * _fact(t + 1) // prod(_fact.(args))
    end
    return Float64(pre * s)
end

J2(x) = Irrep(SU(2), Int(2x))                 # spin j
m_of(j, i) = j - (i - 1)                      # basis index i ↔ m (highest weight first)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function relations_error(rep; basis = :lowering, T = Float64)
    E, F, H = generators(rep; basis, T)
    C = cartan_matrix(rep.algebra)
    err = 0.0
    for i in eachindex(E), j in eachindex(E)
        err = max(err, Float64(norm(E[i] * F[j] - F[j] * E[i] - (i == j ? H[i] : 0 * H[i]), Inf)))
        err = max(err, Float64(norm(H[i] * E[j] - E[j] * H[i] - C[j, i] * E[j], Inf)))
        if i != j
            Z = E[j]
            for _ in 1:(1 - C[j, i])
                Z = E[i] * Z - Z * E[i]
            end
            err = max(err, Float64(norm(Z, Inf)))
        end
    end
    return err
end

# all irreps e ∈ a ⊗ b, d ∈ e ⊗ c, f ∈ b ⊗ c with N^{af}_d > 0
function fusion_channels(a, b, c)
    g = a.algebra
    out = []
    for (el, _) in RK._fusion_rules(a, b), (dl, _) in RK._fusion_rules(Irrep(g, el), c)
        for (fl, _) in RK._fusion_rules(b, c)
            nsymbol(a, Irrep(g, fl), Irrep(g, dl)) > 0 && push!(out, (Irrep(g, el), Irrep(g, dl), Irrep(g, fl)))
        end
    end
    return out
end
fusions(a, b) = [Irrep(a.algebra, l) for l in keys(RK._fusion_rules(a, b))]

# Pentagon and hexagon equations with multiplicities, transcribed from
# TensorKitSectors.pentagon_equation / hexagon_equation (TensorKit index conventions).
function pentagon_error(a, b, c, d)
    err = 0.0
    for f in fusions(a, b), h in fusions(c, d), g in fusions(f, c), i in fusions(b, h)
        for e in intersect(fusions(g, d), fusions(a, i))
            F1 = fsymbol(f, c, d, e, g, h); F2 = fsymbol(a, b, h, e, f, i)
            sz = (size(F1, 1), size(F1, 2), size(F1, 3), size(F2, 1), size(F2, 3), size(F2, 4))
            p1 = zeros(sz); p2 = zeros(sz)
            for idx in CartesianIndices(p1)
                λ, μ, ν, κ, ρ, σ = Tuple(idx)
                p1[idx] = sum(F1[λ, μ, ν, τ] * F2[κ, τ, ρ, σ] for τ in axes(F1, 4); init = 0.0)
            end
            for j in fusions(b, c)
                G1 = fsymbol(a, b, c, g, f, j); G2 = fsymbol(a, j, d, e, g, i); G3 = fsymbol(b, c, d, i, j, h)
                (any(iszero, size(G1)) || any(iszero, size(G2)) || any(iszero, size(G3))) && continue
                for idx in CartesianIndices(p2)
                    λ, μ, ν, κ, ρ, σ = Tuple(idx)
                    p2[idx] += sum(G1[κ, λ, α, β] * G2[β, μ, τ, σ] * G3[α, τ, ν, ρ]
                                   for α in axes(G1, 3), β in axes(G1, 4), τ in axes(G2, 3); init = 0.0)
                end
            end
            err = max(err, maximum(abs.(p1 - p2); init = 0.0))
        end
    end
    return err
end

function hexagon_error(a, b, c)
    err = 0.0
    for e in fusions(c, a), f in fusions(c, b), d in intersect(fusions(e, b), fusions(a, f))
        F = fsymbol(a, c, b, d, e, f)
        for (Rx, Ry, Rz) in ((rsymbol(c, a, e), rsymbol(c, b, f), g -> rsymbol(c, g, d)),
                             (transpose(rsymbol(a, c, e)), transpose(rsymbol(b, c, f)), g -> transpose(rsymbol(g, c, d))))
            sz = (size(Rx, 1), size(F, 2), size(Ry, 2), size(F, 4))
            lhs = zeros(sz); rhs = zeros(sz)
            for idx in CartesianIndices(lhs)
                α, β, μ, ν = Tuple(idx)
                lhs[idx] = sum(Rx[α, λ] * F[λ, β, γ, ν] * Ry[γ, μ] for λ in axes(F, 1), γ in axes(F, 3))
            end
            for g in fusions(a, b)
                G1 = fsymbol(c, a, b, d, e, g); G2 = fsymbol(a, b, c, d, g, f); R = Rz(g)
                (isempty(G1) || isempty(G2)) && continue
                for idx in CartesianIndices(rhs)
                    α, β, μ, ν = Tuple(idx)
                    rhs[idx] += sum(G1[α, β, δ, σ] * R[σ, ψ] * G2[δ, ψ, μ, ν]
                                    for δ in axes(G1, 3), σ in axes(G1, 4), ψ in axes(G2, 2); init = 0.0)
                end
            end
            err = max(err, maximum(abs.(lhs - rhs); init = 0.0))
        end
    end
    return err
end

# Wigner 6j symbol {j1 j2 j12; j3 j j23} from Kiwi's F-symbol
function kiwi6j(j1, j2, j12, j3, j, j23)
    tri(x, y, z) = abs(x - y) <= z <= x + y && isinteger(x + y + z)
    (tri(j1, j2, j12) && tri(j12, j3, j) && tri(j2, j3, j23) && tri(j1, j23, j)) || return 0.0
    F = fsymbol(J2(j1), J2(j2), J2(j3), J2(j), J2(j12), J2(j23))[1]
    return (-1)^Int(j1 + j2 + j3 + j) * F / sqrt((2j12 + 1) * (2j23 + 1))
end

# basis change W (orthogonal) from the lowering basis to the GT basis of an su(N) irrep:
# W e_hw = e_hw and W F_low = F_gt W
function lowering_to_gt(rep)
    L = irrep_module(rep; basis = :lowering); G = irrep_module(rep; basis = :gelfand_tsetlin)
    ws = L.ws
    W = Vector{Matrix{Float64}}(undef, length(ws))
    W[1] = ones(1, 1)
    for k in 2:length(ws)
        ups = [i for i in 1:size(ws.up, 1) if ws.up[i, k] != 0]
        S = reduce(hcat, (L.F[i, k] for i in ups))
        T = reduce(hcat, (G.F[i, k] * W[ws.up[i, k]] for i in ups))
        W[k] = T * RK._right_inverse(S)
    end
    d = sum(ws.mult)
    Wfull = zeros(d, d)
    for k in 1:length(ws)
        r = (ws.offset[k] + 1):(ws.offset[k] + ws.mult[k])
        Wfull[r, r] = W[k]
    end
    return Wfull
end

# multiplicity-space gauge U with C_gt = (W_a ⊗ W_b) C_low W_cᵀ U
function gt_gauge(a, b, c, Wa, Wb, Wc)
    Cl = fusion_tensor(a, b, c); Cg = fusion_tensor(a, b, c; basis = :gelfand_tsetlin)
    da, db, dc, N = size(Cl)
    U = zeros(N, N)
    for μ in 1:N, ν in 1:N
        X = reshape(Wa * reshape(Cl[:, :, :, μ], da, db * dc), da, db, dc)
        X = permutedims(reshape(Wb * reshape(permutedims(X, (2, 1, 3)), db, da * dc), db, da, dc), (2, 1, 3))
        X = reshape(reshape(X, da * db, dc) * transpose(Wc), da, db, dc)
        U[μ, ν] = dot(vec(X), vec(Cg[:, :, :, ν])) / dc
    end
    return U
end

# ---------------------------------------------------------------------------

@testset "Recoupling" begin
    @testset "Modules: Chevalley–Serre relations" begin
        reps = [Irrep(SU(2), 4), Irrep(SU(3), 2, 1), Irrep(SU(4), 1, 0, 1), Irrep(SO(7), 1, 0, 1),
                Irrep(Sp(6), 0, 1, 1), Irrep(SO(8), 1, 0, 1, 0), Irrep(G_series(2), 1, 1),
                Irrep(F_series(4), 0, 0, 0, 1), Irrep(E_series(6), 0, 0, 0, 0, 0, 1)]
        for rep in reps
            @test relations_error(rep) < 1e-11
            E, F, _ = generators(rep)
            @test all(E[i] == sparse(transpose(F[i])) for i in eachindex(E))
            @test size(E[1], 1) == dimension(rep)
        end
        for rep in [Irrep(SU(3), 2, 1), Irrep(G_series(2), 1, 0), Irrep(Sp(4), 1, 1)]
            @test relations_error(rep; T = Rational{BigInt}) == 0
        end
        for rep in [Irrep(SU(3), 2, 1), Irrep(SU(4), 1, 1, 0), Irrep(SU(5), 0, 1, 0, 1)]
            @test relations_error(rep; basis = :gelfand_tsetlin) < 1e-11
            @test length(gt_patterns(rep)) == dimension(rep)
        end
        # su(2): Condon–Shortley, F = J₋
        _, F, _ = generators(J2(3 // 2))
        @test Matrix(F[1]) ≈ [0 0 0 0; √3 0 0 0; 0 2 0 0; 0 0 √3 0]
        # the exact (Shapovalov) basis orthonormalises to the floating-point basis
        for rep in [Irrep(SU(3), 2, 1), Irrep(G_series(2), 1, 1), Irrep(SO(7), 0, 1, 0)]
            Mx = irrep_module(rep; T = Rational{BigInt}); Mf = irrep_module(rep)
            err = 0.0
            for k in 1:length(Mx.ws), i in 1:lie_rank(rep.algebra)
                ku = Mx.ws.up[i, k]
                ku == 0 && continue
                Rk = cholesky(Symmetric(Float64.(Mx.gram[k]))).U
                Ru = cholesky(Symmetric(Float64.(Mx.gram[ku]))).U
                err = max(err, norm(Rk * Float64.(Mx.F[i, k]) / Ru - Mf.F[i, k], Inf))
            end
            @test err < 1e-12
        end
    end

    @testset "Clebsch–Gordan tensors" begin
        cases = [(Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 1)), (Irrep(SU(3), 2, 1), Irrep(SU(3), 1, 1)),
                 (Irrep(G_series(2), 1, 0), Irrep(G_series(2), 0, 1)), (Irrep(SO(7), 0, 0, 1), Irrep(SO(7), 1, 0, 0)),
                 (Irrep(Sp(4), 1, 1), Irrep(Sp(4), 1, 0)), (Irrep(SO(8), 0, 0, 1, 0), Irrep(SO(8), 0, 0, 0, 1)),
                 (Irrep(F_series(4), 0, 0, 0, 1), Irrep(F_series(4), 0, 0, 0, 1))]
        for (a, b) in cases, basis in (a.algebra.series == :A ? (:lowering, :gelfand_tsetlin) : (:lowering,))
            Ea, Fa, _ = generators(a; basis); Eb, Fb, _ = generators(b; basis)
            da, db = dimension(a), dimension(b)
            Ia, Ib = sparse(I, da, da), sparse(I, db, db)
            blocks = Matrix{Float64}[]
            ok = true
            for c in fusions(a, b)
                C = fusion_tensor(a, b, c; basis)
                ok &= size(C, 4) == nsymbol(a, b, c) == tensor_product(a, b)[c]
                Ec, Fc, _ = generators(c; basis)
                for μ in axes(C, 4)
                    M = reshape(C[:, :, :, μ], da * db, size(C, 3))
                    for i in eachindex(Ea)
                        ok &= norm((kron(Ib, Ea[i]) + kron(Eb[i], Ia)) * M - M * Ec[i], Inf) < 1e-11
                        ok &= norm((kron(Ib, Fa[i]) + kron(Fb[i], Ia)) * M - M * Fc[i], Inf) < 1e-11
                    end
                    push!(blocks, M)
                end
            end
            U = reduce(hcat, blocks)
            @test ok
            @test size(U) == (da * db, da * db)
            @test norm(U' * U - I, Inf) < 1e-11            # orthonormal and complete
        end
        # su(2) against the Racah formula
        err = 0.0
        for j1 in 0:1//2:3, j2 in 0:1//2:3, J in abs(j1 - j2):(j1 + j2)
            C = fusion_tensor(J2(j1), J2(j2), J2(J))
            for i in axes(C, 1), k in axes(C, 2), l in axes(C, 3)
                err = max(err, abs(C[i, k, l, 1] - racah_cg(j1, m_of(j1, i), j2, m_of(j2, k), J, m_of(J, l))))
            end
        end
        @test err < 1e-12
    end

    @testset "1j, 2j and 3j symbols" begin
        reps = [J2(1 // 2), J2(1), J2(3 // 2), Irrep(SU(3), 1, 0), Irrep(SU(3), 2, 1), Irrep(SU(3), 1, 1),
                Irrep(SU(4), 0, 1, 0), Irrep(SU(4), 1, 0, 0), Irrep(Sp(4), 1, 0), Irrep(SO(7), 0, 0, 1),
                Irrep(SO(8), 0, 0, 1, 0), Irrep(SO(10), 0, 0, 0, 1, 0), Irrep(G_series(2), 1, 0),
                Irrep(E_series(6), 1, 0, 0, 0, 0, 0), Irrep(E_series(7), 0, 0, 0, 0, 0, 1, 0)]
        for a in reps
            Z = twoj(a); Zb = twoj(conjugate(a))
            κ = onej(a)
            @test norm(Z' * Z - I, Inf) < 1e-11
            @test Z[1, end] > 0                               # Z[hw(a), lw(ā)] > 0
            @test norm(Zb - κ * transpose(Z), Inf) < 1e-11
            @test onej(conjugate(a)) == κ
            is_self_conjugate(a) && @test κ == frobenius_schur_indicator(a)
            u = RK._unit(a.algebra)
            @test κ == sign(fsymbol(a, conjugate(a), a, a, u, u)[1])   # TensorKit's FS phase
        end
        # su(2): 2j symbol (−1)^{j−m} δ_{m',−m} and Wigner's 3j symbol
        for j in 0:1//2:3
            Z = twoj(J2(j))
            @test Z ≈ [(m_of(j, k) == -m_of(j, i)) * (-1)^Int(j - m_of(j, i)) for i in 1:Int(2j + 1), k in 1:Int(2j + 1)]
        end
        err = 0.0
        for j1 in 0:1//2:2, j2 in 0:1//2:2, j3 in abs(j1 - j2):(j1 + j2)
            T3 = threej(J2(j1), J2(j2), J2(j3))
            ph = (-1)^Int(j1 - j2 - j3)
            for i in axes(T3, 1), k in axes(T3, 2), l in axes(T3, 3)
                err = max(err, abs(ph * T3[i, k, l, 1] - wigner3j(j1, m_of(j1, i), j2, m_of(j2, k), j3, m_of(j3, l))))
            end
        end
        @test err < 1e-12
        # 3j symbols are normalised invariants (including outer multiplicity)
        for (a, b, c) in [(Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 1)),
                          (Irrep(SU(3), 1, 0), Irrep(SU(3), 1, 0), Irrep(SU(3), 1, 0)),
                          (Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0), Irrep(G_series(2), 0, 1))]
            T3 = threej(a, b, c)
            Ea, _, _ = generators(a); Eb, _, _ = generators(b); Ec, _, _ = generators(c)
            da, db, dc = size(T3, 1), size(T3, 2), size(T3, 3)
            for μ in axes(T3, 4)
                v = vec(T3[:, :, :, μ])
                @test norm(v) ≈ 1
                for i in eachindex(Ea)
                    Δ = kron(sparse(I, dc, dc), kron(sparse(I, db, db), Ea[i])) +
                        kron(sparse(I, dc, dc), kron(Eb[i], sparse(I, da, da))) +
                        kron(Ec[i], kron(sparse(I, db, db), sparse(I, da, da)))
                    @test norm(Δ * v) < 1e-11
                end
            end
        end
    end

    @testset "F-symbols" begin
        # su(2): F^{j1 j2 j3}_j[j12, j23] = (−1)^{j1+j2+j3+j} √((2j12+1)(2j23+1)) {j1 j2 j12; j3 j j23}
        err = 0.0
        for j1 in 0:1//2:3//2, j2 in 0:1//2:3//2, j3 in 0:1//2:3//2, j12 in abs(j1 - j2):(j1 + j2),
            j23 in abs(j2 - j3):(j2 + j3), j in abs(j12 - j3):(j12 + j3)
            abs(j1 - j23) <= j <= j1 + j23 || continue
            F = fsymbol(J2(j1), J2(j2), J2(j3), J2(j), J2(j12), J2(j23))[1]
            W = (-1)^Int(j1 + j2 + j3 + j) * sqrt((2j12 + 1) * (2j23 + 1)) * wigner6j(j1, j2, j12, j3, j, j23)
            err = max(err, abs(F - W))
        end
        @test err < 1e-12
        # su(2): Biedenharn–Elliott identity and Racah's sum rule, in Wigner-6j form,
        # evaluated with Kiwi's F-symbols
        be = 0.0; rs = 0.0
        js = 0:1//2:1
        for a in js, b in js, c in js, d in js, e in js, f in js, p in js, q in js, r in js
            S = a + b + c + d + e + f + p + q + r
            lhs = sum((-1)^Int(S + x) * (2x + 1) * kiwi6j(a, b, x, c, d, p) * kiwi6j(c, d, x, e, f, q) *
                      kiwi6j(e, f, x, b, a, r) for x in 0:1//2:2 if isinteger(S + x); init = 0.0)
            be = max(be, abs(lhs - kiwi6j(p, q, r, e, a, d) * kiwi6j(p, q, r, f, b, c)))
        end
        for a in js, b in js, c in js, d in js, p in 0:1//2:2, q in 0:1//2:2
            isinteger(a + b + c + d) || continue
            lhs = sum((-1)^Int(p + q + x) * (2x + 1) * kiwi6j(a, b, x, c, d, p) * kiwi6j(a, b, x, d, c, q)
                      for x in 0:1//2:2 if isinteger(p + q + x); init = 0.0)
            rs = max(rs, abs(lhs - kiwi6j(a, c, q, b, d, p)))
        end
        @test be < 1e-12
        @test rs < 1e-12
        # the same F-symbols against the closed form (independent of Kiwi)
        @test kiwi6j(1, 1, 1, 1, 1, 1) ≈ wigner6j(1, 1, 1, 1, 1, 1)
        # pentagon (Biedenharn–Elliott in TensorKit form) for su(2)
        for (a, b, c, d) in [(J2(1 // 2), J2(1), J2(3 // 2), J2(1)), (J2(1), J2(1), J2(1), J2(1)),
                             (J2(3 // 2), J2(1 // 2), J2(1), J2(3 // 2))]
            @test pentagon_error(a, b, c, d) < 1e-12
        end
        # pentagon with outer multiplicities
        su3 = (Irrep(SU(3), 1, 0), Irrep(SU(3), 0, 1), Irrep(SU(3), 1, 1))
        @test pentagon_error(su3[3], su3[3], su3[3], su3[1]) < 1e-12
        @test pentagon_error(su3[1], su3[3], su3[2], su3[3]) < 1e-12
        @test pentagon_error(Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0),
                             Irrep(G_series(2), 0, 1), Irrep(G_series(2), 1, 0)) < 1e-12
        @test pentagon_error(Irrep(Sp(4), 1, 0), Irrep(Sp(4), 0, 1), Irrep(Sp(4), 1, 0), Irrep(Sp(4), 0, 1)) < 1e-12
        # unitarity of the F-move, block by block
        for (a, b, c) in [(su3[3], su3[3], su3[3]), (su3[1], su3[3], su3[2]),
                          (Irrep(SO(7), 0, 0, 1), Irrep(SO(7), 1, 0, 0), Irrep(SO(7), 0, 0, 1))]
            for d in unique(reduce(vcat, [fusions(e, c) for e in fusions(a, b)]))
                es = [e for e in fusions(a, b) if nsymbol(e, c, d) > 0]
                fs = [f for f in fusions(b, c) if nsymbol(a, f, d) > 0]
                blocks = [reshape(fsymbol(a, b, c, d, e, f), :, nsymbol(b, c, f) * nsymbol(a, f, d)) for e in es, f in fs]
                F = reduce(vcat, [reduce(hcat, blocks[i, :]) for i in axes(blocks, 1)])
                @test norm(F' * F - I, Inf) < 1e-12
            end
        end
        # the highest-weight evaluation equals TensorKit's full contraction
        a, b, c = su3[3], su3[3], su3[1]
        for (e, d, f) in fusion_channels(a, b, c)
            A = fusion_tensor(a, b, e); B = fusion_tensor(e, c, d); C = fusion_tensor(b, c, f); D = fusion_tensor(a, f, d)
            F = fsymbol(a, b, c, d, e, f)
            k = dimension(d)                        # any basis vector of d (here the lowest)
            G = zeros(size(F))
            for μ in axes(F, 1), ν in axes(F, 2), κ in axes(F, 3), λ in axes(F, 4)
                G[μ, ν, κ, λ] = sum(D[i, n, k, λ] * C[j, l, n, κ] * A[i, j, m, μ] * B[m, l, k, ν]
                                    for i in axes(A, 1), j in axes(A, 2), m in axes(A, 3), l in axes(C, 2), n in axes(C, 3))
            end
            @test G ≈ F atol = 1e-12
        end
    end

    @testset "R-symbols and hexagon" begin
        su3 = (Irrep(SU(3), 1, 0), Irrep(SU(3), 0, 1), Irrep(SU(3), 1, 1))
        for (a, b) in [(su3[3], su3[3]), (su3[1], su3[3]), (J2(1), J2(3 // 2)), (Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0))]
            for c in fusions(a, b)
                R1 = rsymbol(a, b, c); R2 = rsymbol(b, a, c)
                @test norm(R1 * R2 - I, Inf) < 1e-12        # symmetric braiding
            end
        end
        # su(2): R^{j1 j2}_J = (−1)^{j1+j2−J}
        for j1 in 0:1//2:2, j2 in 0:1//2:2, J in abs(j1 - j2):(j1 + j2)
            @test rsymbol(J2(j1), J2(j2), J2(J))[1] ≈ (-1)^Int(j1 + j2 - J)
        end
        for (a, b, c) in [(J2(1 // 2), J2(1), J2(3 // 2)), (su3[1], su3[1], su3[2]), (su3[3], su3[3], su3[1]),
                          (Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0), Irrep(G_series(2), 0, 1))]
            @test hexagon_error(a, b, c) < 1e-12
        end
    end

    @testset "A and B symbols" begin
        for (a, b) in [(Irrep(SU(3), 1, 0), Irrep(SU(3), 1, 1)), (Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 1)),
                       (J2(1 // 2), J2(3 // 2)), (Irrep(Sp(4), 1, 0), Irrep(Sp(4), 1, 0))]
            ab, bb = conjugate(a), conjugate(b)
            u = RK._unit(a.algebra)
            Za, Zb = twoj(a), twoj(b)
            for c in fusions(a, b)
                # A from fusion tensors (TensorKit's Asymbol_from_fusiontensor)
                C1 = fusion_tensor(a, b, c)[:, 1, :, :]; C2 = fusion_tensor(ab, c, b)[:, :, 1, :]
                A2 = [sqrt(dimension(b) / dimension(c)) * sum(Za[i, j] * C1[i, k, μ] * C2[j, k, ν]
                                                               for i in axes(Za, 1), j in axes(Za, 2), k in axes(C1, 2))
                      for μ in axes(C1, 3), ν in axes(C2, 3)]
                @test asymbol(a, b, c) ≈ A2 atol = 1e-12
                D1 = fusion_tensor(a, b, c)[1, :, :, :]; D2 = fusion_tensor(c, bb, a)[:, :, 1, :]
                B2 = [sqrt(dimension(a) / dimension(c)) * sum(Zb[i, j] * D1[i, k, μ] * D2[k, j, ν]
                                                               for i in axes(Zb, 1), j in axes(Zb, 2), k in axes(D1, 2))
                      for μ in axes(D1, 3), ν in axes(D2, 3)]
                @test bsymbol(a, b, c) ≈ B2 atol = 1e-12
            end
        end
    end

    @testset "Gelfand–Tsetlin and lowering bases agree up to gauge" begin
        a, b, c = Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 1), Irrep(SU(3), 1, 0)
        W = Dict{Vector{Int}, Matrix{Float64}}()
        Wof(r) = get!(() -> lowering_to_gt(r), W, r.dynkin_labels)
        for r in (a, Irrep(SU(3), 2, 1), Irrep(SU(3), 0, 2))
            Wr = Wof(r)
            @test norm(Wr' * Wr - I, Inf) < 1e-12
            El, Fl, _ = generators(r); Eg, Fg, _ = generators(r; basis = :gelfand_tsetlin)
            @test all(norm(Wr * Fl[i] - Fg[i] * Wr, Inf) < 1e-12 for i in 1:2)
        end
        U = Dict{Any, Matrix{Float64}}()
        Uof(x, y, z) = get!(() -> gt_gauge(x, y, z, Wof(x), Wof(y), Wof(z)), U, (x, y, z))
        for (e, d, f) in fusion_channels(a, b, c)
            Ua, Ub, Uc, Ud = Uof(a, b, e), Uof(e, c, d), Uof(b, c, f), Uof(a, f, d)
            for Ux in (Ua, Ub, Uc, Ud)
                @test norm(Ux' * Ux - I, Inf) < 1e-12
            end
            Fl = fsymbol(a, b, c, d, e, f); Fg = fsymbol(a, b, c, d, e, f; basis = :gelfand_tsetlin)
            G = zeros(size(Fl))
            for μ in axes(G, 1), ν in axes(G, 2), κ in axes(G, 3), λ in axes(G, 4)
                G[μ, ν, κ, λ] = sum(Ua[μ′, μ] * Ub[ν′, ν] * Uc[κ′, κ] * Ud[λ′, λ] * Fl[μ′, ν′, κ′, λ′]
                                    for μ′ in axes(G, 1), ν′ in axes(G, 2), κ′ in axes(G, 3), λ′ in axes(G, 4))
            end
            @test G ≈ Fg atol = 1e-12
        end
    end

    @testset "Exact arithmetic" begin
        # exact su(2) F-symbols are the signed square roots of rationals
        for (j1, j2, j3, j, j12, j23) in [(1, 1, 1, 1, 1, 1), (1 // 2, 1, 3 // 2, 1, 1 // 2, 3 // 2), (2, 1, 1, 2, 1, 2)]
            Fx = fsymbol(J2(j1), J2(j2), J2(j3), J2(j), J2(j12), J2(j23); T = Rational{BigInt})[1]
            Ff = fsymbol(J2(j1), J2(j2), J2(j3), J2(j), J2(j12), J2(j23))[1]
            @test Fx isa SqrtRational
            @test float(Fx) ≈ Ff atol = 1e-14
        end
        # multiplicity-free F- and R-symbols of su(3), G₂ and so(7): exact = floating point
        for (a, b, c) in [(Irrep(SU(3), 1, 0), Irrep(SU(3), 1, 0), Irrep(SU(3), 0, 1)),
                          (Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0), Irrep(G_series(2), 1, 0)),
                          (Irrep(SO(7), 0, 0, 1), Irrep(SO(7), 0, 0, 1), Irrep(SO(7), 1, 0, 0))]
            for (e, d, f) in fusion_channels(a, b, c)
                all(n == 1 for n in (nsymbol(a, b, e), nsymbol(e, c, d), nsymbol(b, c, f), nsymbol(a, f, d))) || continue
                @test float(fsymbol(a, b, c, d, e, f; T = Rational{BigInt})[1]) ≈ fsymbol(a, b, c, d, e, f)[1] atol = 1e-13
            end
            for e in fusions(a, b)
                nsymbol(a, b, e) == 1 || continue
                @test float(rsymbol(a, b, e; T = Rational{BigInt})[1]) ≈ rsymbol(a, b, e)[1] atol = 1e-13
            end
        end
        @test SqrtRational(-1, 3 // 4) == -SqrtRational(1, 3 // 4)
        @test float(SqrtRational(-1, 1 // 4)) == -0.5
    end
end
