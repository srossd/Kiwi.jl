using Test
using Kiwi

# Dimension of a Schur functor applied to a d-dimensional space (hook-content formula)
function schur_dim(d::Integer, parts::Vector{Int})
    num = big(1); den = big(1)
    p = Partition(parts)
    for i in 1:length(parts), j in 1:parts[i]
        num *= d + j - i
        den *= hook_length(p, i, j)
    end
    return div(num, den)
end

rep_dim(R) = sum(big(dimension(ir)) * m for (ir, m) in R.components; init = big(0))

@testset "Fast algorithms" begin
    @testset "Conjugation and Frobenius–Schur indicator" begin
        e6 = E_series(6)
        @test conjugate(Irrep(e6, 0, 0, 0, 0, 0, 1)) == Irrep(e6, 0, 0, 0, 0, 0, 1)   # 78
        @test conjugate(Irrep(e6, 1, 0, 0, 0, 0, 0)) == Irrep(e6, 0, 0, 0, 0, 1, 0)   # 27 ↔ 27b
        @test conjugate(Irrep(e6, 0, 1, 0, 0, 0, 0)) == Irrep(e6, 0, 0, 0, 1, 0, 0)
        @test is_self_conjugate(Irrep(e6, 1, 0, 0, 0, 1, 0))                          # 650
        @test conjugate(Irrep(SU(4), 1, 2, 3)) == Irrep(SU(4), 3, 2, 1)
        @test conjugate(Irrep(SO(10), 0, 0, 0, 1, 0)) == Irrep(SO(10), 0, 0, 0, 0, 1)
        @test conjugate(Irrep(SO(8), 0, 0, 1, 0)) == Irrep(SO(8), 0, 0, 1, 0)

        @test frobenius_schur_indicator(Irrep(SU(2), 1)) == -1
        @test frobenius_schur_indicator(Irrep(SU(2), 2)) == 1
        @test frobenius_schur_indicator(Irrep(SU(3), 1, 0)) == 0
        @test frobenius_schur_indicator(Irrep(Sp(4), 1, 0)) == -1
        @test frobenius_schur_indicator(Irrep(Sp(4), 0, 1)) == 1
        @test frobenius_schur_indicator(Irrep(SO(5), 0, 1)) == -1
        @test frobenius_schur_indicator(Irrep(SO(7), 0, 0, 1)) == 1
        @test frobenius_schur_indicator(Irrep(SU(6), 0, 0, 1, 0, 0)) == -1
        @test frobenius_schur_indicator(Irrep(SO(12), 0, 0, 0, 0, 0, 1)) == -1
        @test frobenius_schur_indicator(Irrep(SO(16), 0, 0, 0, 0, 0, 0, 0, 1)) == 1
        @test frobenius_schur_indicator(Irrep(E_series(7), 0, 0, 0, 0, 0, 1, 0)) == -1
        @test is_real(adjoint_irrep(E_series(8)))
    end

    @testset "Dominant characters" begin
        for (g, l) in [(SU(4), [2, 1, 3]), (SO(9), [1, 0, 1, 1]), (Sp(8), [1, 0, 2, 0]),
                       (SO(10), [0, 1, 0, 1, 1]), (G_series(2), [3, 4]), (F_series(4), [1, 0, 1, 0]),
                       (E_series(6), [0, 1, 0, 0, 1, 1]), (E_series(7), [1, 0, 0, 0, 0, 1, 0]),
                       (E_series(8), [0, 0, 0, 0, 0, 0, 2, 0])]
            doms, mults = dominant_character(g, l)
            @test doms[1] == l && mults[1] == 1
            total = sum(big(m) * Kiwi._orbit_size(Kiwi.algebra_data(g), μ) for (μ, m) in zip(doms, mults))
            @test total == dimension(Irrep(g, l))
            # agrees with the full character
            if dimension(Irrep(g, l)) < 200_000
                c = character(Irrep(g, l))
                @test sum(values(c.weights)) == dimension(Irrep(g, l))
                @test all(c[Weight(g, μ)] == m for (μ, m) in zip(doms, mults))
            end
        end
        # E8 adjoint: 240 roots and 8 zero weights
        doms, mults = dominant_character(E_series(8), [0, 0, 0, 0, 0, 0, 1, 0])
        @test sort(mults) == [1, 8]
    end

    @testset "Tensor products (larger cases)" begin
        for ((g, a), b) in [((SU(4), [1, 2, 0]), [2, 0, 1]), ((SO(10), [0, 0, 0, 1, 1]), [1, 0, 0, 1, 0]),
                            ((G_series(2), [2, 3]), [1, 2]), ((F_series(4), [0, 0, 1, 0]), [1, 0, 0, 0]),
                            ((E_series(7), [0, 0, 0, 0, 0, 1, 0]), [0, 0, 0, 0, 0, 1, 0])]
            x = Irrep(g, a); y = Irrep(g, b)
            r = tensor_product(x, y)
            @test rep_dim(r) == big(dimension(x)) * dimension(y)
            @test r.components == tensor_product(y, x).components
            @test all(m > 0 for m in values(r.components))
        end
        # 56 ⊗ 56 of E7 = 1 + 133 + 1463 + 1539
        e7 = E_series(7)
        r = tensor_product(Irrep(e7, 0, 0, 0, 0, 0, 1, 0), Irrep(e7, 0, 0, 0, 0, 0, 1, 0))
        @test sort([dimension(k) for k in keys(r.components)]) == [1, 133, 1463, 1539]
        # Reducible ⊗ reducible agrees with bilinearity
        g = SU(3)
        A = Rep(Irrep(g, 1, 0)) + Rep(Irrep(g, 1, 1))
        B = Rep(Irrep(g, 0, 2)) + 2 * Rep(Irrep(g, 1, 0))
        AB = tensor_product(A, B)
        expected = Rep(g, Dict{Irrep, Int}())
        for (x, m) in A.components, (y, n) in B.components
            expected = expected + (m * n) * tensor_product(x, y)
        end
        @test AB.components == expected.components
    end

    @testset "Plethysms" begin
        # Dimensions from the hook-content formula
        for ((g, l), p) in [((SU(3), [1, 1]), [3, 1]), ((SO(8), [1, 0, 0, 0]), [2, 2, 1]),
                            ((G_series(2), [1, 0]), [5]), ((E_series(6), [1, 0, 0, 0, 0, 0]), [2, 1]),
                            ((Sp(6), [1, 0, 0]), [1, 1, 1, 1]), ((SU(2), [3]), [4, 2])]
            V = Irrep(g, l)
            R = plethysm(V, SymmetricIrrep(p))
            @test rep_dim(R) == schur_dim(dimension(V), p)
            @test all(m > 0 for m in values(R.components))
        end
        # Λ^k of the vector of SO(10)
        g = SO(10)
        @test antisymmetric_power(2, Irrep(g, 1, 0, 0, 0, 0)).components == Dict(Irrep(g, 0, 1, 0, 0, 0) => 1)
        @test antisymmetric_power(5, Irrep(g, 1, 0, 0, 0, 0)).components ==
              Dict(Irrep(g, 0, 0, 0, 2, 0) => 1, Irrep(g, 0, 0, 0, 0, 2) => 1)
        # Plethysm of a reducible representation: S²(V ⊕ W) = S²V ⊕ V⊗W ⊕ S²W
        g = SU(3)
        V = Irrep(g, 1, 0); W = Irrep(g, 1, 1)
        lhs = plethysm(Rep(V) + Rep(W), SymmetricIrrep([2]))
        rhs = symmetric_power(2, V) + tensor_product(V, W) + symmetric_power(2, W)
        @test lhs.components == rhs.components
        # Sum over all Schur functors of degree n recovers V^{⊗n}
        V = Irrep(G_series(2), 1, 0)
        total = Rep(V.algebra, Dict{Irrep, Int}())
        for p in all_partitions(3)
            total = total + dimension(SymmetricIrrep(p)) * plethysm(V, SymmetricIrrep(p))
        end
        @test total.components == tensor_product(tensor_product(V, V), Rep(V)).components
    end

    @testset "Symmetric group characters" begin
        # Column orthogonality: Σ_λ χ^λ(μ)² = |centraliser of μ|
        n = 7
        parts = all_partitions(n)
        for μ in parts
            s = sum(character(SymmetricIrrep(λ), ConjugacyClass(μ))^2 for λ in parts)
            @test s == factorial(n) ÷ conjugacy_class_size(ConjugacyClass(μ))
        end
    end

    @testset "cosh_product" begin
        # Kostant: adjoint ↦ 2^rank copies of V_ρ
        for g in [SU(2), SU(3), SU(5), SO(7), Sp(6), SO(8), SO(10), G_series(2), F_series(4), E_series(6)]
            R = cosh_product(adjoint_irrep(g))
            @test R.components == Dict(Irrep(g, ones(Int, lie_rank(g))) => 2^lie_rank(g))
        end
        # Vector of SO(N): the spinor(s), times 2 for odd N (one zero weight)
        @test cosh_product(Irrep(SO(7), 1, 0, 0)).components == Dict(Irrep(SO(7), 0, 0, 1) => 2)
        @test cosh_product(Irrep(SO(10), 1, 0, 0, 0, 0)).components ==
              Dict(Irrep(SO(10), 0, 0, 0, 1, 0) => 1, Irrep(SO(10), 0, 0, 0, 0, 1) => 1)
        # Principal SU(2) in SO(5): spinor 4 = [3]
        @test cosh_product(Irrep(SU(2), 4)).components == Dict(Irrep(SU(2), 3) => 2)
        # Product over all weights = exterior algebra
        for r in [Irrep(SU(3), 1, 1), Irrep(G_series(2), 1, 0), Irrep(SO(5), 0, 1)]
            L = Rep(Irrep(r.algebra, zeros(Int, lie_rank(r.algebra))))
            for k in 1:dimension(r)
                L = L + antisymmetric_power(k, r)
            end
            @test cosh_product(r; half = false).components == L.components
        end
        # χ² = 2^{m₀} ch Λ•(r), and dim χ = 2^{|Φ⁺(r)| + m₀}
        for r in [Irrep(SU(3), 2, 2), Irrep(G_series(2), 2, 0), Irrep(SO(7), 0, 0, 2)]
            m0 = character(r)[Weight(r.algebra, zeros(Int, lie_rank(r.algebra)))]
            R = cosh_product(r)
            @test rep_dim(R) == big(2)^((dimension(r) - m0) ÷ 2 + m0)
            L = cosh_product(r; half = false)
            @test tensor_product(R, R).components == Dict(k => v * 2^m0 for (k, v) in L.components)
        end
        # Multiplicities beyond 64 bits
        R = cosh_product(Irrep(SU(3), 4, 4))
        @test rep_dim(R) == big(2)^(60 + 5)    # 60 weight pairs, 5 zero weights
        # Every choice of Levi subalgebra gives the same answer
        for r in [Irrep(SU(4), 2, 0, 2), Irrep(SO(9), 0, 0, 1, 0), Irrep(F_series(4), 0, 0, 0, 1),
                  Irrep(G_series(2), 1, 1)]
            n = lie_rank(r.algebra)
            ref = cosh_product(r; levi = :none).components
            @test all(cosh_product(r; levi = setdiff(1:n, k)).components == ref for k in 1:n)
            @test cosh_product(r; levi = Int[]).components == ref
            @test cosh_product(r).components == ref
        end
        # Symmetric spaces: E8/Spin(16) needs a large orbit of 128 weights; check a smaller one,
        # E6/F4: the 26 of F4 gives a single irrep
        @test cosh_product(Irrep(F_series(4), 0, 0, 0, 1)).components ==
              Dict(Irrep(F_series(4), 0, 0, 1, 1) => 4)
        # Pseudo-real representations do not give an integral character
        @test_throws ErrorException cosh_product(Irrep(SU(2), 1))
        @test_throws ErrorException cosh_product(Irrep(SU(3), 1, 0))
    end

    @testset "Multi-modular arithmetic" begin
        T = Kiwi._exact_int_type(300)
        a = big(3)^100 - big(7)^50
        b = -(big(5)^50) + 12345
        x = T(a) * T(b) + T(b) - T(a)
        @test Kiwi._to_bigint(x) == a * b + b - a
        @test Kiwi._to_bigint(-T(a)) == -a
    end
end
