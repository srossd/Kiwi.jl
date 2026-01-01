using Test
using Kiwi

@testset "Tensor Products" begin
    @testset "Rep Type" begin
        # Test Rep construction from Irrep
        g = A_series(2)
        irrep = Irrep(g, 1, 0)
        rep = Rep(irrep)
        
        @test length(rep.components) == 1
        @test rep.components[irrep] == 1
        
        # Test Rep addition
        irrep2 = Irrep(g, 0, 1)
        rep2 = Rep(irrep2)
        sum_rep = rep + rep2
        @test length(sum_rep.components) == 2
        @test sum_rep.components[irrep] == 1
        @test sum_rep.components[irrep2] == 1
        
        # Test scalar multiplication
        doubled = 2 * rep
        @test doubled.components[irrep] == 2
    end
    
    @testset "A Series Tensor Products" begin
        g = A_series(2)
        
        # [1,0] ⊗ [1,0] = [2,0] ⊕ [0,1]
        fund = Irrep(g, 1, 0)
        result = tensor_product(fund, fund)
        @test length(result.components) == 2
        @test result.components[Irrep(g, 2, 0)] == 1
        @test result.components[Irrep(g, 0, 1)] == 1
        
        # Check dimension
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(fund)^2
        
        # [1,0] ⊗ [0,1]
        fund1 = Irrep(g, 1, 0)
        fund2 = Irrep(g, 0, 1)
        result = tensor_product(fund1, fund2)
        
        # Check dimension is correct
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(fund1) * dimension(fund2)
        
        # Adjoint ⊗ adjoint
        adj = adjoint_irrep(g)
        result = tensor_product(adj, adj)
        
        # Check dimension adds up
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(adj)^2
        
        # Should contain trivial rep
        @test haskey(result.components, Irrep(g, 0, 0))
        @test result.components[Irrep(g, 0, 0)] == 1
        
        # Should contain adjoint with multiplicity 2
        @test result.components[adj] == 2
    end
    
    @testset "B Series Tensor Products" begin
        g = B_series(2)
        
        # Vector ⊗ vector
        vec = Irrep(g, 1, 0)
        result = tensor_product(vec, vec)
        
        # Check dimension
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(vec)^2
        
        # Should decompose into symmetric and antisymmetric parts
        @test length(result.components) >= 2
    end
    
    @testset "C Series Tensor Products" begin
        g = C_series(2)
        
        # Fundamental ⊗ fundamental
        fund = Irrep(g, 1, 0)
        result = tensor_product(fund, fund)
        
        # Check dimension
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(fund)^2
        
        # For Sp(4), [1,0] ⊗ [1,0] should give symmetric and antisymmetric parts
        @test length(result.components) >= 2
    end
    
    @testset "G_2 Tensor Products" begin
        g = G_series(2)
        
        # 7-dimensional fundamental
        fund = Irrep(g, 1, 0)
        result = tensor_product(fund, fund)
        
        # Check dimension
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(fund)^2
    end
    
    @testset "Tensor Product with Trivial" begin
        g = A_series(2)
        
        # Check dimension is preserved
        rep = Irrep(g, 2, 1)
        trivial = Irrep(g, 0, 0)
        result = tensor_product(rep, trivial)
        
        total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
        @test total_dim == dimension(rep)
    end
    
    @testset "Symmetry" begin
        # Tensor product should be commutative
        g = A_series(2)
        rep1 = Irrep(g, 1, 0)
        rep2 = Irrep(g, 0, 1)
        
        result1 = tensor_product(rep1, rep2)
        result2 = tensor_product(rep2, rep1)
        
        # Should give the same decomposition
        @test result1.components == result2.components
    end
    
    @testset "Reducible Tensor Products" begin
        g = A_series(2)
        
        # Create a reducible rep: [1,0] ⊕ [0,1]
        irrep1 = Irrep(g, 1, 0)
        irrep2 = Irrep(g, 0, 1)
        rep = Rep(irrep1) + Rep(irrep2)
        
        # Tensor with fundamental
        fund = Irrep(g, 1, 0)
        result = tensor_product(rep, fund)
        
        # Should be the sum of individual tensor products
        expected = tensor_product(irrep1, fund) + tensor_product(irrep2, fund)
        @test result.components == expected.components
    end
    
    @testset "Dimension Consistency" begin
        # Test various tensor products have correct dimensions
        test_cases = [
            (A_series(2), [1, 0], [1, 0]),
            (A_series(2), [1, 1], [1, 1]),
            (A_series(3), [1, 0, 0], [1, 0, 0]),
            (B_series(2), [1, 0], [1, 0]),
            (C_series(2), [1, 0], [1, 0]),
            (G_series(2), [1, 0], [0, 1]),
        ]
        
        for (g, labels1, labels2) in test_cases
            rep1 = Irrep(g, labels1...)
            rep2 = Irrep(g, labels2...)
            result = tensor_product(rep1, rep2)
            
            total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
            expected_dim = dimension(rep1) * dimension(rep2)
            
            @test total_dim == expected_dim
        end
    end
    
    @testset "irreps Function" begin
        g = A_series(2)
        adj = adjoint_irrep(g)
        result = tensor_product(adj, adj)
        
        # Test irreps function returns correct format
        components = irreps(result)
        @test all(x -> x[1] isa Irrep && x[2] isa Int, components)
        
        # Check total count
        @test length(components) == length(result.components)
        
        # Check each component
        for (irrep, mult) in components
            @test result.components[irrep] == mult
        end
    end
end

@testset "Plethysms" begin
    @testset "Symmetric Powers" begin
        g = A_series(2)  # SU(3)
        fund = Irrep(g, 1, 0)
        
        # S²(3) = [2,0] (symmetric square)
        sym2 = symmetric_power(2, fund)
        @test length(sym2.components) == 1
        @test haskey(sym2.components, Irrep(g, 2, 0))
        @test dimension(sym2) == 6  # dim(3) * (dim(3) + 1) / 2
        
        # S³(3) = [3,0] (symmetric cube)
        sym3 = symmetric_power(3, fund)
        @test length(sym3.components) == 1
        @test haskey(sym3.components, Irrep(g, 3, 0))
        @test dimension(sym3) == 10  # binomial(3+3-1, 3)
        
        # S¹(3) = [1,0] (identity)
        sym1 = symmetric_power(1, fund)
        @test length(sym1.components) == 1
        @test haskey(sym1.components, fund)
        @test dimension(sym1) == 3
    end
    
    @testset "Antisymmetric Powers" begin
        g = A_series(2)  # SU(3)
        fund = Irrep(g, 1, 0)
        
        # Λ²(3) = [0,1] (antisymmetric square)
        alt2 = antisymmetric_power(2, fund)
        @test length(alt2.components) == 1
        @test haskey(alt2.components, Irrep(g, 0, 1))
        @test dimension(alt2) == 3  # dim(3) * (dim(3) - 1) / 2
        
        # Λ³(3) = [0,0] (top form, trivial)
        alt3 = antisymmetric_power(3, fund)
        @test length(alt3.components) == 1
        @test haskey(alt3.components, Irrep(g, 0, 0))
        @test dimension(alt3) == 1
        
        # Λ¹(3) = [1,0] (identity)
        alt1 = antisymmetric_power(1, fund)
        @test length(alt1.components) == 1
        @test haskey(alt1.components, fund)
        @test dimension(alt1) == 3
    end
    
    @testset "Symmetric Power of Adjoint" begin
        g = A_series(2)
        adj = Irrep(g, 1, 1)
        
        # S²(adj) should decompose into multiple irreps
        sym2 = symmetric_power(2, adj)
        @test dimension(sym2) == 36  # 8 * 9 / 2
        
        # Check it contains the trivial and adjoint
        @test haskey(sym2.components, Irrep(g, 0, 0))
        @test haskey(sym2.components, adj)
        
        # S³(adj) - more complex decomposition
        sym3 = symmetric_power(3, adj)
        @test dimension(sym3) == 120  # 8 * 9 * 10 / 6
        @test length(sym3.components) >= 5  # Should have multiple components
    end
    
    @testset "General Plethysm" begin
        g = A_series(2)
        fund = Irrep(g, 1, 0)
        
        # [2,1] representation of S₃ applied to fundamental
        rho = SymmetricIrrep([2, 1])
        result = plethysm(fund, rho)
        
        # This should give the adjoint [1,1]
        @test length(result.components) == 1
        @test haskey(result.components, Irrep(g, 1, 1))
        @test dimension(result) == 8
        
        # [1,1,1] is same as Λ³
        rho_alt3 = SymmetricIrrep([1, 1, 1])
        result = plethysm(fund, rho_alt3)
        alt3_direct = antisymmetric_power(3, fund)
        @test result.components == alt3_direct.components
    end
    
    @testset "Plethysm Dimension Consistency" begin
        # For SU(2), test several cases
        g = A_series(1)
        fund = Irrep(g, 1)
        
        # S²(2) should have dimension 3
        sym2 = symmetric_power(2, fund)
        @test dimension(sym2) == 3
        
        # Λ²(2) should have dimension 1
        alt2 = antisymmetric_power(2, fund)
        @test dimension(alt2) == 1
        
        # S³(2) should have dimension 4
        sym3 = symmetric_power(3, fund)
        @test dimension(sym3) == 4
    end
    
    @testset "Different Lie Algebras" begin
        # Test plethysms work for different Lie algebra types
        
        # B_2 (SO(5))
        g_b = B_series(2)
        vec_b = Irrep(g_b, 1, 0)
        sym2_b = symmetric_power(2, vec_b)
        @test dimension(sym2_b) == 15  # 5 * 6 / 2
        
        # C_2 (Sp(4))
        g_c = C_series(2)
        fund_c = Irrep(g_c, 1, 0)
        sym2_c = symmetric_power(2, fund_c)
        @test dimension(sym2_c) == 10  # 4 * 5 / 2
        
        # G_2
        g_g = G_series(2)
        fund_g = Irrep(g_g, 1, 0)
        sym2_g = symmetric_power(2, fund_g)
        @test dimension(sym2_g) == 28  # 7 * 8 / 2
    end
    
    @testset "Plethysm Properties" begin
        g = A_series(2)
        fund = Irrep(g, 1, 0)
        
        # S⁰ should give trivial (not implemented, but conceptually)
        # S¹ should be identity
        sym1 = symmetric_power(1, fund)
        @test sym1.components[fund] == 1
        
        # Λ⁰ should give trivial (conceptually, dimension 1)
        # Λ¹ should be identity
        alt1 = antisymmetric_power(1, fund)
        @test alt1.components[fund] == 1
        
        # Dimension formula: binomial(n+k-1, k) for Sᵏ(n)
        # For S²(3): binomial(3+2-1, 2) = binomial(4, 2) = 6
        sym2 = symmetric_power(2, fund)
        @test dimension(sym2) == binomial(3+2-1, 2)
        
        # Dimension formula: binomial(n, k) for Λᵏ(n)
        # For Λ²(3): binomial(3, 2) = 3
        alt2 = antisymmetric_power(2, fund)
        @test dimension(alt2) == binomial(3, 2)
    end
end
