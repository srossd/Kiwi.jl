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
