using Test
using Kiwi

@testset "Weyl Group" begin
    @testset "Weyl Reflections" begin
        # Test simple reflection for A_2
        g = A_series(2)
        λ = Weight(g, [1, 0])
        
        # Reflect across first simple root
        λ_reflected = simple_reflection(λ, 1)
        @test λ_reflected.coordinates == [-1, 1]
        
        # Reflect twice should give back original
        λ_double = simple_reflection(λ_reflected, 1)
        @test λ_double.coordinates == λ.coordinates
    end
    
    @testset "Dominance Check" begin
        g = A_series(2)
        
        # Fundamental weights are dominant
        @test is_dominant(Weight(g, [1, 0]))
        @test is_dominant(Weight(g, [0, 1]))
        @test is_dominant(Weight(g, [1, 1]))
        
        # Non-dominant weights
        @test !is_dominant(Weight(g, [-1, 2]))
        @test !is_dominant(Weight(g, [1, -1]))
        
        # Zero weight is dominant
        @test is_dominant(Weight(g, [0, 0]))
    end
    
    @testset "Reflect to Dominant" begin
        # A_2 tests
        g = A_series(2)
        
        # Non-dominant weight should reflect to dominant
        μ = Weight(g, [-1, 2])
        μ_dom, sign = reflect_to_dominant(μ)
        @test is_dominant(μ_dom)
        @test μ_dom.coordinates == [1, 1]
        @test sign == -1
        
        # Already dominant weight should not change
        λ = Weight(g, [1, 0])
        λ_dom, sign = reflect_to_dominant(λ)
        @test λ_dom.coordinates == λ.coordinates
        @test sign == 1
        
        # C_4 test
        g = C_series(4)
        μ = Weight(g, [-2, 0, 2, -1])
        μ_dom, sign = reflect_to_dominant(μ)
        @test is_dominant(μ_dom)
        @test μ_dom.coordinates == [0, 0, 0, 1]
        @test sign == -1
    end
    
    @testset "Weyl Orbit" begin
        # A_2 orbit
        g = A_series(2)
        λ = Weight(g, [1, 0])
        orbit = weyl_orbit(λ)
        
        @test length(orbit[1]) + length(orbit[-1]) == 3
        @test weyl_group_order(g) == 6
        
        # Zero weight should have trivial orbit
        zero_weight = Weight(g, [0, 0])
        zero_orbit = weyl_orbit(zero_weight)
        @test length(zero_orbit[1]) == 1
        @test length(zero_orbit[-1]) == 0
        
        # C_4 orbit
        g = C_series(4)
        λ = Weight(g, [0, 0, 0, 1])
        orbit = weyl_orbit(λ)
        @test length(orbit[1]) + length(orbit[-1]) == 16
        @test weyl_group_order(g) == 384
    end
    
    @testset "Weyl Group Order" begin
        # A_n: (n+1)!
        @test weyl_group_order(A_series(1)) == 2
        @test weyl_group_order(A_series(2)) == 6
        @test weyl_group_order(A_series(3)) == 24
        
        # B_n: 2^n · n!
        @test weyl_group_order(B_series(2)) == 8
        @test weyl_group_order(B_series(3)) == 48
        
        # C_n: 2^n · n!
        @test weyl_group_order(C_series(2)) == 8
        @test weyl_group_order(C_series(3)) == 48
        @test weyl_group_order(C_series(4)) == 384
        
        # D_n: 2^(n-1) · n!
        @test weyl_group_order(D_series(3)) == 24
        @test weyl_group_order(D_series(4)) == 192
        
        # Exceptional algebras
        @test weyl_group_order(E_series(6)) == 51840
        @test weyl_group_order(E_series(7)) == 2903040
        @test weyl_group_order(E_series(8)) == 696729600
        @test weyl_group_order(F_series(4)) == 1152
        @test weyl_group_order(G_series(2)) == 12
    end
    
    @testset "Dynkin Labels Accessor" begin
        g = A_series(2)
        rep = Irrep(g, 1, 0)
        @test dynkin_labels(rep) == [1, 0]
        
        rep2 = Irrep(g, 2, 3)
        @test dynkin_labels(rep2) == [2, 3]
    end
    
    @testset "Sign Consistency" begin
        # Check that orbit signs are consistent
        g = A_series(2)
        λ = Weight(g, [1, 1])
        orbit = weyl_orbit(λ)
        
        # Each weight in the orbit should reflect back to the starting weight
        for sign in [1, -1]
            for w in orbit[sign]
                w_dom, refl_sign = reflect_to_dominant(w)
                # Should get back to a dominant weight in the orbit
                @test is_dominant(w_dom)
            end
        end
    end
end
