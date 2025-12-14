using Test
using Kiwi

@testset "Root Systems" begin
    
    @testset "Cartan Matrices" begin
        # A series
        @test cartan_matrix(A_series(2)) == [2 -1; -1 2]
        @test cartan_matrix(A_series(3)) == [2 -1 0; -1 2 -1; 0 -1 2]
        
        # B series
        @test cartan_matrix(B_series(2)) == [2 -2; -1 2]
        @test cartan_matrix(B_series(3)) == [2 -1 0; -1 2 -2; 0 -1 2]
        
        # C series
        @test cartan_matrix(C_series(2)) == [2 -1; -2 2]
        @test cartan_matrix(C_series(3)) == [2 -1 0; -1 2 -1; 0 -2 2]
        
        # D series
        @test cartan_matrix(D_series(3)) == [2 -1 -1; -1 2 0; -1 0 2]
        @test cartan_matrix(D_series(4)) == [2 -1 0 0; -1 2 -1 -1; 0 -1 2 0; 0 -1 0 2]
        
        # G series
        @test cartan_matrix(G_series(2)) == [2 -1; -3 2]
    end
    
    @testset "Root Squared Lengths" begin
        # Simply-laced: all roots have same length
        @test Kiwi.simple_root_squared_lengths(A_series(2)) == [2, 2]
        @test Kiwi.simple_root_squared_lengths(A_series(3)) == [2, 2, 2]
        @test Kiwi.simple_root_squared_lengths(D_series(4)) == [2, 2, 2, 2]
        @test Kiwi.simple_root_squared_lengths(E_series(6)) == fill(Rational{Int}(2), 6)
        
        # Non-simply-laced: different root lengths
        @test Kiwi.simple_root_squared_lengths(B_series(2)) == [2, 1]
        @test Kiwi.simple_root_squared_lengths(B_series(3)) == [2, 2, 1]
        @test Kiwi.simple_root_squared_lengths(C_series(2)) == [1, 2]
        @test Kiwi.simple_root_squared_lengths(C_series(3)) == [1, 1, 2]
        @test Kiwi.simple_root_squared_lengths(G_series(2)) == [2//3, 2]
        @test Kiwi.simple_root_squared_lengths(F_series(4)) == [2, 2, 1, 1]
    end
    
    @testset "Positive Root Counts" begin
        # A series: n(n+1)/2 roots
        @test length(Kiwi.positive_roots(A_series(2))) == 3
        @test length(Kiwi.positive_roots(A_series(3))) == 6
        @test length(Kiwi.positive_roots(A_series(4))) == 10
        @test length(Kiwi.positive_roots(A_series(5))) == 15
        
        # B series: n^2 roots
        @test length(Kiwi.positive_roots(B_series(2))) == 4
        @test length(Kiwi.positive_roots(B_series(3))) == 9
        @test length(Kiwi.positive_roots(B_series(4))) == 16
        @test length(Kiwi.positive_roots(B_series(5))) == 25
        
        # C series: n^2 roots
        @test length(Kiwi.positive_roots(C_series(2))) == 4
        @test length(Kiwi.positive_roots(C_series(3))) == 9
        @test length(Kiwi.positive_roots(C_series(4))) == 16
        @test length(Kiwi.positive_roots(C_series(5))) == 25
        
        # D series: n(n-1) roots
        @test length(Kiwi.positive_roots(D_series(3))) == 6
        @test length(Kiwi.positive_roots(D_series(4))) == 12
        @test length(Kiwi.positive_roots(D_series(5))) == 20
        @test length(Kiwi.positive_roots(D_series(6))) == 30
        
        # Exceptional algebras
        @test length(Kiwi.positive_roots(G_series(2))) == 6
        @test length(Kiwi.positive_roots(F_series(4))) == 24
        @test length(Kiwi.positive_roots(E_series(6))) == 36
        @test length(Kiwi.positive_roots(E_series(7))) == 63
        @test length(Kiwi.positive_roots(E_series(8))) == 120
    end
    
    @testset "Inner Products" begin
        # Test B_2 inner product
        g = B_series(2)
        w1 = Kiwi.Weight(g, [1, 2])
        w2 = Kiwi.Weight(g, [3, 4])
        @test Kiwi.inner_product(w1, w2) == 12
        
        # Test G_2 inner product
        g = G_series(2)
        w1 = Kiwi.Weight(g, [1, 2])
        w2 = Kiwi.Weight(g, [3, 4])
        @test Kiwi.inner_product(w1, w2) == 28
        
        # Test that inner products are symmetric
        g = A_series(3)
        w1 = Kiwi.Weight(g, [1, 2, 3])
        w2 = Kiwi.Weight(g, [4, 5, 6])
        @test Kiwi.inner_product(w1, w2) == Kiwi.inner_product(w2, w1)
    end
    
    @testset "Cartan Matrix Verification" begin
        # For any algebra, the inner product should reproduce the Cartan matrix
        # C_ij = 2(α_i, α_j) / (α_j, α_j)
        for g in [A_series(2), B_series(2), C_series(2), G_series(2)]
            C = cartan_matrix(g)
            sr = Kiwi.simple_roots(g)
            lengths_sq = Kiwi.simple_root_squared_lengths(g)
            
            for i in 1:g.rank
                for j in 1:g.rank
                    expected = C[i, j]
                    computed = 2 * Kiwi.inner_product(sr[i], sr[j]) / lengths_sq[j]
                    @test computed == expected
                end
            end
        end
    end
end
