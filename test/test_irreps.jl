using Test
using Kiwi

@testset "Irreducible Representations" begin
    
    @testset "A Series Dimensions" begin
        # SU(2) = A_1
        @test dimension(Irrep(A_series(1), 0)) == 1  # trivial
        @test dimension(Irrep(A_series(1), 1)) == 2  # fundamental
        @test dimension(Irrep(A_series(1), 2)) == 3  # spin-1
        @test dimension(Irrep(A_series(1), 3)) == 4  # spin-3/2
        
        # SU(3) = A_2
        @test dimension(Irrep(A_series(2), 0, 0)) == 1   # trivial
        @test dimension(Irrep(A_series(2), 1, 0)) == 3   # fundamental
        @test dimension(Irrep(A_series(2), 0, 1)) == 3   # antifundamental
        @test dimension(Irrep(A_series(2), 1, 1)) == 8   # adjoint
        @test dimension(Irrep(A_series(2), 2, 0)) == 6   # symmetric square
        @test dimension(Irrep(A_series(2), 0, 2)) == 6   # antisymmetric square
        @test dimension(Irrep(A_series(2), 3, 0)) == 10
        @test dimension(Irrep(A_series(2), 2, 1)) == 15
        
        # SU(4) = A_3
        @test dimension(Irrep(A_series(3), 1, 0, 0)) == 4   # fundamental
        @test dimension(Irrep(A_series(3), 0, 1, 0)) == 6   # 2-index antisym
        @test dimension(Irrep(A_series(3), 0, 0, 1)) == 4   # antifundamental
        @test dimension(Irrep(A_series(3), 1, 0, 1)) == 15  # adjoint
    end
    
    @testset "B Series Dimensions" begin
        # SO(5) = B_2
        @test dimension(Irrep(B_series(2), 0, 0)) == 1   # trivial
        @test dimension(Irrep(B_series(2), 1, 0)) == 5   # vector
        @test dimension(Irrep(B_series(2), 0, 1)) == 4   # spinor
        @test dimension(Irrep(B_series(2), 2, 0)) == 14  # symmetric traceless 2-tensor
        @test dimension(Irrep(B_series(2), 1, 1)) == 16
        
        # SO(7) = B_3
        @test dimension(Irrep(B_series(3), 1, 0, 0)) == 7   # vector
        @test dimension(Irrep(B_series(3), 0, 0, 1)) == 8   # spinor
        @test dimension(Irrep(B_series(3), 0, 1, 0)) == 21  # 2-index antisym
    end
    
    @testset "C Series Dimensions" begin
        # Sp(4) = C_2
        @test dimension(Irrep(C_series(2), 0, 0)) == 1   # trivial
        @test dimension(Irrep(C_series(2), 1, 0)) == 4   # fundamental
        @test dimension(Irrep(C_series(2), 0, 1)) == 5   # adjoint component
        @test dimension(Irrep(C_series(2), 2, 0)) == 10
        @test dimension(Irrep(C_series(2), 1, 1)) == 16
        
        # Sp(6) = C_3
        @test dimension(Irrep(C_series(3), 1, 0, 0)) == 6   # fundamental
        @test dimension(Irrep(C_series(3), 0, 1, 0)) == 14
        @test dimension(Irrep(C_series(3), 0, 0, 1)) == 14
    end
    
    @testset "D Series Dimensions" begin
        # SO(6) = D_3 = SU(4)
        @test dimension(Irrep(D_series(3), 1, 0, 0)) == 6  # vector (same as A_3 [0,1,0])
        @test dimension(Irrep(D_series(3), 0, 1, 0)) == 4  # spinor+
        @test dimension(Irrep(D_series(3), 0, 0, 1)) == 4  # spinor-
        
        # SO(8) = D_4 (triality)
        @test dimension(Irrep(D_series(4), 0, 0, 0, 0)) == 1   # trivial
        @test dimension(Irrep(D_series(4), 1, 0, 0, 0)) == 8   # vector
        @test dimension(Irrep(D_series(4), 0, 0, 1, 0)) == 8   # spinor+
        @test dimension(Irrep(D_series(4), 0, 0, 0, 1)) == 8   # spinor-
        @test dimension(Irrep(D_series(4), 0, 1, 0, 0)) == 28  # 2-index antisym
        @test dimension(Irrep(D_series(4), 2, 0, 0, 0)) == 35  # symmetric traceless
        
        # SO(10) = D_5
        @test dimension(Irrep(D_series(5), 1, 0, 0, 0, 0)) == 10  # vector
        @test dimension(Irrep(D_series(5), 0, 0, 0, 1, 0)) == 16  # spinor+
        @test dimension(Irrep(D_series(5), 0, 0, 0, 0, 1)) == 16  # spinor-
    end
    
    @testset "G_2 Dimensions" begin
        @test dimension(Irrep(G_series(2), 0, 0)) == 1   # trivial
        @test dimension(Irrep(G_series(2), 1, 0)) == 7   # fundamental
        @test dimension(Irrep(G_series(2), 0, 1)) == 14  # adjoint
        @test dimension(Irrep(G_series(2), 2, 0)) == 27
        @test dimension(Irrep(G_series(2), 1, 1)) == 64
        @test dimension(Irrep(G_series(2), 3, 0)) == 77
        @test dimension(Irrep(G_series(2), 0, 2)) == 77
    end
    
    @testset "F_4 Dimensions" begin
        @test dimension(Irrep(F_series(4), 0, 0, 0, 0)) == 1    # trivial
        @test dimension(Irrep(F_series(4), 1, 0, 0, 0)) == 52   # fundamental
        @test dimension(Irrep(F_series(4), 0, 0, 0, 1)) == 26
        @test dimension(Irrep(F_series(4), 0, 1, 0, 0)) == 1274 # adjoint
    end
    
    @testset "E Series Dimensions" begin
        # E_6
        @test dimension(Irrep(E_series(6), 0, 0, 0, 0, 0, 0)) == 1   # trivial
        @test dimension(Irrep(E_series(6), 1, 0, 0, 0, 0, 0)) == 27  # fundamental
        @test dimension(Irrep(E_series(6), 0, 0, 0, 0, 0, 1)) == 78  # adjoint
        @test dimension(Irrep(E_series(6), 0, 0, 0, 0, 1, 0)) == 27  # conjugate fund
        
        # E_7
        @test dimension(Irrep(E_series(7), 0, 0, 0, 0, 0, 0, 0)) == 1   # trivial
        @test dimension(Irrep(E_series(7), 1, 0, 0, 0, 0, 0, 0)) == 133 # adjoint (was 56 but reversed)
        @test dimension(Irrep(E_series(7), 0, 0, 0, 0, 0, 1, 0)) == 56  # fundamental (was 133 but reversed)
        
        # E_8
        @test dimension(Irrep(E_series(8), 0, 0, 0, 0, 0, 0, 0, 0)) == 1    # trivial
        @test dimension(Irrep(E_series(8), 0, 0, 0, 0, 0, 0, 1, 0)) == 248  # adjoint
        @test dimension(Irrep(E_series(8), 1, 0, 0, 0, 0, 0, 0, 0)) == 3875
    end
    
    @testset "Adjoint Representations" begin
        # Adjoint dimension = number of generators = rank + number of positive roots
        @test dimension(adjoint_irrep(A_series(2))) == 8    # SU(3)
        @test dimension(adjoint_irrep(A_series(3))) == 15   # SU(4)
        @test dimension(adjoint_irrep(B_series(2))) == 10   # SO(5)
        @test dimension(adjoint_irrep(B_series(3))) == 21   # SO(7)
        @test dimension(adjoint_irrep(C_series(2))) == 10   # Sp(4)
        @test dimension(adjoint_irrep(C_series(3))) == 21   # Sp(6)
        @test dimension(adjoint_irrep(D_series(3))) == 15   # SO(6)
        @test dimension(adjoint_irrep(D_series(4))) == 28   # SO(8)
        @test dimension(adjoint_irrep(G_series(2))) == 14
        @test dimension(adjoint_irrep(F_series(4))) == 52
        @test dimension(adjoint_irrep(E_series(6))) == 78
        @test dimension(adjoint_irrep(E_series(7))) == 133
        @test dimension(adjoint_irrep(E_series(8))) == 248
    end
    
    @testset "Quadratic Casimir" begin
        # A_2 (SU(3)) tests
        @test quadratic_casimir(adjoint_irrep(A_series(2))) == 3
        @test quadratic_casimir(Irrep(A_series(2), 1, 0)) == 4//3
        @test quadratic_casimir(Irrep(A_series(2), 0, 1)) == 4//3
        
        # Verify E_6 representation
        @test quadratic_casimir(Irrep(E_series(6), 1, 3, 0, 0, 2, 1)) == 368//3
        
        # Verify some other cases
        @test quadratic_casimir(Irrep(A_series(1), 1)) == 3//4  # SU(2) fundamental
        @test quadratic_casimir(adjoint_irrep(G_series(2))) == 4
    end
    
    @testset "Dynkin Index" begin
        # A_2 (SU(3)) tests
        @test dynkin_index(Irrep(A_series(2), 1, 0)) == 1//2
        @test dynkin_index(adjoint_irrep(A_series(2))) == 3
        @test dynkin_index(adjoint_irrep(A_series(3))) == 4
        @test dynkin_index(adjoint_irrep(G_series(2))) == 4
        
        # Verify E_6 representation
        @test dynkin_index(Irrep(E_series(6), 1, 3, 0, 0, 2, 1)) == 57164267160
    end
    
    @testset "Congruency Classes" begin
        # A_3 (SU(4)) - Z_4 center
        su4 = A_series(3)
        @test congruency_class(Irrep(su4, 0, 0, 0)).value == 0  # trivial
        @test congruency_class(Irrep(su4, 1, 0, 0)).value == 1  # fundamental
        @test congruency_class(Irrep(su4, 0, 1, 0)).value == 2  
        @test congruency_class(Irrep(su4, 0, 0, 1)).value == 3  # antifundamental
        @test congruency_class(Irrep(su4, 1, 0, 1)).value == 0  # adjoint (self-conjugate)
        # Three dimension-20 irreps with different congruency classes
        @test congruency_class(Irrep(su4, 0, 1, 1)).value == 1  # 20
        @test congruency_class(Irrep(su4, 1, 1, 0)).value == 3  # 20-bar
        @test congruency_class(Irrep(su4, 0, 2, 0)).value == 0  # 20-prime
        
        # B_3 (SO(7)) - Z_2 center
        so7 = B_series(3)
        @test congruency_class(Irrep(so7, 1, 0, 0)).value == 0  # vector
        @test congruency_class(Irrep(so7, 0, 1, 0)).value == 0  
        @test congruency_class(Irrep(so7, 0, 0, 1)).value == 1  # spinor
        
        # C_3 (Sp(6)) - Z_2 center
        sp6 = C_series(3)
        @test congruency_class(Irrep(sp6, 1, 0, 0)).value == 1  
        @test congruency_class(Irrep(sp6, 0, 1, 0)).value == 0  
        @test congruency_class(Irrep(sp6, 0, 0, 1)).value == 1  
        @test congruency_class(Irrep(sp6, 1, 0, 1)).value == 0  
        
        # D_4 (SO(8)) - Z_2 × Z_2 center (triality)
        so8 = D_series(4)
        @test congruency_class(Irrep(so8, 1, 0, 0, 0)).value == (0, 1)  # vector
        @test congruency_class(Irrep(so8, 0, 1, 0, 0)).value == (0, 0)  # adjoint
        @test congruency_class(Irrep(so8, 0, 0, 1, 0)).value == (1, 1)  # spinor_s
        @test congruency_class(Irrep(so8, 0, 0, 0, 1)).value == (1, 0)  # spinor_c
        
        # D_5 (SO(10)) - Z_4 center
        so10 = D_series(5)
        @test congruency_class(Irrep(so10, 1, 0, 0, 0, 0)).value == 2  # vector
        @test congruency_class(Irrep(so10, 0, 0, 0, 1, 0)).value == 3  
        @test congruency_class(Irrep(so10, 0, 0, 0, 0, 1)).value == 1  # spinor
        
        # E_6 - Z_3 center
        e6 = E_series(6)
        @test congruency_class(Irrep(e6, 1, 0, 0, 0, 0, 0)).value == 1
        @test congruency_class(Irrep(e6, 0, 0, 0, 0, 0, 1)).value == 0  # adjoint
        @test congruency_class(Irrep(e6, 0, 1, 0, 0, 0, 0)).value == 2
        
        # E_7 - Z_2 center
        e7 = E_series(7)
        @test congruency_class(Irrep(e7, 1, 0, 0, 0, 0, 0, 0)).value == 0  # adjoint
        @test congruency_class(Irrep(e7, 0, 0, 0, 1, 0, 0, 0)).value == 1
        
        # E_8 - trivial center
        e8 = E_series(8)
        @test congruency_class(Irrep(e8, 1, 0, 0, 0, 0, 0, 0, 0)).value === nothing
        
        # G_2 - trivial center
        g2 = G_series(2)
        @test congruency_class(Irrep(g2, 1, 0)).value === nothing
        
        # F_4 - trivial center
        f4 = F_series(4)
        @test congruency_class(Irrep(f4, 1, 0, 0, 0)).value === nothing
    end
    
    @testset "Finding Irreps by Dimension" begin
        # SU(3) - small test
        su3 = A_series(2)
        irreps_10 = irreps_up_to_dim(su3, 10)
        @test length(irreps_10) == 8
        @test all(dimension(rep) <= 10 for rep in irreps_10)
        @test Irrep(su3, 0, 0) in irreps_10  # trivial
        @test Irrep(su3, 1, 0) in irreps_10  # fundamental
        @test Irrep(su3, 0, 1) in irreps_10  # antifundamental
        @test Irrep(su3, 1, 1) in irreps_10  # adjoint
        @test Irrep(su3, 3, 0) in irreps_10  # dim 10
        @test Irrep(su3, 0, 3) in irreps_10  # dim 10
        
        # SU(4) - test caching
        su4 = A_series(3)
        irreps_50 = irreps_up_to_dim(su4, 50)
        @test all(dimension(rep) <= 50 for rep in irreps_50)
        
        # This should use cache
        irreps_20 = irreps_up_to_dim(su4, 20)
        @test all(dimension(rep) <= 20 for rep in irreps_20)
        @test length(irreps_20) < length(irreps_50)
        
        # Verify all irreps in irreps_20 are in irreps_50
        @test all(rep in irreps_50 for rep in irreps_20)
        
        # SO(8) - test with different algebra
        so8 = D_series(4)
        irreps_so8 = irreps_up_to_dim(so8, 30)
        @test all(dimension(rep) <= 30 for rep in irreps_so8)
        @test Irrep(so8, 0, 0, 0, 0) in irreps_so8  # trivial
        @test Irrep(so8, 1, 0, 0, 0) in irreps_so8  # vector (dim 8)
        @test Irrep(so8, 0, 0, 1, 0) in irreps_so8  # spinor_s (dim 8)
        @test Irrep(so8, 0, 0, 0, 1) in irreps_so8  # spinor_c (dim 8)
        @test Irrep(so8, 0, 1, 0, 0) in irreps_so8  # adjoint (dim 28)
        
        # Test with maxDim = 1 (only trivial)
        g2 = G_series(2)
        irreps_1 = irreps_up_to_dim(g2, 1)
        @test length(irreps_1) == 1
        @test irreps_1[1] == Irrep(g2, 0, 0)
        
        # Test cache clearing
        clear_irreps_cache!()
        irreps_20_again = irreps_up_to_dim(su4, 20)
        @test length(irreps_20_again) == length(irreps_20)
    end
    
    @testset "Dimension Labels" begin
        # SU(3) tests
        su3 = A_series(2)
        @test dimension_label(Irrep(su3, 0, 0)) == "1"
        @test dimension_label(Irrep(su3, 1, 0)) == "3"
        @test dimension_label(Irrep(su3, 0, 1)) == "3b"
        @test dimension_label(Irrep(su3, 2, 0)) == "6b"
        @test dimension_label(Irrep(su3, 0, 2)) == "6"
        @test dimension_label(Irrep(su3, 1, 1)) == "8"
        @test dimension_label(Irrep(su3, 3, 0)) == "10"
        @test dimension_label(Irrep(su3, 0, 3)) == "10b"
        
        # SU(4) tests - dimension 20 irreps
        su4 = A_series(3)
        @test dimension_label(Irrep(su4, 0, 1, 1)) == "20"
        @test dimension_label(Irrep(su4, 1, 1, 0)) == "20b"
        @test dimension_label(Irrep(su4, 0, 2, 0)) == "20′"
        @test dimension_label(Irrep(su4, 3, 0, 0)) == "20b′′"
        @test dimension_label(Irrep(su4, 0, 0, 3)) == "20′′"
        
        # SU(4) fundamental and antifundamental
        @test dimension_label(Irrep(su4, 1, 0, 0)) == "4"
        @test dimension_label(Irrep(su4, 0, 0, 1)) == "4b"
        @test dimension_label(Irrep(su4, 1, 0, 1)) == "15"  # adjoint
        
        # SO(8) triality tests
        so8 = D_series(4)
        @test dimension_label(Irrep(so8, 0, 0, 0, 0)) == "1"
        @test dimension_label(Irrep(so8, 1, 0, 0, 0)) == "8_v"
        @test dimension_label(Irrep(so8, 0, 0, 1, 0)) == "8_c"
        @test dimension_label(Irrep(so8, 0, 0, 0, 1)) == "8_s"
        @test dimension_label(Irrep(so8, 0, 1, 0, 0)) == "28"
        @test dimension_label(Irrep(so8, 2, 0, 0, 0)) == "35_v"
        @test dimension_label(Irrep(so8, 0, 0, 2, 0)) == "35_c"
        @test dimension_label(Irrep(so8, 0, 0, 0, 2)) == "35_s"
        @test dimension_label(Irrep(so8, 1, 0, 1, 0)) == "56_s"
        @test dimension_label(Irrep(so8, 1, 0, 0, 1)) == "56_c"
        @test dimension_label(Irrep(so8, 0, 0, 1, 1)) == "56_v"
    end
end
