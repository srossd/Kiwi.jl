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
end
