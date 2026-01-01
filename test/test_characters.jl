using Test
using Kiwi

@testset "Characters" begin
    @testset "A Series Characters" begin
        # A_1 (SU(2)) spin-1/2
        g = A_series(1)
        rep = Irrep(g, 1)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 2
        @test length(char) == 2
        
        # A_1 (SU(2)) spin-1
        rep = Irrep(g, 2)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 3
        @test length(char) == 3
        
        # A_2 (SU(3)) fundamental
        g = A_series(2)
        rep = Irrep(g, 1, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 3
        @test length(char) == 3
        
        # A_2 (SU(3)) adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 8
        @test length(char) == 7  # 6 weights mult 1, 1 weight mult 2
        
        # A_3 [1,1,1]
        g = A_series(3)
        rep = Irrep(g, 1, 1, 1)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 64
        @test length(char) == 38
    end
    
    @testset "B Series Characters" begin
        # B_2 (SO(5)) vector
        g = B_series(2)
        rep = Irrep(g, 1, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 5
        
        # B_2 (SO(5)) adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 10
    end
    
    @testset "C Series Characters" begin
        # C_2 (Sp(4)) fundamental
        g = C_series(2)
        rep = Irrep(g, 1, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 4
        
        # C_2 (Sp(4)) adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 10
    end
    
    @testset "D Series Characters" begin
        # D_3 (SO(6)) vector
        g = D_series(3)
        rep = Irrep(g, 1, 0, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 6
        
        # D_4 (SO(8)) vector
        g = D_series(4)
        rep = Irrep(g, 1, 0, 0, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 8
    end
    
    @testset "G_2 Characters" begin
        # G_2 fundamental
        g = G_series(2)
        rep = Irrep(g, 1, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 7
        
        # G_2 adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 14
    end
    
    @testset "F_4 Characters" begin
        # F_4 fundamental
        g = F_series(4)
        rep = Irrep(g, 1, 0, 0, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 52
        
        # F_4 adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 52
    end
    
    @testset "E Series Characters" begin
        # E_6 fundamental
        g = E_series(6)
        rep = Irrep(g, 1, 0, 0, 0, 0, 0)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 27
        
        # E_6 adjoint
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 78
        
        # E_7 adjoint
        g = E_series(7)
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 133
        
        # E_8 adjoint
        g = E_series(8)
        rep = adjoint_irrep(g)
        char = character(rep)
        @test dimension(char) == dimension(rep)
        @test dimension(rep) == 248
    end
    
    @testset "Utility Functions" begin
        # Test total_weights and distinct_weights
        g = A_series(2)
        rep = adjoint_irrep(g)
        char = character(rep)
        
        @test dimension(char) == 8
        @test length(char) == 7
        @test dimension(char) == dimension(rep)
    end
    
    @testset "Weight Multiplicity Structure" begin
        # Verify that WeightMultiplicity objects are created correctly
        g = A_series(1)
        rep = Irrep(g, 1)
        char = character(rep)
        
        @test all(wm -> wm isa WeightMultiplicity, char)
        @test all(wm -> wm.multiplicity > 0, char)
        @test all(wm -> wm.weight isa Weight, char)
    end
    
    @testset "Lazy Character Evaluation" begin
        # Test lazy character evaluation against eager evaluation
        
        # Small representation: SU(3) fundamental
        g1 = A_series(2)
        rep1 = Irrep(g1, [1, 0])
        char1 = character(rep1)
        lazy1 = character(rep1; lazy=true)
        
        # Test that lazy and full give same multiplicities
        for wm in char1
            @test lazy1[wm.weight] == wm.multiplicity
        end
        
        # Test dimension
        @test dimension(lazy1) == dimension(char1)
        @test dimension(lazy1) == dimension(rep1)
        
        # Medium representation: SU(3) adjoint
        g2 = A_series(2)
        rep2 = adjoint_irrep(g2)
        char2 = character(rep2)
        lazy2 = character(rep2; lazy=true)
        
        # Test all weights match
        for wm in char2
            @test lazy2[wm.weight] == wm.multiplicity
        end
        
        # Test querying non-existent weight returns 0
        w_invalid = Weight(g2, [10, 10])
        @test lazy2[w_invalid] == 0
        @test char2[w_invalid] == 0
        
        # Larger representation: SU(5) [1,1,1,1]
        g3 = A_series(4)
        rep3 = Irrep(g3, [1, 1, 1, 1])
        char3 = character(rep3)
        lazy3 = character(rep3; lazy=true)
        
        # Test highest weight
        hw = highest_weight(rep3)
        @test lazy3[hw] == 1
        @test char3[hw] == 1
        
        # Test a sample of weights
        sample_weights = [
            Weight(g3, [1, 1, 1, 1]),  # highest weight
            Weight(g3, [0, 1, 1, 0]),
            Weight(g3, [0, 0, 0, 0]),
            Weight(g3, [1, 0, 1, 1]),
        ]
        
        for w in sample_weights
            @test lazy3[w] == char3[w]
        end
        
        @test dimension(lazy3) == dimension(char3)
        @test dimension(lazy3) == 1024
        
        # Test B_2 (SO(5))
        g4 = B_series(2)
        rep4 = Irrep(g4, [1, 0])
        char4 = character(rep4)
        lazy4 = character(rep4; lazy=true)
        
        for wm in char4
            @test lazy4[wm.weight] == wm.multiplicity
        end
        
        @test dimension(lazy4) == dimension(char4)
        @test dimension(lazy4) == 5
    end
end
