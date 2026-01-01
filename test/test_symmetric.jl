using Test
using Kiwi

@testset "Symmetric Group" begin
    
    @testset "Partitions" begin
        # Basic construction
        p1 = Partition([3, 2, 1])
        @test p1.parts == [3, 2, 1]
        
        # Automatic sorting
        p2 = Partition([1, 3, 2])
        @test p2.parts == [3, 2, 1]
        @test p1 == p2
        
        # Remove zeros
        p3 = Partition([3, 0, 2, 0, 1])
        @test p3.parts == [3, 2, 1]
        
        # Convenience constructors
        p4 = Partition(3, 2, 1)
        @test p4 == p1
        
        p5 = Partition(5)
        @test p5.parts == [5]
        
        # Empty partition (edge case)
        p_empty = Partition(Int[])
        @test partition_size(p_empty) == 0
        
        # Partition size
        @test partition_size(p1) == 6
        @test partition_size(p5) == 5
    end
    
    @testset "Symmetric Irreps" begin
        # Basic construction
        irrep1 = SymmetricIrrep(Partition([3, 2]))
        @test partition_size(irrep1) == 5
        
        irrep2 = SymmetricIrrep(Partition([5]))
        irrep3 = SymmetricIrrep(Partition([1, 1, 1, 1, 1]))
        
        @test irrep1 != irrep2
        @test irrep2 != irrep3
    end
    
    @testset "Conjugacy Classes" begin
        # Basic construction
        cc1 = ConjugacyClass(Partition([3, 2, 1]))
        @test partition_size(cc1) == 6
        
        # Conjugacy class sizes in S₆
        # [6]: one 6-cycle, size = 6!/6 = 120
        @test conjugacy_class_size(ConjugacyClass(Partition([6]))) == 120
        
        # [1,1,1,1,1,1]: identity, size = 1
        @test conjugacy_class_size(ConjugacyClass(Partition([1,1,1,1,1,1]))) == 1
        
        # [2,2,2]: three 2-cycles, size = 6!/(2³·3!) = 15
        @test conjugacy_class_size(ConjugacyClass(Partition([2,2,2]))) == 15
        
        # [3,3]: two 3-cycles, size = 6!/(3²·2!) = 40
        @test conjugacy_class_size(ConjugacyClass(Partition([3,3]))) == 40
        
        # [4,2]: one 4-cycle and one 2-cycle, size = 6!/(4·2) = 90
        @test conjugacy_class_size(ConjugacyClass(Partition([4,2]))) == 90
    end
    
    @testset "Hook Lengths" begin
        # Partition [3,2,1]
        #  □ □ □
        #  □ □
        #  □
        p = Partition([3, 2, 1])
        
        # Hook lengths:
        # 5 3 1
        # 3 1
        # 1
        @test hook_length(p, 1, 1) == 5
        @test hook_length(p, 1, 2) == 3
        @test hook_length(p, 1, 3) == 1
        @test hook_length(p, 2, 1) == 3
        @test hook_length(p, 2, 2) == 1
        @test hook_length(p, 3, 1) == 1
        
        # Partition [4]
        p2 = Partition([4])
        @test hook_length(p2, 1, 1) == 4
        @test hook_length(p2, 1, 2) == 3
        @test hook_length(p2, 1, 3) == 2
        @test hook_length(p2, 1, 4) == 1
    end
    
    @testset "Dimensions via Hook Length Formula" begin
        # S₃ irreps
        # [3]: trivial rep, dim = 1
        @test dimension(SymmetricIrrep(Partition([3]))) == 1
        
        # [1,1,1]: sign rep, dim = 1
        @test dimension(SymmetricIrrep(Partition([1,1,1]))) == 1
        
        # [2,1]: standard rep, dim = 2
        @test dimension(SymmetricIrrep(Partition([2,1]))) == 2
        
        # S₄ irreps
        # [4]: dim = 1
        @test dimension(SymmetricIrrep(Partition([4]))) == 1
        
        # [3,1]: dim = 3
        @test dimension(SymmetricIrrep(Partition([3,1]))) == 3
        
        # [2,2]: dim = 2
        @test dimension(SymmetricIrrep(Partition([2,2]))) == 2
        
        # [2,1,1]: dim = 3
        @test dimension(SymmetricIrrep(Partition([2,1,1]))) == 3
        
        # [1,1,1,1]: dim = 1
        @test dimension(SymmetricIrrep(Partition([1,1,1,1]))) == 1
        
        # S₅ irreps
        # [5]: dim = 1
        @test dimension(SymmetricIrrep(Partition([5]))) == 1
        
        # [4,1]: dim = 4
        @test dimension(SymmetricIrrep(Partition([4,1]))) == 4
        
        # [3,2]: dim = 5
        @test dimension(SymmetricIrrep(Partition([3,2]))) == 5
        
        # [3,1,1]: dim = 6
        @test dimension(SymmetricIrrep(Partition([3,1,1]))) == 6
        
        # Check sum of squares equals |G|
        dims_S5 = [dimension(SymmetricIrrep(p)) for p in all_partitions(5)]
        @test sum(d^2 for d in dims_S5) == factorial(5)
    end
    
    @testset "All Partitions" begin
        # S₁
        @test length(all_partitions(1)) == 1
        @test Partition([1]) in all_partitions(1)
        
        # S₂
        ps2 = all_partitions(2)
        @test length(ps2) == 2
        @test Partition([2]) in ps2
        @test Partition([1,1]) in ps2
        
        # S₃
        ps3 = all_partitions(3)
        @test length(ps3) == 3
        @test Partition([3]) in ps3
        @test Partition([2,1]) in ps3
        @test Partition([1,1,1]) in ps3
        
        # S₄
        ps4 = all_partitions(4)
        @test length(ps4) == 5
        
        # S₅
        ps5 = all_partitions(5)
        @test length(ps5) == 7
    end
    
    @testset "Murnaghan-Nakayama Rule" begin
        # S₃ character table
        # Characters on identity [1,1,1]
        @test character(SymmetricIrrep(Partition([3])), ConjugacyClass(Partition([1,1,1]))) == 1
        @test character(SymmetricIrrep(Partition([2,1])), ConjugacyClass(Partition([1,1,1]))) == 2
        @test character(SymmetricIrrep(Partition([1,1,1])), ConjugacyClass(Partition([1,1,1]))) == 1
        
        # Characters on [3] (3-cycle)
        @test character(SymmetricIrrep(Partition([3])), ConjugacyClass(Partition([3]))) == 1
        @test character(SymmetricIrrep(Partition([2,1])), ConjugacyClass(Partition([3]))) == -1
        @test character(SymmetricIrrep(Partition([1,1,1])), ConjugacyClass(Partition([3]))) == 1
        
        # Characters on [2,1] (2-cycle)
        @test character(SymmetricIrrep(Partition([3])), ConjugacyClass(Partition([2,1]))) == 1
        @test character(SymmetricIrrep(Partition([2,1])), ConjugacyClass(Partition([2,1]))) == 0
        @test character(SymmetricIrrep(Partition([1,1,1])), ConjugacyClass(Partition([2,1]))) == -1
        
        # S₄ examples
        # [4,1] on identity [1,1,1,1,1]
        @test character(SymmetricIrrep(Partition([4,1])), ConjugacyClass(Partition([1,1,1,1,1]))) == 4
        
        # [3,2] on identity
        @test character(SymmetricIrrep(Partition([3,2])), ConjugacyClass(Partition([1,1,1,1,1]))) == 5
        
        # [2,2] on [2,2,1]
        @test character(SymmetricIrrep(Partition([2,2,1])), ConjugacyClass(Partition([2,2,1]))) == 1
    end
    
    @testset "Orthogonality Relations" begin
        # First orthogonality: Σ_g χ^λ(g) χ^μ(g) = |G| δ_λμ
        # For S₅
        n = 5
        partitions = all_partitions(n)
        
        for λ in partitions
            for μ in partitions
                irrep_λ = SymmetricIrrep(λ)
                irrep_μ = SymmetricIrrep(μ)
                
                # Sum over conjugacy classes weighted by class size
                total = 0
                for ρ in partitions
                    cc = ConjugacyClass(ρ)
                    χ_λ = character(irrep_λ, cc)
                    χ_μ = character(irrep_μ, cc)
                    class_size = conjugacy_class_size(cc)
                    total += χ_λ * χ_μ * class_size
                end
                
                if λ == μ
                    @test total == factorial(n)
                else
                    @test total == 0
                end
            end
        end
    end
    
    @testset "Character Table Properties" begin
        # Build full character table for S₄
        table = character_table(4)
        
        # Check that it has the right number of entries
        n_partitions = length(all_partitions(4))
        @test length(table) == n_partitions^2
        
        # Check dimensions sum correctly
        partitions = all_partitions(4)
        for λ in partitions
            irrep = SymmetricIrrep(λ)
            # χ(1) = dimension
            id_class = ConjugacyClass(Partition([1,1,1,1]))
            @test table[(irrep, id_class)] == dimension(irrep)
        end
    end
    
    @testset "SymmetricCharacter - Eager Evaluation" begin
        # Test eager evaluation (default)
        irrep = SymmetricIrrep([2, 1])
        char = SymmetricCharacter(irrep)
        
        # Should be precomputed (not lazy)
        @test !char.lazy
        @test length(char.values) == length(all_partitions(3))  # All conjugacy classes of S₃
        
        # Test indexing by ConjugacyClass
        @test char[ConjugacyClass([3])] == -1
        @test char[ConjugacyClass([2, 1])] == -0
        @test char[ConjugacyClass([1, 1, 1])] == 2
        
        # Test indexing by Vector{Int}
        @test char[[3]] == -1
        @test char[[2, 1]] == 0
        @test char[[1, 1, 1]] == 2
        
        # Test dimension
        @test dimension(char) == 2
        
        # Character values should match direct computation
        for cc_partition in all_partitions(3)
            cc = ConjugacyClass(cc_partition)
            @test char[cc] == character(irrep, cc)
        end
    end
    
    @testset "SymmetricCharacter - Lazy Evaluation" begin
        # Test lazy evaluation
        irrep = SymmetricIrrep([2, 2])
        char = SymmetricCharacter(irrep; lazy=true)
        
        # Should be lazy
        @test char.lazy
        @test length(char.values) == 0  # Nothing computed yet
        
        # Access a value - should compute and cache
        val1 = char[[4]]
        @test val1 == 0
        @test length(char.values) == 1  # One value now cached
        
        # Access same value again - should be cached
        val2 = char[ConjugacyClass([4])]
        @test val2 == val1
        @test length(char.values) == 1  # Still just one value
        
        # Access different value
        val3 = char[[2, 2]]
        @test val3 == 2
        @test length(char.values) == 2  # Two values now
        
        # Test dimension (doesn't require computing all values)
        @test dimension(char) == 2
        
        # All computed values should match direct computation
        for (cc, value) in char.values
            @test value == character(irrep, cc)
        end
    end
    
    @testset "SymmetricCharacter - Consistency" begin
        # Eager and lazy should give same results
        irrep = SymmetricIrrep([3, 2, 1])
        char_eager = SymmetricCharacter(irrep)
        char_lazy = SymmetricCharacter(irrep; lazy=true)
        
        # Test on all conjugacy classes
        for cc_partition in all_partitions(6)
            cc = ConjugacyClass(cc_partition)
            @test char_eager[cc] == char_lazy[cc]
        end
        
        # After lazy evaluation, should have computed all values
        @test length(char_lazy.values) == length(all_partitions(6))
    end
    
    @testset "SymmetricCharacter - S₄ Examples" begin
        # Standard representation [3, 1]
        irrep_std = SymmetricIrrep([3, 1])
        char_std = SymmetricCharacter(irrep_std)
        
        @test dimension(char_std) == 3
        @test char_std[[1,1,1,1]] == 3  # identity: dimension
        @test char_std[[4]] == -1        # 4-cycle
        @test char_std[[2,2]] == -1      # product of two 2-cycles
        @test char_std[[3,1]] == 0       # 3-cycle times 1-cycle
        @test char_std[[2,1,1]] == 1     # 2-cycle times identity
        
        # Sign representation [1,1,1,1]
        irrep_sign = SymmetricIrrep([1,1,1,1])
        char_sign = SymmetricCharacter(irrep_sign)
        
        @test dimension(char_sign) == 1
        @test char_sign[[1,1,1,1]] == 1
        @test char_sign[[4]] == -1
        @test char_sign[[2,2]] == 1
        @test char_sign[[3,1]] == 1
        @test char_sign[[2,1,1]] == -1
    end
    
end

