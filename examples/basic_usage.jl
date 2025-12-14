# Example usage of Kiwi.jl - showcasing representation theory capabilities

using Kiwi

println("=== Kiwi.jl: Lie Algebra Representations ===\n")

# ============================================================================
# Example 1: Creating Lie Algebras and Irreducible Representations
# ============================================================================
println("Example 1: Creating Lie Algebras and Irreps")
println("-" ^ 60)

# Create classical and exceptional Lie algebras
su3 = A_series(2)     # SU(3)
so5 = B_series(2)     # SO(5) 
sp4 = C_series(2)     # Sp(4)
so8 = D_series(4)     # SO(8)
e6 = E_series(6)      # E₆
g2 = G_series(2)      # G₂

println("Created algebras: SU(3), SO(5), Sp(4), SO(8), E₆, G₂")
println("SU(3) rank: ", lie_rank(su3), ", dual Coxeter: ", dual_coxeter_number(su3))
println()

# Create irreducible representations using Dynkin labels
fund = Irrep(su3, 1, 0)        # Fundamental rep of SU(3)
adj = adjoint_irrep(su3)       # Adjoint representation
high_weight = Irrep(su3, 3, 2) # Higher weight rep

println("Fundamental [1,0]: dimension = ", dimension(fund))
println("Adjoint [1,1]: dimension = ", dimension(adj))
println("High weight [3,2]: dimension = ", dimension(high_weight))
println()

# ============================================================================
# Example 2: Representation Properties
# ============================================================================
println("\nExample 2: Representation Properties")
println("-" ^ 60)

println("Fundamental [1,0]:")
println("  Dimension: ", dimension(fund))
println("  Quadratic Casimir: ", quadratic_casimir(fund))
println("  Dynkin index: ", dynkin_index(fund))
println("  Is self-conjugate? ", is_self_conjugate(fund))
println("  Conjugate: ", conjugate(fund))
println()

anti_fund = conjugate(fund)
println("Checking 3 ⊗ 3̄ contains trivial:")
prod = fund ⊗ anti_fund
trivial = Irrep(su3, 0, 0)
println("  ", prod)
println("  Contains trivial? ", prod[trivial] > 0)
println("  Multiplicity of trivial: ", prod[trivial])
println()

# ============================================================================
# Example 3: Characters - Complete Weight Decomposition
# ============================================================================
println("\nExample 3: Computing Characters")
println("-" ^ 60)

# Character gives all weights with their multiplicities
println("Character of SU(3) fundamental [1,0]:")
char = character(fund)
println("  Total number of weights: ", length(char))
println("  First few weights:")
for i in 1:min(5, length(char))
    wm = char[i]
    println("    Weight ", wm.weight.coordinates, " with multiplicity ", wm.multiplicity)
end
println()

# Can index characters to get multiplicity of specific weights
w = Weight(su3, [1, 0])
println("Multiplicity of weight [1,0] in fundamental: ", char[w])
println()

# Larger example
println("Character of SU(3) adjoint [1,1]:")
adj_char = character(adj)
println("  Total weights: ", length(adj_char), " (dimension = ", dimension(adj), ")")
println("  Number of distinct weights: ", distinct_weights(adj_char))
println()

# ============================================================================
# Example 4: Tensor Products using Klimyk's Formula
# ============================================================================
println("\nExample 4: Tensor Products")
println("-" ^ 60)

println("SU(3) tensor products:")
println("  3 ⊗ 3 = ", fund ⊗ fund)
println("  3 ⊗ 3̄ = ", fund ⊗ anti_fund)
println("  8 ⊗ 8 = ", adj ⊗ adj)
println()

# Verify dimensions
result = adj ⊗ adj
total_dim = sum(mult * dimension(irrep) for (irrep, mult) in irreps(result))
println("Dimension check for 8 ⊗ 8:")
println("  Sum of component dimensions: ", total_dim)
println("  Expected (8 × 8): ", dimension(adj)^2)
println("  Match: ", total_dim == dimension(adj)^2)
println()

# ============================================================================
# Example 5: Reducible Representations
# ============================================================================
println("\nExample 5: Reducible Representations")
println("-" ^ 60)

# Create reducible representations
rep1 = Rep(fund)           # Convert irrep to Rep type
rep2 = Rep(anti_fund)
reducible = rep1 ⊕ rep2    # Direct sum

println("3 ⊕ 3̄ = ", reducible)
println()

# Tensor products with reducible reps
println("(3 ⊕ 3̄) ⊗ 8:")
result = reducible ⊗ Rep(adj)
println("  ", result)
println()

# Can also build reducible reps from decompositions
println("Building 3 ⊗ 3 ⊗ 3:")
triple = fund ⊗ fund ⊗ fund
println("  ", triple)
println("  Number of irrep components: ", length(irreps(triple)))
println()

# Access individual multiplicities
println("Checking multiplicities in 8 ⊗ 8:")
for (irrep, mult) in sort(collect(irreps(adj ⊗ adj)), by=x->dimension(x[1]))
    labels = dynkin_labels(irrep)
    println("  [", join(labels, ","), "]: multiplicity ", mult, 
            ", dimension ", dimension(irrep))
end
println()

# ============================================================================
# Example 6: Weyl Group Operations
# ============================================================================
println("\nExample 6: Weyl Group")
println("-" ^ 60)

weight = Weight(su3, [2, 1])
println("Original weight: ", weight.coordinates)
println("Is dominant? ", is_dominant(weight))
println()

# Weyl reflections
reflected = weyl_reflection(weight, simple_roots(su3)[1])
println("After reflection by α₁: ", reflected.coordinates)
println()

# Reflect to dominant chamber
weight2 = Weight(su3, [-1, 2])
dom, sign = reflect_to_dominant(weight2)
println("Weight ", weight2.coordinates, " → dominant: ", dom.coordinates, 
        " (sign: ", sign, ")")
println()

# Weyl orbit
orbit = weyl_orbit(weight)
println("Weyl orbit of [2,1]:")
println("  Weights with sign +1: ", length(orbit[1]))
println("  Weights with sign -1: ", length(orbit[-1]))
println("  Total: ", length(orbit[1]) + length(orbit[-1]))
println("  Weyl group order: ", weyl_group_order(su3))
println()

# ============================================================================
# Example 7: Different Lie Algebra Series
# ============================================================================
println("\nExample 7: Various Lie Algebras")
println("-" ^ 60)

println("SO(5) = B₂:")
vector_so5 = Irrep(so5, 1, 0)
spinor_so5 = Irrep(so5, 0, 1)
println("  Vector [1,0]: dim = ", dimension(vector_so5))
println("  Spinor [0,1]: dim = ", dimension(spinor_so5))
println("  Vector ⊗ Vector = ", vector_so5 ⊗ vector_so5)
println()

println("Sp(4) = C₂:")
fund_sp4 = Irrep(sp4, 1, 0)
println("  Fundamental [1,0]: dim = ", dimension(fund_sp4))
println("  Self-conjugate? ", is_self_conjugate(fund_sp4))
println()

println("G₂:")
fund_g2 = Irrep(g2, 1, 0)
adj_g2 = adjoint_irrep(g2)
println("  7-dimensional [1,0]: dim = ", dimension(fund_g2))
println("  14-dimensional adjoint: dim = ", dimension(adj_g2))
println("  7 ⊗ 7 = ", fund_g2 ⊗ fund_g2)
println()

println("E₆:")
fund_e6 = Irrep(e6, 1, 0, 0, 0, 0, 0)
adj_e6 = adjoint_irrep(e6)
println("  27-dimensional [1,0,0,0,0,0]: dim = ", dimension(fund_e6))
println("  78-dimensional adjoint: dim = ", dimension(adj_e6))
println()

println("SO(8) = D₄ (Triality):")
vector_8 = Irrep(so8, 1, 0, 0, 0)
spinor_8s = Irrep(so8, 0, 0, 1, 0)
spinor_8c = Irrep(so8, 0, 0, 0, 1)
println("  Vector [1,0,0,0]: dim = ", dimension(vector_8))
println("  Spinor [0,0,1,0]: dim = ", dimension(spinor_8s))
println("  Conjugate spinor [0,0,0,1]: dim = ", dimension(spinor_8c))
println("  These three 8-dimensional reps are related by triality!")
println()

# ============================================================================
# Example 8: Performance - Large Representations
# ============================================================================
println("\nExample 8: Performance on Large Representations")
println("-" ^ 60)

# F₄ example with over 1 million dimensions
f4 = F_series(4)
large_rep = Irrep(f4, 2, 0, 1, 1)
println("F₄ representation [2,0,1,1]:")
println("  Dimension: ", dimension(large_rep))

# Character computation
println("  Computing character...")
@time char_large = character(large_rep)
println("  Total weights: ", length(char_large))
println()

println("=== End of Examples ===")
