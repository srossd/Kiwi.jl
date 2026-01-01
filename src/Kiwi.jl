module Kiwi

using LinearAlgebra

# Export main types
export LieAlgebra, Irrep, TensorProduct, WeightMultiplicity, Weight, Rep, Character
export CongruencyClass, WeylWord
export A_series, B_series, C_series, D_series, E_series, F_series, G_series
export SU, SO, Sp

# Export symmetric group types and functions
export Partition, SymmetricIrrep, ConjugacyClass, SymmetricCharacter
export partition_size, conjugacy_class_size, all_partitions, character_table
export hook_length

# Export utility functions
export lie_rank, adjoint_irrep, dual_coxeter_number, cartan_matrix, weyl_vector
export dimension, quadratic_casimir, dynkin_index, dimension_label
export conjugate, is_trivial, is_self_conjugate
export congruency_class
export highest_weight
export tensor_product, plethysm, symmetric_power, antisymmetric_power
export character, dominant_weights
export weyl_reflection, simple_reflection, is_dominant, dynkin_labels, inner_product
export reflect_to_dominant, weyl_orbit, weyl_group_order, simple_roots
export longest_weyl_word
export irreps, irreps_up_to_dim, clear_irreps_cache!
export ⊕, ⊗

# Include source files
include("lie_algebras.jl")
include("weights.jl")
include("irreps.jl")
include("invariants.jl")
include("characters.jl")
include("weyl.jl")
include("reducible.jl")
include("symmetric.jl")
include("plethysms.jl")

end # module
