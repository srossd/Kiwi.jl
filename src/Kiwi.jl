module Kiwi

using LinearAlgebra

# Export main types
export LieAlgebra, Irrep, TensorProduct, WeightMultiplicity, Weight, Rep, Character
export AbstractCharacter, LazyCharacter
export A_series, B_series, C_series, D_series, E_series, F_series, G_series
export SU, SO, Sp

# Export utility functions
export lie_rank, adjoint_irrep, dual_coxeter_number, cartan_matrix, weyl_vector
export dimension, quadratic_casimir, dynkin_index
export conjugate, is_trivial, is_self_conjugate
export highest_weight
export tensor_product
export character
export weyl_reflection, simple_reflection, is_dominant, dynkin_labels, inner_product
export reflect_to_dominant, weyl_orbit, weyl_group_order, simple_roots
export irreps
export ⊕, ⊗

# Include source files
include("lie_algebras.jl")
include("weights.jl")
include("irreps.jl")
include("invariants.jl")
include("characters.jl")
include("weyl.jl")
include("reducible.jl")

end # module
