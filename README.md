# Kiwi.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://srossd.github.io/Kiwi.jl/)

A Julia package for computations with irreducible representations of Lie algebras.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/srossd/Kiwi.jl")
```

Or in development mode:

```julia
using Pkg
Pkg.develop(path="/path/to/Kiwi.jl")
```

## Quick Start

```julia
using Kiwi

su3 = SU(3)

# Create irreducible representations by their Dynkin labels
trivial = Irrep(su3, [0, 0])
fundamental = Irrep(su3, [1, 0])
anti_fundamental = Irrep(su3, [0, 1])
adjoint = Irrep(su3, [1, 1])

# Compute representation properties
dimension(fundamental)        # 3
dimension(adjoint)            # 8
quadratic_casimir(adjoint)    # 3
dynkin_index(adjoint)         # 3

# Check conjugation
conjugate(fundamental) == anti_fundamental  # true
is_self_conjugate(adjoint)                   # true

# Tensor product decomposition
product = fundamental ⊗ anti_fundamental
# or equivalently:
product = tensor_product(fundamental, anti_fundamental)

# Compute characters (weight decomposition)
char = character(adjoint)
dimension(char)              # 8 (total dimension)
length(char)                 # 7 (number of distinct weights)

# For large representations, use lazy characters
e6 = E_series(6)
large_rep = Irrep(e6, [2, 1, 0, 0, 0, 1])
lazy = character(large_rep; lazy=true)
lazy[highest_weight(large_rep)]  # Query specific weight multiplicities
```
