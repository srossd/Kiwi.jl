# Kiwi.jl

A Julia package for computations with irreducible representations of Lie algebras.

## Features

- **Lie Algebra Types**: Support for classical Lie algebras (A, B, C, D series)
- **Irreducible Representations**: Work with irreps specified by Dynkin labels
- **Representation Invariants**: 
  - Compute dimensions using Weyl's dimension formula
  - Calculate quadratic Casimir eigenvalues
  - Compute Dynkin indices
- **Tensor Products**: Decompose tensor products of irreps into irreducible components

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

# Create a Lie algebra (SU(3))
su3 = A_series(2)

# Create irreducible representations by their Dynkin labels
trivial = Irrep(su3, [0, 0])
fundamental = Irrep(su3, [1, 0])
anti_fundamental = Irrep(su3, [0, 1])
adjoint = Irrep(su3, [1, 1])

# Compute representation properties
dimension(fundamental)        # 3
dimension(adjoint)            # 8
quadratic_casimir(adjoint)    # Casimir eigenvalue
dynkin_index(adjoint)         # Dynkin index

# Check conjugation
conjugate(fundamental) == anti_fundamental  # true
is_self_conjugate(adjoint)                   # true

# Tensor product decomposition
product = fundamental ⊗ anti_fundamental
# or equivalently:
product = tensor_product(fundamental, anti_fundamental)
```

## Supported Lie Algebras

- **A_n**: `A_series(n)` - Special unitary algebras su(n+1)
- **B_n**: `B_series(n)` - Odd orthogonal algebras so(2n+1)
- **C_n**: `C_series(n)` - Symplectic algebras sp(2n)
- **D_n**: `D_series(n)` - Even orthogonal algebras so(2n)

## Examples

### SU(3) Representations

```julia
su3 = A_series(2)

# The 3 representation
fund = Irrep(su3, [1, 0])
dimension(fund)  # 3

# The 6 representation (symmetric square)
sym2 = Irrep(su3, [2, 0])
dimension(sym2)  # 6

# The 15 representation
rep_15 = Irrep(su3, [2, 1])
dimension(rep_15)  # 15
```

### SO(5) Representations

```julia
so5 = B_series(2)

# Vector representation
vector = Irrep(so5, [1, 0])
dimension(vector)  # 5

# Spinor representation
spinor = Irrep(so5, [0, 1])
dimension(spinor)  # 4
```

### Tensor Products

```julia
su3 = A_series(2)
fund = Irrep(su3, [1, 0])
anti_fund = Irrep(su3, [0, 1])

# Decompose 3 ⊗ 3̄
product = fund ⊗ anti_fund

# Iterate over components
for (irrep, multiplicity) in product
    println("$(irrep): multiplicity $(multiplicity)")
end
```

## Roadmap

Future features planned:
- Full Littlewood-Richardson rule implementation for tensor products
- 6j-symbols and recoupling coefficients
- Support for exceptional Lie algebras (E, F, G series)
- Branching rules for subgroup embeddings
- Character formulas

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("Kiwi")
```

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
