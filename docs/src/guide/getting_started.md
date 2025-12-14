# Getting Started

## Installation

Install Kiwi.jl using Julia's package manager:

```julia
using Pkg
Pkg.add("Kiwi")
```

## Basic Usage

```@example reps
using Kiwi

# Create a Lie algebra
SU3 = SU(3)

# Get properties of representations
fund = Irrep(SU3, 1, 0)
println("Properties of fundamental representation:")
println("\tDimension: ", dimension(fund))
println("\tCasimir:   ", quadratic_casimir(fund))
println("\tIndex:     ", dynkin_index(fund))
```

## Key Concepts

### Lie Algebras

Kiwi.jl supports all classical and exceptional simple Lie algebras:
- **A series**: SU(n+1)
- **B series**: SO(2n+1)
- **C series**: Sp(2n)
- **D series**: SO(2n)
- **E series**: E₆, E₇, E₈
- **F series**: F₄
- **G series**: G₂

### Irreducible Representations

Representations are specified by their highest weight in terms of Dynkin labels:

```@example reps
# SU(3) representations
anti_fund = Irrep(SU(3), 0, 1)   # 3̄
adjoint = Irrep(SU(3), 1, 1)     # 8

# SO(5) representations
vector = Irrep(SO(5), 1, 0)      # 5
spinor = Irrep(SO(5), 0, 1)      # 4

print([rep => dimension(rep) for rep in [fund, anti_fund, adjoint, vector, spinor]])
```

### Characters

The character of a representation is the complete list of weights with their multiplicities:

```@example reps
print(character(fund))
```

### Tensor Products

Compute tensor product decompositions using the ⊗ operator:

```@example reps
print(fund ⊗ anti_fund)
```

## Next Steps

- [Lie Algebras](lie_algebras.md)
- [Representations](representations.md)
- [Tensor Products](tensor_products.md)
- [Characters](characters.md)
- [Weyl groups](weyl_groups.md)
