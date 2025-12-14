# Kiwi.jl

*A Julia package for Lie algebra representation theory*

Kiwi.jl provides efficient tools for working with Lie algebras and their representations.

## Quick Example

```@example
using Kiwi

# Create SU(3) and fundamental representations
fund = Irrep(SU(3), 1, 0)
anti_fund = conjugate(fund)

# Compute tensor products
println("3 ⊗ 3= ", fund ⊗ anti_fund)

# Check if trivial representation appears
trivial = Irrep(SU(3), 0, 0)
result = fund ⊗ anti_fund
println("Multiplicity of trivial: ", result[trivial])

# Compute characters
char = character(fund)
println("Dimension: ", sum(wm.multiplicity for wm in char))
```

## Installation

```julia
using Pkg
Pkg.add("Kiwi")
```

## Contents

```@contents
Pages = [
    "guide/getting_started.md",
    "guide/lie_algebras.md",
    "guide/representations.md",
    "guide/characters.md",
    "guide/tensor_products.md",
    "guide/weyl_groups.md",
]
Depth = 2
```

## Index

```@index
```
