# Tensor Products

Compute the decomposition of tensor products of representations into irreducible components using the `⊗` operator.

## Basic Examples

Let's use SU(3) to see how tensor products work:

```@example su3
using Kiwi

# SU(3) = A_2
su3 = SU(3)

# Define the fundamental representations
fund = Irrep(su3, [1, 0])       # 3
anti_fund = Irrep(su3, [0, 1])  # 3̄
adj = Irrep(su3, [1, 1])        # 8
nothing # hide
```

### Fundamental × Fundamental

```@example su3
# 3 ⊗ 3 = 3̄ ⊕ 6
result = fund ⊗ fund
println("[1,0] ⊗ [1,0] = ", result)
println("Dimension of result: ", dimension(result))
```

### Adjoint × Adjoint

```@example su3
# 8 ⊗ 8 = 1 ⊕ 8 ⊕ 8 ⊕ 10 ⊕ 10̄ ⊕ 27
result = adj ⊗ adj
println("[1,1] ⊗ [1,1] = ", result)
println("Dimension of result: ", dimension(result))
```

## Accessing Components

```@example su3
# Get multiplicity of a specific irrep in 3 ⊗ 3̄
result = fund ⊗ anti_fund
trivial = Irrep(su3, [0, 0])
result[trivial]  # 1
```

## Plethysms

Plethysms compute symmetric and antisymmetric powers of representations. Use the `symmetric_power` and `antisymmetric_power` helper functions, or the general `plethysm` function with symmetric group irreps.

### Symmetric Powers

```@example su3
# S²([1,0]) - symmetric square of the fundamental
result = symmetric_power(2, fund)
println("S²([1,0]) = ", result)
println("Dimension: ", dimension(result))  # 6
```

### Antisymmetric Powers

```@example su3
# Λ²([1,0]) - antisymmetric square of the fundamental
result = antisymmetric_power(2, fund)
println("Λ²([1,0]) = ", result)
println("Dimension: ", dimension(result))  # 3
```

### Higher Powers

```@example su3
# S³([1,1]) - symmetric cube of the adjoint
result = symmetric_power(3, adj)
println("S³([1,1]) = ", result)
println("Total dimension: ", dimension(result))  # 120
```

### General Plethysms

For more general plethysms, use the `plethysm` function with symmetric group representations:

```@example su3
# Mixed symmetry: [2,1] representation of S₃
rho = SymmetricIrrep([2, 1])
result = plethysm(fund, rho)
println("[2,1]([1,0]) = ", result)
```

## See Also

- [API Reference: Tensor Products](../api/tensor_products.md)
