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

## See Also

- [API Reference: Tensor Products](../api/tensor_products.md)
