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

## Plethysms of reducible representations

`plethysm` also accepts a reducible `Rep`:

```@example su3
V = Rep(fund) + Rep(adj)
plethysm(V, SymmetricIrrep([2]))   # S²(3 ⊕ 8)
```

## Spinor-type products over weights

For a self-conjugate representation `r`, [`cosh_product`](@ref) decomposes

```math
\chi = 2^{m_0} \prod_{w \in \Phi^+(r)} \left(e^{w/2} + e^{-w/2}\right),
```

where ``\Phi^+(r)`` contains one weight from each pair ``\pm w`` of non-zero weights
and ``m_0`` is the multiplicity of the zero weight.  For a real `r` this is (up to a
power of 2) the spinor representation of ``\mathfrak{so}(\dim r)`` restricted to the
algebra; for the adjoint representation it is ``2^{\mathrm{rk}}`` copies of the irrep
with highest weight ``\rho`` (Kostant).

```@example su3
cosh_product(adj)               # 4 × [1,1]
cosh_product(Irrep(su3, [2, 2]))
```

With `half = false` the product runs over all weights and gives the character of the
exterior algebra ``\Lambda^\bullet r``.  Use [`frobenius_schur_indicator`](@ref) /
[`is_real`](@ref) to check whether an irrep is real.
