# Representations

Representations in Kiwi.jl are specified by their Dynkin labels, which correspond to the highest weight in the fundamental weight basis.

## Defining Representations

Let's use SO(6) as an example:

```@example so6
using Kiwi

# SO(6) = D_3
so6 = SO(6)

# Define representations by Dynkin labels
vector = Irrep(so6, [1, 0, 0])      # 6-dimensional vector
spinor_plus = Irrep(so6, [0, 1, 0])  # 4-dimensional spinor
spinor_minus = Irrep(so6, [0, 0, 1]) # 4-dimensional spinor

# Using individual arguments also works
vector_alt = Irrep(so6, 1, 0, 0)
nothing # hide
```

## Properties

### Dimension

```@example so6
[dimension(vector), dimension(spinor_plus), dimension(spinor_minus)]  # [6, 4, 4]
```

### Quadratic Casimir

```@example so6
quadratic_casimir(vector)  # 5/4
```

### Dynkin Index

The Dynkin index is the trace normalization:

```@example so6
[dynkin_index(vector), dynkin_index(spinor_plus)]  # [1, 1/2]
```

### Conjugation

```@example so6
[is_self_conjugate(vector), conjugate(spinor_plus) == spinor_minus]  # [true, true]
```

## Special Representations

```@example so6
# Adjoint representation
adj = adjoint_irrep(so6)
[dimension(adj), is_trivial(Irrep(so6, [0, 0, 0]))]  # [15, true]
```

## See Also

- [Tensor Products](tensor_products.md): Decomposing tensor products
- [Characters](characters.md): Weight decompositions
- [API Reference: Representations](../api/representations.md)

