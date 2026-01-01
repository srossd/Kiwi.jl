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

You can find all irreps up to a given dimension:

```@example so6
irreps_up_to_dim(so6, 10)  # All SO(6) irreps with dimension ≤ 10
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

### Congruency Class

The congruency class specifies how a representation acts on the center of the group:

```@example so6
# SO(8) triality: three distinct 8-dimensional irreps
so8 = D_series(4)
[congruency_class(Irrep(so8, 1, 0, 0, 0)).value,  # vector: (0, 1)
 congruency_class(Irrep(so8, 0, 0, 1, 0)).value,  # spinor: (1, 1)
 congruency_class(Irrep(so8, 0, 0, 0, 1)).value]  # conjugate spinor: (1, 0)
```

### Conjugation

```@example so6
[is_self_conjugate(vector), conjugate(spinor_plus) == spinor_minus]  # [true, true]
```

### Dimension Label

The dimension label provides a compact notation that distinguishes representations with the same dimension:

```@example so6
# SU(4) has three dimension-20 representations
su4 = A_series(3)
[dimension_label(Irrep(su4, 0, 1, 1)),  # "20"
 dimension_label(Irrep(su4, 1, 1, 0)),  # "20b"
 dimension_label(Irrep(su4, 0, 2, 0))]  # "20′"
```

## Special Representations

```@example so6
# Adjoint representation
adj = adjoint_irrep(so6)
[dimension(adj), is_trivial(Irrep(so6, [0, 0, 0]))]  # [15, true]
```

## Representation Tables

You can generate comprehensive tables of representations:

```@example tables
using Kiwi

# SU(3) table
println("SU(3) Representations (dimension ≤ 15)\n")
su3 = A_series(2)
irreps = sort(irreps_up_to_dim(su3, 15), by=r->dimension(r))

println("┌──────────────┬─────────────────┬──────────────┬──────────────┬──────────┐")
println("│ Dim Label    │ Dynkin Labels   │ C₂           │ Dynkin Index │ Cong. Cl │")
println("├──────────────┼─────────────────┼──────────────┼──────────────┼──────────┤")
for rep in irreps
    labels = dynkin_labels(rep)
    dim_label = dimension_label(rep)
    cas = round(quadratic_casimir(rep), digits=3)
    dyn_idx = dynkin_index(rep)
    cc = congruency_class(rep).value
    
    println("│ ", rpad(dim_label, 12), " │ ",
            rpad("[$(join(labels, ","))]", 15), " │ ",
            rpad(string(cas), 12), " │ ",
            rpad(string(dyn_idx), 12), " │ ",
            rpad(string(cc), 8), " │")
end
println("└──────────────┴─────────────────┴──────────────┴──────────────┴──────────┘")
nothing # hide
```

```@example tables
# SO(8) table showing triality
println("\nSO(8) Representations (dimension ≤ 60)\n")
so8 = D_series(4)
irreps = sort(irreps_up_to_dim(so8, 60), by=r->dimension(r))

println("┌──────────────┬─────────────────┬──────────────┬──────────────┬──────────┐")
println("│ Dim Label    │ Dynkin Labels   │ C₂           │ Dynkin Index │ Cong. Cl │")
println("├──────────────┼─────────────────┼──────────────┼──────────────┼──────────┤")
for rep in irreps
    labels = dynkin_labels(rep)
    dim_label = dimension_label(rep)
    cas = round(quadratic_casimir(rep), digits=3)
    dyn_idx = dynkin_index(rep)
    cc = congruency_class(rep).value
    cc_str = cc === nothing ? "trivial" : string(cc)
    
    println("│ ", rpad(dim_label, 12), " │ ",
            rpad("[$(join(labels, ","))]", 15), " │ ",
            rpad(string(cas), 12), " │ ",
            rpad(string(dyn_idx), 12), " │ ",
            rpad(cc_str, 8), " │")
end
println("└──────────────┴─────────────────┴──────────────┴──────────────┴──────────┘")
nothing # hide
```

## See Also

- [Tensor Products](tensor_products.md): Decomposing tensor products
- [Characters](characters.md): Weight decompositions
- [API Reference: Representations](../api/representations.md)

