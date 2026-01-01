# Representations API

## Creating Representations

```@autodocs
Modules = [Kiwi]
Pages = ["irreps.jl"]
Order = [:function]
Filter = t -> t == Kiwi.adjoint_irrep
```

## Finding Representations

```@autodocs
Modules = [Kiwi]
Pages = ["irreps.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.irreps_up_to_dim, Kiwi.clear_irreps_cache!]
```

## Representation Properties

```@autodocs
Modules = [Kiwi]
Pages = ["irreps.jl", "invariants.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.dimension, Kiwi.quadratic_casimir, Kiwi.dynkin_index, Kiwi.conjugate, Kiwi.is_self_conjugate, Kiwi.is_trivial, Kiwi.dynkin_labels, Kiwi.congruency_class, Kiwi.dimension_label]
```