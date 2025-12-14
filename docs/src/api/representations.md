# Representations API

## Creating Representations

```@autodocs
Modules = [Kiwi]
Pages = ["irreps.jl"]
Order = [:function]
Filter = t -> t == Kiwi.adjoint_irrep
```

## Representation Properties

```@autodocs
Modules = [Kiwi]
Pages = ["irreps.jl", "invariants.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.dimension, Kiwi.quadratic_casimir, Kiwi.dynkin_index, Kiwi.conjugate, Kiwi.is_self_conjugate, Kiwi.is_trivial, Kiwi.dynkin_labels]
```