# Lie Algebras API

## Creating Algebras

```@autodocs
Modules = [Kiwi]
Pages = ["lie_algebras.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.A_series, Kiwi.B_series, Kiwi.C_series, Kiwi.D_series, Kiwi.E_series, Kiwi.F_series, Kiwi.G_series]
```

## Algebra Properties

```@autodocs
Modules = [Kiwi]
Pages = ["lie_algebras.jl", "cartan.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.lie_rank, Kiwi.dual_coxeter_number, Kiwi.cartan_matrix]
```

## Root System

```@autodocs
Modules = [Kiwi]
Pages = ["weights.jl"]
Order = [:function]
Filter = t -> t == Kiwi.simple_roots
```
