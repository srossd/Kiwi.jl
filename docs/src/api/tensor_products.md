# Tensor Products API

## Computing Tensor Products

```@autodocs
Modules = [Kiwi]
Pages = ["reducible.jl"]
Order = [:function]
Filter = t -> t == Kiwi.tensor_product
```

## Plethysms

```@autodocs
Modules = [Kiwi]
Pages = ["plethysms.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.plethysm, Kiwi.symmetric_power, Kiwi.antisymmetric_power]
```

## Unicode Operator

```julia
rep1 ⊗ rep2  # Same as tensor_product(rep1, rep2)
```

Type `\otimes<TAB>` in Julia REPL to get ⊗.
