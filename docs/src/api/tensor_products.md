# Tensor Products API

## Computing Tensor Products

```@autodocs
Modules = [Kiwi]
Pages = ["reducible.jl"]
Order = [:function]
Filter = t -> t == Kiwi.tensor_product
```

## Unicode Operator

```julia
rep1 ⊗ rep2  # Same as tensor_product(rep1, rep2)
```

Type `\otimes<TAB>` in Julia REPL to get ⊗.
