# Characters API

## Computing Characters

```@autodocs
Modules = [Kiwi]
Pages = ["characters.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.character, Kiwi.dimension, Kiwi.highest_weight, Kiwi.dominant_weights]
```

## Dominant characters

```@autodocs
Modules = [Kiwi]
Pages = ["core.jl"]
Order = [:function]
Filter = t -> t == Kiwi.dominant_character
```
