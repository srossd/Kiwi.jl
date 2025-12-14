# Weyl Groups API

## Weyl Reflections

```@autodocs
Modules = [Kiwi]
Pages = ["weyl.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.weyl_reflection, Kiwi.simple_reflection]
```

## Dominant Weights

```@autodocs
Modules = [Kiwi]
Pages = ["weyl.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.is_dominant, Kiwi.reflect_to_dominant]
```

## Weyl Orbits

```@autodocs
Modules = [Kiwi]
Pages = ["weyl.jl"]
Order = [:function]
Filter = t -> t in [Kiwi.weyl_orbit, Kiwi.weyl_group_order]
```
