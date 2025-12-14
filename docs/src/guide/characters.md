# Characters

Characters provide the complete weight decomposition of a representation, giving all weights with their multiplicities.

## Computing Characters

Let's start with some SU(3) examples:

```@example su3chars
using Kiwi

su3 = SU(3)
fund = Irrep(su3, [1, 0])      # 3
anti_fund = Irrep(su3, [0, 1])  # 3̄
adj = Irrep(su3, [1, 1])       # 8

# Compute characters
char_fund = character(fund)
char_anti = character(anti_fund)
char_adj = character(adj)
```

Check dimensions and number of distinct weights:

```@example su3chars
[dimension(char_fund), length(char_fund)]  # [3, 3]
```

```@example su3chars
[dimension(char_anti), length(char_anti)]  # [3, 3]
```

```@example su3chars
[dimension(char_adj), length(char_adj)]    # [8, 7]
```

### Accessing Weight Multiplicities

```@example su3chars
# Query a specific weight in the adjoint
w = Weight(su3, [0, 0])
char_adj[w]  # 2 (the zero weight has multiplicity 2 in the adjoint)
```

## Lazy Characters

For large representations, `LazyCharacter` computes multiplicities on demand without computing the full character:

```@example lazy
using Kiwi

# Large E₆ representation: [1,0,1,1,0,1] has dimension 138,881,925
e6 = E_series(6)
rep = Irrep(e6, [2, 1, 0, 0, 0, 1])

println("Dimension: ", dimension(rep))

# Create lazy character (instant)
lazy_char = LazyCharacter(rep)
```

Query specific weights quickly:

```@repl lazy
# Highest weight
@time lazy_char[highest_weight(rep)]
length(lazy_char.computed)
@time lazy_char[Weight(e6, [-1, -2, 1, 1, -2, 1])]
length(lazy_char.computed)
```

For comparison, computing the full character takes much longer:

```@repl lazy
@time char_full = character(rep);
```

The lazy character only computes what you need, making it efficient for querying specific multiplicities in large representations.

## See Also

- [API Reference: Characters](../api/characters.md)
