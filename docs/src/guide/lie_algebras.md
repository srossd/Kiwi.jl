# Lie Algebras

Kiwi.jl provides comprehensive support for all simple Lie algebras.

## Creating Lie Algebras

### Classical Algebras

You can create Lie algebras using either the series notation or standard mathematical notation:

```@example alg
using Kiwi # hide

# A series: SU(n+1)
su3 = A_series(2)    # SU(3)
su5 = A_series(4)    # SU(5)

# Or use the SU helper
su3_alt = SU(3)      # Same as A_series(2)

# B series: SO(2n+1)
so5 = B_series(2)    # SO(5)
so9 = B_series(4)    # SO(9)

# Or use the SO helper
so5_alt = SO(5)      # Same as B_series(2)

# C series: Sp(2n)
sp4 = C_series(2)    # Sp(4)
sp6 = C_series(3)    # Sp(6)

# Or use the Sp helper
sp4_alt = Sp(4)      # Same as C_series(2)

# D series: SO(2n)
so8 = D_series(4)    # SO(8)
so10 = D_series(5)   # SO(10)

# Or use the SO helper
so8_alt = SO(8)      # Same as D_series(4)
nothing # hide
```

The helper functions `SU(n)`, `SO(n)`, and `Sp(n)` provide more intuitive notation:
- `SU(n)` creates su(n), which is A_{n-1} (requires n ≥ 2)
- `SO(n)` creates so(n), which is B_k for n=2k+1 or D_k for n=2k (requires n ≥ 3)
- `Sp(n)` creates sp(n), which is C_{n/2} (requires even n ≥ 2)

### Exceptional Algebras

```@example alg
e6 = E_series(6)
e7 = E_series(7)
e8 = E_series(8)
f4 = F_series(4)
g2 = G_series(2)
nothing # hide
```

## Algebra Properties

### Basic Properties

```@example alg
println("SU(3) properties:")
println("--------------")
println("Dimension: ", dimension(su3))
println("Rank: ", lie_rank(su3))
println("Dual Coxeter number: ", dual_coxeter_number(su3))
println("Cartan matrix:\n", cartan_matrix(su3))
```

### Root System

Access the simple roots:

```@example alg
println("Simple roots of SU(3): ", simple_roots(su3))
```

The Weyl vector ρ (half-sum of positive roots):

```@example alg
println("Weyl vector of SU(3): ", weyl_vector(su3))
```

## See Also

- [API Reference: Lie Algebras](../api/algebras.md): Complete function documentation
