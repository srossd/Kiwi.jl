# Kiwi.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://srossd.github.io/Kiwi.jl/)

A Julia package for computations with irreducible representations of Lie algebras.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/srossd/Kiwi.jl")
```

Or in development mode:

```julia
using Pkg
Pkg.develop(path="/path/to/Kiwi.jl")
```

## Quick Start

```julia
using Kiwi

su3 = SU(3)

# Create irreducible representations by their Dynkin labels
trivial = Irrep(su3, [0, 0])
fundamental = Irrep(su3, [1, 0])
anti_fundamental = Irrep(su3, [0, 1])
adjoint = Irrep(su3, [1, 1])

# Compute representation properties
dimension(fundamental)        # 3
dimension(adjoint)            # 8
quadratic_casimir(adjoint)    # 3
dynkin_index(adjoint)         # 3

# Check conjugation
conjugate(fundamental) == anti_fundamental  # true
is_self_conjugate(adjoint)                   # true

# Tensor product decomposition
product = fundamental ⊗ anti_fundamental
# or equivalently:
product = tensor_product(fundamental, anti_fundamental)

# Compute characters (weight decomposition)
char = character(adjoint)
dimension(char)              # 8 (total dimension)
length(char)                 # 7 (number of distinct weights)

# For large representations, use lazy characters
e6 = E_series(6)
large_rep = Irrep(e6, [2, 1, 0, 0, 0, 1])
lazy = character(large_rep; lazy=true)
lazy[highest_weight(large_rep)]  # Query specific weight multiplicities
```

## Performance

All computations run on integer weight vectors with per-algebra data (Cartan
matrix, scaled metric, positive roots, coroots) computed once and cached:

* **Characters** use Freudenthal's formula on *dominant* weights only, with the
  Moody–Patera grouping of roots into orbits of the weight's stabiliser; the full
  weight system is obtained by expanding Weyl orbits.
* **Tensor products** use the Racah–Speiser/Klimyk algorithm with the character of
  the smaller factor.
* **Plethysms** use the recursion
  `n·s_λ[V] = Σ_k ψᵏ(V) ⊗ Σ_{λ/μ k-border strip} (-1)^{ht} s_μ[V]`, which only ever
  multiplies genuine representations by Adams operations (no sums over conjugacy
  classes); any (reducible) `Rep` can be plethysm'd.

Compared with the previous rational-arithmetic implementation, characters are 10²–10³× faster, tensor
products 10²–10⁵× and plethysms 10³–10⁶× (e.g. the E₇ character of dimension 10⁸
in 0.5 s, `56 ⊗ 56`-type E₇ products in milliseconds).

```julia
# 56 ⊗ 56 of E7, characters of large E8 irreps, plethysms of E-type reps: all fast
tensor_product(Irrep(E_series(7), [0,0,0,0,0,1,0]), Irrep(E_series(7), [0,0,0,0,0,1,0]))
character(Irrep(E_series(8), [1,0,0,0,0,0,1,0]))
plethysm(Irrep(E_series(6), [1,0,0,0,0,0]), SymmetricIrrep([3, 2]))
```

As a representative comparison with [LieART](https://lieart.hepforge.org/) 2.1.0
(Mathematica), the full character of the E₆ irrep [1,1,0,0,1,1] (dimension 4.2·10⁶)
takes 284 s in LieART and 0.018 s in Kiwi (≈16000×), and the E₆ tensor product
[1,1,0,0,0,1] ⊗ [0,1,0,0,1,1] takes 46.6 s in LieART and 0.0039 s in Kiwi (≈12000×)
(timings from two different machines, so the ratios are indicative).

## Spinor-type products over weights

`cosh_product(r)` decomposes `χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2})` for a
self-conjugate representation `r` (one weight from each `±w` pair; `m₀` = number of
zero weights).  For real `r` this is the spinor module of `so(dim r)` restricted to
the algebra; for the adjoint it is `2^rank` copies of `[1,1,…,1]`.

```julia
cosh_product(adjoint_irrep(F_series(4)))      # 16 × [1,1,1,1]
cosh_product(Irrep(SO(16), [0,0,0,0,0,0,0,1]))  # 135 irreps (E8/Spin(16))
frobenius_schur_indicator(Irrep(SO(7), [0,0,1]))  # 1 (real)
```


## Recoupling coefficients

Explicit representation matrices, Clebsch–Gordan coefficients and recoupling
symbols for every compact simple Lie algebra, including outer multiplicities.
The conventions are stated in full in the documentation (*Recoupling Coefficients*);
layouts follow [TensorKitSectors.jl](https://github.com/QuantumKitHub/TensorKitSectors.jl).

```julia
a = Irrep(SU(3), 1, 1)
E, F, H = generators(a)              # Chevalley generators (sparse, orthonormal basis)
C = fusion_tensor(a, a, a)           # ⟨8 i; 8 j | 8 k; μ⟩, μ = 1, 2
F = fsymbol(a, a, a, a, a, a)        # F^{888}_8[8, 8], a 2×2×2×2 array
R = rsymbol(a, a, a)                 # R^{88}_8
twoj(a), onej(a), threej(a, a, a)    # 2j, 1j (Frobenius–Schur) and 3j symbols
j(x) = Irrep(SU(2), x)               # spin x/2
fsymbol(j(1), j(1), j(2), j(2), j(0), j(1); T = Rational{BigInt})   # exact: -√(1/3)
```

Clebsch–Gordan coefficients are built from highest-weight vectors (via the
Parthasarathy–Ranga Rao–Varadarajan description of the multiplicity space) and
lowered on demand; F- and R-symbols are evaluated on highest-weight vectors
without forming full Clebsch–Gordan tensors.  For su(N) a Gelfand–Tsetlin basis is
available as well (`basis = :gelfand_tsetlin`).  With TensorKitSectors loaded,
`tensorkit_sector(g)` makes the irreps of `g` a TensorKit sector; they pass
TensorKitSectors' sector test suite.

Compared with [SUNRepresentations.jl](https://github.com/QuantumKitHub/SUNRepresentations.jl)
0.4 (no caches, same machine): SU(5) `1024 ⊗ 24 → 1024` Clebsch–Gordan
coefficients take 0.09 s in Kiwi and 1.9 s in SUNRepresentations, SU(6)
`3675 ⊗ 35 → 3675` 0.75 s and 33 s, and the F-symbol
F^{1024, 24, 24}_{1024}[1024, 24] 0.006 s and 4.4 s.
