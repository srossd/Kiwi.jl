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

See [`benchmark/`](benchmark/README.md) for timings: compared with the previous
rational-arithmetic implementation, characters are 10²–10³× faster, tensor
products 10²–10⁵× and plethysms 10³–10⁶× (e.g. the E₇ character of dimension 10⁸
in 0.5 s, `56 ⊗ 56`-type E₇ products in milliseconds).

```julia
# 56 ⊗ 56 of E7, characters of large E8 irreps, plethysms of E-type reps: all fast
tensor_product(Irrep(E_series(7), [0,0,0,0,0,1,0]), Irrep(E_series(7), [0,0,0,0,0,1,0]))
character(Irrep(E_series(8), [1,0,0,0,0,0,1,0]))
plethysm(Irrep(E_series(6), [1,0,0,0,0,0]), SymmetricIrrep([3, 2]))
```

### Comparison with LieART

Wall-clock seconds for the same computations in Kiwi and in
[LieART](https://lieart.hepforge.org/) 2.1.0 (Mathematica). Kiwi times are the
best of 5 calls on one thread, compilation excluded; LieART times are one call
after a warm-up. **The two columns come from different machines**: Kiwi ran on a
4-core cloud VM and LieART on a different computer, so the ratios are indicative only.
The two programs agree on every irrep dimension and on the number of irreducible
components of every tensor product below. For characters, LieART returns the full
list of weights (one entry per dimension) while Kiwi returns the distinct weights
with multiplicities. LieART has no plethysm function, and the E₇ [0,0,1,0,0,0,1]
character (dim 10⁸), which exceeded LieART's 300 s cut-off, is omitted.
Scripts: `benchmark/benchmarks.jl` and `benchmark/lieart_benchmark.wl`.

| task | case | LieART (s) | Kiwi (s) | ratio |
|---|---|---:|---:|---:|
| character | A₂ [10,10] (dim 1331) | 0.12 | 0.00025 | 480× |
| character | A₄ [2,1,1,2] (dim 6125) | 0.26 | 0.0025 | 100× |
| character | D₅ [1,1,0,1,1] (dim 36750) | 1.49 | 0.0010 | 1500× |
| character | C₄ [1,1,1,1] (dim 65536) | 2.08 | 0.00067 | 3100× |
| character | G₂ [6,6] (dim 117649) | 3.50 | 0.00037 | 9400× |
| character | F₄ [1,1,0,1] (dim 379848) | 19.4 | 0.0015 | 13000× |
| character | E₆ [1,1,0,0,1,1] (dim 4.2·10⁶) | 284 | 0.018 | 16000× |
| character | E₈ [1,0,0,0,0,0,1,0] (dim 779247) | 198 | 0.031 | 6400× |
| tensor product | A₂ [8,5] ⊗ [6,7] | 10.1 | 0.000070 | 1.4·10⁵× |
| tensor product | A₄ [1,1,1,1] ⊗ [2,1,0,1] | 0.028 | 0.000056 | 500× |
| tensor product | D₅ [0,1,0,1,1] ⊗ [1,0,1,0,0] | 0.095 | 0.000097 | 980× |
| tensor product | G₂ [3,3] ⊗ [4,2] | 0.45 | 0.000087 | 5100× |
| tensor product | F₄ [0,0,1,1] ⊗ [1,0,0,1] | 0.18 | 0.000081 | 2200× |
| tensor product | E₆ [1,1,0,0,0,1] ⊗ [0,1,0,0,1,1] | 46.6 | 0.0039 | 12000× |
| tensor product | E₇ [0,0,0,1,0,0,0] ⊗ [0,0,1,0,0,0,0] | 16.6 | 0.0022 | 7500× |
| tensor product | E₈ [0,…,0,1,1] ⊗ [1,0,…,0,1,0] | 107 | 0.022 | 4800× |

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

A table of results for real irreps is in
[`benchmark/cosh_product_table.md`](benchmark/cosh_product_table.md).

