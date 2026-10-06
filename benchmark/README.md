# Benchmarks

* `benchmarks.jl`: timings for characters, tensor products and plethysms.
  Run with `julia --project benchmark/benchmarks.jl` (add `--once` to time only
  the first call, compilation included).
* `cosh_table.jl`: decomposes `χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2})` for a
  list of real irreps `r` with `cosh_product`. It writes
  `cosh_product_table.md` (summary) and `cosh_product_full.md` (complete
  decompositions). Run with `julia --project -t auto benchmark/cosh_table.jl`.

Every result in the table is checked against `dim R = 2^{|Φ⁺(r)| + m₀}`. The
test suite also checks the method against independent identities: the adjoint
gives `2^rk · V_ρ`, `R ⊗ R = 2^{m₀} Λ•(r)` (via `tensor_product` and `plethysm`),
and all Levi-subalgebra strategies agree.

## Library timings: before and after

Single thread, seconds (the "before" column is the code at commit `8cc4eb0`,
compilation excluded, cut off at 300 s). Produced with `benchmarks.jl`. Times
below 1 ms are at the resolution of the measurement, so those speed-ups are approximate.

| task | case | before | after | speed-up |
|---|---|---:|---:|---:|
| character | A₂ [10,10] (dim 1331) | 0.011 | 0.0003 | 36× |
| character | A₄ [2,1,1,2] (dim 6125) | 0.29 | 0.0004 | 700× |
| character | D₅ [1,1,0,1,1] (dim 36750) | 2.35 | 0.0011 | 2100× |
| character | C₄ [1,1,1,1] (dim 65536) | 0.82 | 0.0007 | 1200× |
| character | G₂ [6,6] (dim 117649) | 0.038 | 0.0004 | 95× |
| character | F₄ [1,1,0,1] (dim 379848) | 2.18 | 0.0017 | 1300× |
| character | E₆ [1,1,0,0,1,1] (dim 4.2·10⁶) | 95.6 | 0.069 | 1400× |
| character | E₇ [0,0,1,0,0,0,1] (dim 10⁸) | > 300 | 0.47 | > 600× |
| character | E₈ [1,0,0,0,0,0,1,0] (dim 779247) | > 300 | 0.15 | > 2000× |
| tensor product | A₂ [8,5] ⊗ [6,7] | 0.0078 | 0.0001 | 78× |
| tensor product | A₄ [1,1,1,1] ⊗ [2,1,0,1] | 0.16 | 0.0001 | 1600× |
| tensor product | D₅ [0,1,0,1,1] ⊗ [1,0,1,0,0] | 1.20 | 0.0001 | 12000× |
| tensor product | G₂ [3,3] ⊗ [4,2] | 0.041 | 0.0001 | 400× |
| tensor product | F₄ [0,0,1,1] ⊗ [1,0,0,1] | 0.41 | 0.0001 | 4000× |
| tensor product | E₆ [1,1,0,0,0,1] ⊗ [0,1,0,0,1,1] | 43.0 | 0.0042 | 10000× |
| tensor product | E₇ [0,0,0,1,0,0,0] ⊗ [0,0,1,0,0,0,0] | 179 | 0.0022 | 80000× |
| tensor product | E₈ [0,…,0,1,1] ⊗ [1,0,…,0,1,0] | > 300 | 0.042 | > 7000× |
| plethysm | A₂ [1,1], s₍₄,₂₎ | 0.62 | 0.0002 | 3000× |
| plethysm | A₂ [2,1], s₍₃,₂,₁₎ | 2.36 | 0.0003 | 8000× |
| plethysm | A₃ [1,0,1], s₍₆₎ | 14.7 | 0.0001 | 10⁵× |
| plethysm | D₄ [1,0,0,0], s₍₄,₃,₁₎ | 183 | 0.0002 | 10⁶× |
| plethysm | G₂ [1,0], s₍₈₎ | 3.86 | 0.0001 | 4·10⁴× |
| plethysm | F₄ [0,0,0,1], s₍₄₎ | 6.14 | < 0.0001 | > 6·10⁴× |
| plethysm | E₆ [1,0,0,0,0,0], s₍₃,₂₎ | > 300 | 0.0001 | > 10⁶× |
| plethysm | E₇ [0,0,0,0,0,1,0], s₍₄₎ | > 300 | 0.0001 | > 10⁶× |
| plethysm | E₈ [0,…,0,1,0], s₍₃₎ | > 300 | 0.0002 | > 10⁶× |
| plethysm | A₁ [1], s₍₁₀,₈,₆,₂₎ | > 300 | 0.0042 | > 7·10⁴× |

Results were checked against the old implementation on 31 fixed and 117
random cases (characters, tensor products, plethysms): all identical.

## `cosh_product` timings

See `cosh_product_table.md`. With 4 threads every listed case except the
SO(16) spinor takes under 3 s. The SO(16) spinor (a single Weyl orbit of 128
weights in rank 8; the answer is the 135 irreps predicted for E₈/Spin(16)) takes
about 220 s. E₆ 650 and the larger E₇/E₈ cases are out of reach. For every Levi
choice, the estimated Racah–Speiser work is ≳ 10¹¹ steps.
