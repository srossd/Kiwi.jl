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
