# Timing benchmarks for characters, tensor products and plethysms.
#
#   julia --project benchmark/benchmarks.jl
#
# Each case is run once to compile and once more for the timing.  Pass
# `--once` to time the first call only (useful for very slow implementations).
using Kiwi
using Printf

const ONCE = "--once" in ARGS

const CHARACTER_CASES = [
    (SU(3), [10, 10]), (SU(5), [2, 1, 1, 2]), (SO(10), [1, 1, 0, 1, 1]), (Sp(8), [1, 1, 1, 1]),
    (G_series(2), [6, 6]), (F_series(4), [1, 1, 0, 1]), (E_series(6), [1, 1, 0, 0, 1, 1]),
    (E_series(7), [0, 0, 1, 0, 0, 0, 1]), (E_series(8), [1, 0, 0, 0, 0, 0, 1, 0]),
]

const TENSOR_CASES = [
    ((SU(3), [8, 5]), (SU(3), [6, 7])), ((SU(5), [1, 1, 1, 1]), (SU(5), [2, 1, 0, 1])),
    ((SO(10), [0, 1, 0, 1, 1]), (SO(10), [1, 0, 1, 0, 0])), ((G_series(2), [3, 3]), (G_series(2), [4, 2])),
    ((F_series(4), [0, 0, 1, 1]), (F_series(4), [1, 0, 0, 1])),
    ((E_series(6), [1, 1, 0, 0, 0, 1]), (E_series(6), [0, 1, 0, 0, 1, 1])),
    ((E_series(7), [0, 0, 0, 1, 0, 0, 0]), (E_series(7), [0, 0, 1, 0, 0, 0, 0])),
    ((E_series(8), [0, 0, 0, 0, 0, 0, 1, 1]), (E_series(8), [1, 0, 0, 0, 0, 0, 1, 0])),
]

const PLETHYSM_CASES = [
    ((SU(3), [1, 1]), [4, 2]), ((SU(3), [2, 1]), [3, 2, 1]), ((SU(4), [1, 0, 1]), [6]),
    ((SO(8), [1, 0, 0, 0]), [4, 3, 1]), ((G_series(2), [1, 0]), [8]), ((F_series(4), [0, 0, 0, 1]), [4]),
    ((E_series(6), [1, 0, 0, 0, 0, 0]), [3, 2]), ((E_series(7), [0, 0, 0, 0, 0, 1, 0]), [4]),
    ((E_series(8), [0, 0, 0, 0, 0, 0, 1, 0]), [3]), ((SU(2), [1]), [10, 8, 6, 2]),
]

function timeit(f)
    if ONCE
        t = @elapsed r = f()
    else
        f()
        t = @elapsed r = f()
    end
    return t, r
end

Kiwi_reset() = isdefined(Kiwi, :_DOMCHAR_CACHE) && empty!(Kiwi._DOMCHAR_CACHE)

println("## Characters (full weight system)\n")
println("| algebra | irrep | dim | # weights | time (s) |\n|---|---|---:|---:|---:|")
for (g, l) in CHARACTER_CASES
    rep = Irrep(g, l)
    t, c = timeit(() -> (Kiwi_reset(); character(rep)))
    @printf("| %s | %s | %d | %d | %.4f |\n", g, l, dimension(rep), length(c.weights), t); flush(stdout)
end

println("\n## Tensor products\n")
println("| algebra | irreps | # components | time (s) |\n|---|---|---:|---:|")
for ((g, a), (_, b)) in TENSOR_CASES
    t, r = timeit(() -> (Kiwi_reset(); tensor_product(Irrep(g, a), Irrep(g, b))))
    @printf("| %s | %s ⊗ %s | %d | %.4f |\n", g, a, b, length(r.components), t); flush(stdout)
end

println("\n## Plethysms\n")
println("| algebra | irrep | partition | # components | time (s) |\n|---|---|---|---:|---:|")
for ((g, a), p) in PLETHYSM_CASES
    t, r = timeit(() -> (Kiwi_reset(); plethysm(Irrep(g, a), SymmetricIrrep(p))))
    @printf("| %s | %s | %s | %d | %.4f |\n", g, a, p, length(r.components), t); flush(stdout)
end
