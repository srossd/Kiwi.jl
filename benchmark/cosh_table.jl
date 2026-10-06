# Benchmark: decompose χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2}) for real irreps r.
#
#   julia --project -t auto benchmark/cosh_table.jl [output_dir [cache_file]]
#
# If `cache_file` is given, computed results are stored there and reused on the
# next run (handy when only the formatting changes).
#
# Writes `cosh_product_table.md` (summary table) and `cosh_product_full.md`
# (complete decompositions for results with at most MAX_FULL components).
using Kiwi
using Printf
using Serialization

const CASES = [
    # (algebra, Dynkin labels)
    (SU(2), [4]), (SU(2), [6]), (SU(2), [8]), (SU(2), [10]), (SU(2), [12]), (SU(2), [16]), (SU(2), [20]),
    (SU(3), [2, 2]), (SU(3), [3, 3]), (SU(3), [4, 4]), (SU(3), [5, 5]),
    (SU(4), [0, 1, 0]), (SU(4), [0, 2, 0]), (SU(4), [2, 0, 2]), (SU(4), [0, 3, 0]), (SU(4), [1, 2, 1]),
    (SU(5), [0, 1, 1, 0]), (SU(5), [2, 0, 0, 2]),
    (SU(6), [0, 0, 2, 0, 0]), (SU(6), [0, 1, 0, 1, 0]),
    (SO(5), [1, 0]), (SO(5), [2, 0]), (SO(5), [0, 4]), (SO(5), [1, 2]),
    (SO(7), [1, 0, 0]), (SO(7), [0, 0, 1]), (SO(7), [2, 0, 0]), (SO(7), [0, 0, 2]), (SO(7), [1, 0, 1]),
    (SO(9), [1, 0, 0, 0]), (SO(9), [0, 0, 0, 1]), (SO(9), [2, 0, 0, 0]), (SO(9), [0, 0, 1, 0]), (SO(9), [1, 0, 0, 1]),
    (SO(11), [1, 0, 0, 0, 0]), (SO(11), [2, 0, 0, 0, 0]),
    (Sp(6), [0, 1, 0]), (Sp(6), [0, 2, 0]), (Sp(6), [1, 0, 1]), (Sp(6), [0, 0, 2]),
    (Sp(8), [0, 1, 0, 0]), (Sp(8), [0, 0, 0, 1]),
    (SO(8), [1, 0, 0, 0]), (SO(8), [2, 0, 0, 0]), (SO(8), [0, 0, 1, 1]), (SO(8), [1, 1, 0, 0]),
    (SO(10), [1, 0, 0, 0, 0]), (SO(10), [2, 0, 0, 0, 0]), (SO(10), [0, 0, 1, 0, 0]), (SO(10), [0, 0, 0, 1, 1]),
    (SO(12), [1, 0, 0, 0, 0, 0]), (SO(12), [2, 0, 0, 0, 0, 0]),
    (SO(16), [0, 0, 0, 0, 0, 0, 0, 1]),
    (G_series(2), [1, 0]), (G_series(2), [2, 0]), (G_series(2), [1, 1]), (G_series(2), [3, 0]),
    (G_series(2), [0, 2]), (G_series(2), [4, 0]), (G_series(2), [2, 1]),
    (F_series(4), [0, 0, 0, 1]), (F_series(4), [0, 0, 1, 0]), (F_series(4), [0, 0, 0, 2]),
]

# Hand-picked Levi subalgebras (nodes kept) where `levi = :auto` is not the best choice
const LEVI = Dict((SO(16), [0, 0, 0, 0, 0, 0, 0, 1]) => [1, 2, 3, 4, 5, 7, 8])

const MAX_FULL = 400        # list complete decompositions up to this many components
const MAX_INLINE = 4        # show the decomposition in the summary table up to this many

fmt_labels(l) = "[" * join(l, ",") * "]"

function fmt_mult(m::Integer)
    m < 10^6 && return string(m)
    k = trailing_zeros(m)
    o = m >> k
    return o == 1 ? "2^$k" : "2^$k·$o"
end

# (irrep, multiplicity, dimension), sorted by decreasing dimension
function sorted_components(R)
    items = [(ir, m, dimension_big(ir)) for (ir, m) in R.components]
    sort!(items, by = t -> (-t[3], t[1].dynkin_labels))
end

fmt_term((ir, m, d)) = (m == 1 ? "" : fmt_mult(m) * "·") * fmt_labels(ir.dynkin_labels) * "<sub>$d</sub>"

fmt_rep(R) = join(fmt_term.(sorted_components(R)), " + ")

dimension_big(ir::Irrep) = Kiwi._dimension(Kiwi.algebra_data(ir.algebra), ir.dynkin_labels)

function run_case(g, l)
    r = Irrep(g, l)
    @assert frobenius_schur_indicator(r) == 1 "$g $l is not real"
    doms, mults = dominant_character(g, l)
    m0 = sum((m for (μ, m) in zip(doms, mults) if all(iszero, μ)); init = 0)
    npos = (dimension(r) - m0) ÷ 2
    t = @elapsed R = cosh_product(r; levi = get(LEVI, (g, l), :auto))
    E = npos + m0
    total = sum(dimension_big(ir) * big(m) for (ir, m) in R.components)
    @assert total == big(2)^E "dimension check failed for $g $l"
    return (g = g, l = l, dim = dimension(r), m0 = m0, npos = npos, E = E, R = R, time = t)
end

function main(outdir, cachefile = nothing)
    cache = (cachefile !== nothing && isfile(cachefile)) ? deserialize(cachefile) : Dict()
    results = []
    for (g, l) in CASES
        key = (string(g), l)
        res = haskey(cache, key) ? cache[key] : run_case(g, l)
        cache[key] = res
        cachefile !== nothing && serialize(cachefile, cache)
        @printf("%-5s %-22s dim=%-5d components=%-7d %.3fs\n", string(g), fmt_labels(l), res.dim,
                length(res.R.components), res.time)
        flush(stdout)
        push!(results, res)
    end

    open(joinpath(outdir, "cosh_product_table.md"), "w") do io
        println(io, "# `cosh_product` for real irreps\n")
        println(io, "R is the representation with character χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2}), ",
                    "where Φ⁺(r) holds one weight from each ±pair of non-zero weights of r ",
                    "(so #Φ⁺(r) = (dim r − m₀)/2) and m₀ is the multiplicity of the zero weight. ",
                    "Components are written `multiplicity·[Dynkin labels]` with the dimension as a subscript; ",
                    "for long decompositions the two largest components are shown and the complete list (up to ",
                    "$MAX_FULL components) is in `cosh_product_full.md`. Times are wall-clock seconds with 4 threads ",
                    "(first call, i.e. including compilation for the first row).\n")
        println(io, "| G | r | dim r | #Φ⁺(r) | m₀ | dim R | # irreps in R | R (or its largest components) | time (s) |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---|---:|")
        for res in results
            R = res.R
            ncomp = length(R.components)
            if ncomp <= MAX_INLINE
                desc = fmt_rep(R)
            else
                # the two largest components
                desc = join(fmt_term.(sorted_components(R)[1:2]), " + ") * " + …"
            end
            @printf(io, "| %s | %s | %d | %d | %d | 2^%d | %d | %s | %.3f |\n", string(res.g), fmt_labels(res.l),
                    res.dim, res.npos, res.m0, res.E, ncomp, desc, res.time)
        end
    end

    open(joinpath(outdir, "cosh_product_full.md"), "w") do io
        println(io, "# Complete decompositions\n")
        println(io, "χ = 2^{m₀} ∏_{w ∈ Φ⁺(r)} (e^{w/2} + e^{-w/2}); components are written ",
                    "`mult·[Dynkin labels]` with the dimension as a subscript.\n")
        for res in results
            ncomp = length(res.R.components)
            println(io, "## ", res.g, " ", fmt_labels(res.l), " (dim ", res.dim, ")\n")
            if ncomp <= MAX_FULL
                println(io, fmt_rep(res.R), "\n")
            else
                println(io, "$ncomp components (omitted).\n")
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(length(ARGS) >= 1 ? ARGS[1] : @__DIR__, length(ARGS) >= 2 ? ARGS[2] : nothing)
end
