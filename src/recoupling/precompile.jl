# Precompilation workload: compiles the recoupling layer for Float64, Rational{BigInt}
# and the reference prime field, so that first calls do not pay for compilation.
# All caches filled here are emptied again, so nothing computed is stored in the image.

@setup_workload begin
    @compile_workload begin
        g = SU(3)
        q, qb, a = Irrep(g, 1, 0), Irrep(g, 0, 1), Irrep(g, 1, 1)
        for T in (Float64, Rational{BigInt}, Fp{_P61})
            fsymbol(a, a, a, a, a, a; T, normalized = false)
            fsymbol(q, q, qb, q, qb, a; T, normalized = false)
            rsymbol(a, a, a; T, normalized = false)
            fusion_tensor(q, qb, a; T, normalized = false)
            twoj(q; T, normalized = false)
            asymbol(q, a, q; T, normalized = false)
            bsymbol(q, a, q; T, normalized = false)
        end
        fsymbol(a, a, a, a, a, a)
        fsymbol(q, q, qb, q, qb, a; T = Rational{BigInt})
        rsymbol(a, a, a)
        fusion_tensor(a, a, a)
        fusion_tensor(a, a, a; sparse = true)
        twoj(a); onej(a); threej(q, q, q); asymbol(q, a, q); bsymbol(q, a, q)
        generators(a)
        fusion_tensor(q, qb, a; basis = :gelfand_tsetlin)
        s = Irrep(SU(2), 1)
        fsymbol(s, s, s, s, Irrep(SU(2), 0), Irrep(SU(2), 0))
        clear_recoupling_cache!()
        lock(() -> empty!(_NSYMBOL_CACHE), _NSYMBOL_LOCK)
        lock(() -> empty!(_DOMCHAR_CACHE), _DOMCHAR_LOCK)
        lock(() -> empty!(_ALGEBRA_DATA), _ALGEBRA_DATA_LOCK)
        empty!(_MN_CACHE)
        empty!(_IRREPS_CACHE)
    end
end
