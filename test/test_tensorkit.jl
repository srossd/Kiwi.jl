# Kiwi irreps as TensorKitSectors sectors, checked with TensorKitSectors' own test suite
# (pentagon, hexagon, triangle, unitarity, fusion tensors vs F/A/B/R/Frobenius–Schur data).
using Test
using Kiwi
import TensorKitSectors
using TensorKitSectors: Fsymbol, Nsymbol, dual, unit

include(joinpath(dirname(dirname(pathof(TensorKitSectors))), "test", "testsuite.jl"))

@testset "TensorKitSectors extension" begin
    I = tensorkit_sector(Kiwi.SU(3))
    s = tensorkit_sector(Irrep(Kiwi.SU(3), 1, 1))
    @test s isa I
    @test Irrep(s) == Irrep(Kiwi.SU(3), 1, 1)
    @test Nsymbol(s, s, s) == 2
    @test Fsymbol(s, s, s, s, s, s) ≈ fsymbol(Irrep(s), Irrep(s), Irrep(s), Irrep(s), Irrep(s), Irrep(s))
    @test dual(tensorkit_sector(Irrep(Kiwi.SU(3), 1, 0))) == tensorkit_sector(Irrep(Kiwi.SU(3), 0, 1))
    for g in (Kiwi.SU(3), Kiwi.Sp(4), Kiwi.G_series(2))
        SectorTestSuite.test_sector(tensorkit_sector(g))
    end
end
