# Hooks for the TensorKitSectors.jl extension (ext/KiwiTensorKitSectorsExt.jl).

"""
    tensorkit_sector(g::LieAlgebra) -> Type
    tensorkit_sector(r::Irrep) -> sector

The TensorKitSectors.jl sector type of the irreps of `g`, and the sector of an irrep.
Its `Nsymbol`, `Fsymbol`, `Rsymbol`, `fusiontensor`, `dim`, `dual` and
`frobenius_schur_phase` are Kiwi's [`nsymbol`](@ref), [`fsymbol`](@ref),
[`rsymbol`](@ref), [`fusion_tensor`](@ref), `dimension`, `conjugate` and
[`onej`](@ref) (lowering basis, `Float64`).  `Irrep(s)` converts back.

Available after `using TensorKitSectors` (or TensorKit), through a package extension.
"""
function tensorkit_sector end
