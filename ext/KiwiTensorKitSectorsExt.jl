module KiwiTensorKitSectorsExt

# Kiwi irreps as TensorKitSectors.jl sectors.  All topological data are Kiwi's,
# in the lowering basis with Float64 scalars; see `recoupling/symbols.jl` for the
# conventions (which follow TensorKitSectors).

using Kiwi
using Kiwi: Irrep, LieAlgebra
using TensorKitSectors
using TensorKitSectors: SectorValues, SectorSet, GenericFusion, Bosonic
import TensorKitSectors: Nsymbol, Fsymbol, Rsymbol, dim, fusiontensor, dual, unit, ⊗,
                         FusionStyle, BraidingStyle, sectorscalartype, fusionscalartype,
                         braidingscalartype, frobenius_schur_indicator, frobenius_schur_phase,
                         findindex, type_repr

"""
    LieIrrep{S, R} <: Sector

Irrep of the simple Lie algebra `LieAlgebra(S, R)` with Dynkin labels `labels`.
"""
struct LieIrrep{S, R} <: Sector
    labels::NTuple{R, Int}
end

_algebra(::Type{LieIrrep{S, R}}) where {S, R} = LieAlgebra(S, R)
_irrep(s::LieIrrep{S, R}) where {S, R} = Irrep(LieAlgebra(S, R), collect(s.labels))
_sector(::Type{I}, labels::AbstractVector{Int}) where {I <: LieIrrep} = I(Tuple(labels))

Kiwi.tensorkit_sector(g::LieAlgebra) = LieIrrep{g.series, g.rank}
Kiwi.tensorkit_sector(r::Irrep) = LieIrrep{r.algebra.series, r.algebra.rank}(Tuple(r.dynkin_labels))
Kiwi.Irrep(s::LieIrrep) = _irrep(s)

# printing that parses back (in any module where `Kiwi` is bound)
_typestring(::Type{LieIrrep{S, R}}) where {S, R} = "Kiwi.tensorkit_sector(Kiwi.LieAlgebra(:$S, $R))"
Base.show(io::IO, ::Type{T}) where {T <: LieIrrep} = print(io, _typestring(T))
type_repr(::Type{T}) where {T <: LieIrrep} = _typestring(T)
Base.show(io::IO, s::LieIrrep{S, R}) where {S, R} =
    print(io, "Kiwi.tensorkit_sector(Kiwi.Irrep(Kiwi.LieAlgebra(:$S, $R), ", collect(s.labels), "))")

# enumeration of all irreps (Manhattan order of Dynkin labels)
Base.IteratorSize(::Type{SectorValues{I}}) where {I <: LieIrrep} = Base.IsInfinite()
Base.iterate(v::SectorValues{I}, i::Int = 1) where {I <: LieIrrep} = (v[i], i + 1)
function Base.getindex(::SectorValues{LieIrrep{S, R}}, i::Int) where {S, R}
    sz = ntuple(Returns(typemax(Int)), R)
    J = TensorKitSectors.manhattan_to_multidimensional_index(i, sz)
    return LieIrrep{S, R}(map(x -> x - 1, Tuple(J)))
end
function findindex(::SectorValues{LieIrrep{S, R}}, s::LieIrrep{S, R}) where {S, R}
    sz = ntuple(Returns(typemax(Int)), R)
    return TensorKitSectors.to_manhattan_index(map(x -> x + 1, s.labels), sz)
end
Base.isless(a::I, b::I) where {I <: LieIrrep} =
    isless(findindex(values(I), a), findindex(values(I), b))

unit(::Type{LieIrrep{S, R}}) where {S, R} = LieIrrep{S, R}(ntuple(Returns(0), R))
dual(s::I) where {I <: LieIrrep} = _sector(I, Kiwi.conjugate(_irrep(s)).dynkin_labels)
dim(s::LieIrrep) = Int(Kiwi.dimension(_irrep(s)))

FusionStyle(::Type{<:LieIrrep}) = GenericFusion()
BraidingStyle(::Type{<:LieIrrep}) = Bosonic()
sectorscalartype(::Type{<:LieIrrep}) = Float64
fusionscalartype(::Type{<:LieIrrep}) = Float64
braidingscalartype(::Type{<:LieIrrep}) = Float64

function ⊗(a::I, b::I) where {I <: LieIrrep}
    labels = sort!([I(Tuple(l)) for l in keys(Kiwi._fusion_rules(_irrep(a), _irrep(b)))])
    return SectorSet{I}(labels)
end

Nsymbol(a::I, b::I, c::I) where {I <: LieIrrep} = Kiwi.nsymbol(_irrep(a), _irrep(b), _irrep(c))::Int
Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: LieIrrep} =
    Kiwi.fsymbol(_irrep(a), _irrep(b), _irrep(c), _irrep(d), _irrep(e), _irrep(f))::Array{Float64, 4}
Rsymbol(a::I, b::I, c::I) where {I <: LieIrrep} =
    Kiwi.rsymbol(_irrep(a), _irrep(b), _irrep(c))::Matrix{Float64}
fusiontensor(a::I, b::I, c::I) where {I <: LieIrrep} =
    Kiwi.fusion_tensor(_irrep(a), _irrep(b), _irrep(c))::Array{Float64, 4}
frobenius_schur_indicator(s::LieIrrep) = Kiwi.frobenius_schur_indicator(_irrep(s))::Int
frobenius_schur_phase(s::LieIrrep) = Kiwi.onej(_irrep(s))::Int

end # module
