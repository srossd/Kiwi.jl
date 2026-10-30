using Test
using Kiwi

@testset "Kiwi.jl" begin
    include("test_root_systems.jl")
    include("test_irreps.jl")
    include("test_characters.jl")
    include("test_weyl.jl")
    include("test_tensor_products.jl")
    include("test_symmetric.jl")
    include("test_fast_algorithms.jl")
    include("test_recoupling.jl")
    include("test_tensorkit.jl")
end