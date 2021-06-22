using GeophysicalModelGenerator
using Test

@testset "Data import.jl" begin
    include("test_data_import.jl")
end
@testset "Data types.jl" begin
    include("test_data_types.jl")
end
@testset "Paraview" begin
    include("test_paraview.jl")
end
@testset "Gravity model" begin
    include("test_voxel_gravity.jl")
end
@testset "Utils" begin
    include("test_utils.jl")
end