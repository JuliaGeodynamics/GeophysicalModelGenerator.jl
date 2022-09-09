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
@testset "Transformations" begin
    include("test_transformation.jl")
end

@testset "LaMEM" begin
    include("test_lamem.jl")
end

@testset "SetupGeometry" begin
    include("test_setup_geometry.jl")
end

@testset "STL" begin
    include("test_stl.jl")
end

@testset "Dislocation Models" begin
    include("test_dislocation_models.jl")
end

# Cleanup 
foreach(rm, filter(endswith(".vts"), readdir()))
foreach(rm, filter(endswith(".vtu"), readdir()))
rm("./markers/",recursive=true)

