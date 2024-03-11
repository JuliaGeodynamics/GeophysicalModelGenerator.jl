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
@testset "Paraview collection" begin
    include("test_paraview_collection.jl")
end
@testset "Gravity model" begin
    include("test_voxel_gravity.jl")
end
@testset "Nearest points" begin
    include("test_nearest_points.jl")
end
@testset "Utils" begin
    include("test_utils.jl")
end
@testset "Transformations" begin
    include("test_transformation.jl")
end
@testset "Surfaces" begin
    include("test_surfaces.jl")
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

@testset "IO" begin
    include("test_IO.jl")
end

@testset "ProfileProcessing" begin
    include("test_ProfileProcessing.jl")
end

@testset "GMT integration" begin
    include("test_GMT.jl")
end

@testset "Event counts" begin
    include("test_event_counts.jl")
end
@testset "Create movie" begin
    include("test_create_movie.jl")
end

# Cleanup
foreach(rm, filter(endswith(".vts"), readdir()))
foreach(rm, filter(endswith(".vtu"), readdir()))
rm("./markers/",recursive=true)
