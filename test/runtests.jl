using GeophysicalModelGenerator
using Test

@testset "GeophysicalModelGenerator.jl" begin
    #include("test_data_import.jl")
    include("test_data_types.jl")
    include("test_paraview.jl")
end
