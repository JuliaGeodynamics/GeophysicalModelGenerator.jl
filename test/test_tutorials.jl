#this runs the tutorials that are part of the JOSS paper to ensure that they keep working moving forward

@testset "Basic tutorial" begin
    include("../tutorials/Tutorial_Basic.jl")
end

#@testset "Jura tutorial" begin
#    include("../tutorials/Tutorial_Jura.jl")
#end

@testset "LaPalma tutorial" begin
    include("../tutorials/Tutorial_LaPalma.jl")
end

# Deactivating this one as it plots in the tutorial
#@testset "AlpineData tutorial" begin
#    include("../tutorials/Tutorial_AlpineData.jl")
#end

@testset "2D Numerical Model tutorial" begin
    include("../tutorials/Tutorial_NumericalModel_2D.jl")
end
#
#@testset "3D Numerical Model tutorial" begin
#    include("../tutorials/Tutorial_NumericalModel_3D.jl")
#end