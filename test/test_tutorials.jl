#this runs the tutorials that are part of the JOSS paper to ensure that they keep working moving forward

@testset "Basic tutorial" begin
    include("../tutorials/Tutorial_Basic.jl")
end

#@testset "Jura tutorial" begin
#    include("../tutorials/Tutorial_Jura.jl")
#end

# This one downloads topography, which topo_prefetch.jl has put in the tile cache. If the
# prefetch could not reach the GMT data server, skip on a push/PR run; the scheduled run
# insists, so a real breakage is caught.
topo_server_ok = get(ENV, "GMG_TOPO_PREFETCH_OK", "true") == "true" ||
    get(ENV, "GITHUB_EVENT_NAME", "") == "schedule"

@testset "LaPalma tutorial" begin
    if !topo_server_ok
        @test_skip "GMT data server unreachable: LaPalma tutorial skipped"
    else
        include("../tutorials/Tutorial_LaPalma.jl")
    end
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
