module GeophysicalModelGenerator

using Base: String, show_index, Tuple, FieldDescStorage
using Requires


# Load & export some useful commands/functions from GeoParams:
import GeoParams
using .GeoParams
export 
        @u_str, uconvert, upreffered, unit, ustrip, NoUnits,  #  Units 
        GeoUnit, GEO_units, SI_units, NO_units, AbstratGeoUnits, 
        Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
        superscript, upreferred, GEO, SI, NONE, isDimensional, 
        km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, kg, mol,
        isDimensional, Value, NumValue, Unit, UnitValue

export ReadCSV_LatLon, meshgrid, voxGrav

# julia standard library packages
using DelimitedFiles, Statistics            

# other packages
using   WriteVTK, Colors, MeshIO, FileIO, Interpolations, Geodesy

export vtk_multiblock, vtk_save         # Simplifies writing multiblock files
export LLA
export load

# add files for specific tasks
include("data_types.jl")
include("data_import.jl")
include("coord_conversion.jl")
include("utils.jl")
include("Paraview_output.jl")
include("transformation.jl")
include("voxel_gravity.jl")
include("LaMEM_io.jl")
include("LaMEM_geometry.jl")
include("stl.jl")

# Add optional routines (only activated when the packages are loaded)
function __init__()
        @require GMT = "5752ebe1-31b9-557e-87aa-f909b540aa54" begin
                @eval include("./GMT_utils.jl")
        end
        @require GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
                println("Loading GLMakie plotting routines within GMG")
                @eval include("./Visualisation.jl")
        end
        
end



end
