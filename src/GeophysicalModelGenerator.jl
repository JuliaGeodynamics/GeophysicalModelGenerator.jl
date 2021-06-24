module GeophysicalModelGenerator

using Base: String, show_index, Tuple, FieldDescStorage


# Load & export some useful commands/functions from GeoParams:
import GeoParams
using .GeoParams
export 
        @u_str, uconvert, upreffered, unit, ustrip, NoUnits,  #  Units 
        GeoUnit, GEO_units, SI_units, NO_units, AbstractGeoUnits, 
        Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
        superscript, upreferred, GEO, SI, NONE, isDimensional, 
        km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, kg, mol

export ReadCSV_LatLon, meshgrid, voxGrav, Screenshot_To_GeoData


# julia standard library packages
using DelimitedFiles            

# other packages
using   WriteVTK, Colors, FileIO, Interpolations, Geodesy

# add files for specific tasks
include("data_types.jl")
include("data_import.jl")
include("coord_conversion.jl")
include("utils.jl")
include("Paraview_output.jl")
include("voxel_gravity.jl")

end
