module GeophysicalModelGenerator

using Base: String, show_index, Tuple, FieldDescStorage

# Load & export some useful commands/functions from GeoParams:
import GeoParams
using .GeoParams
export
    @u_str, uconvert, upreffered, unit, ustrip, NoUnits, #  Units
    GeoUnit, GEO_units, SI_units, NO_units, AbstractGeoUnits,
    Nondimensionalize, Nondimensionalize!, Dimensionalize, Dimensionalize!,
    superscript, upreferred, GEO, SI, NONE, isDimensional,
    km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, kg, mol,
    isDimensional, Value, NumValue, Unit, UnitValue

export ReadCSV_LatLon, meshgrid, voxel_grav

abstract type AbstractGeneralGrid end                                    # general grid types

export AbstractGeneralGrid

# julia standard library packages
using DelimitedFiles, Statistics

# other packages
using WriteVTK, Colors, MeshIO, FileIO, Interpolations, Geodesy

export vtk_multiblock, vtk_save         # Simplifies writing multiblock files
export LLA
export load

# add files for specific tasks
include("data_types.jl")
include("data_import.jl")
include("nearest_points.jl")
include("utils.jl")
include("Paraview_output.jl")
include("Paraview_collection.jl")
include("transformation.jl")
include("voxel_gravity.jl")
include("LaMEM_io.jl")
include("pTatin_IO.jl")
include("Setup_geometry.jl")
include("stl.jl")
include("ProfileProcessing.jl")
include("IO.jl")
include("IO_ASAGI.jl")
include("event_counts.jl")
include("surface_functions.jl")
include("movies_from_pics.jl")
include("sea_lvl.jl")
include("WaterFlow.jl")

# Add optional routines (only activated when the packages are loaded)

# GMT routines

"""
        import_topo
Optional routine that imports topography. It requires you to load `GMT`
"""
function import_topo end

"""
        import_GeoTIFF
Optional routine that imports GeoTIFF images. It requires you to load `GMT`
"""
function import_GeoTIFF end
export import_topo, import_GeoTIFF

# GLMakie routines

"""
        visualise
Interactive widget that allows you to explore a 3D data set `DataSet` in an interactive manner.
It requires you to load `GLMakie`.
""";
function visualise end
export visualise


"""
    import_Gmsh(fname::String)

Reads a Gmsh file. Requires loading `GridapGmsh`.
"""
function import_Gmsh end
export import_Gmsh


end
