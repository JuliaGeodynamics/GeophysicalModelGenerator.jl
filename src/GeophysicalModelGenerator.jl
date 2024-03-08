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
        km, m, cm, mm, Myrs, yr, s, MPa, Pa, Pas, K, C, kg, mol,
        isDimensional, Value, NumValue, Unit, UnitValue

export ReadCSV_LatLon, meshgrid, voxGrav

abstract type AbstractGeneralGrid end                                    # general grid types

export AbstractGeneralGrid

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
include("nearest_points.jl")
include("utils.jl")
include("Paraview_output.jl")
include("Paraview_collection.jl")
include("transformation.jl")
include("voxel_gravity.jl")
include("LaMEM_io.jl")
include("Setup_geometry.jl")
include("stl.jl")
include("ProfileProcessing.jl")
include("IO.jl")
include("event_counts.jl")
include("surface_functions.jl")
include("movies_from_pics.jl")

# Add optional routines (only activated when the packages are loaded)

# GMT routines

"""
        ImportTopo
Optional routine that imports topography. It requires you to load `GMT`
"""
function ImportTopo end

"""
        ImportGeoTIFF
Optional routine that imports GeoTIFF images. It requires you to load `GMT`
"""
function ImportGeoTIFF end
export ImportTopo, ImportGeoTIFF

# GLMakie routines

"""
        Visualise
Interactive widget that allows you to explore a 3D data set `DataSet` in an interactive manner.
It requires you to load `GLMakie`.
"""
function Visualise end
export Visualise



end
