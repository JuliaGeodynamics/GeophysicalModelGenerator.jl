module GeophysicalModelGenerator

# Load & export some useful commands/functions from GeoParams:
import GeoParams
using .GeoParams
export
    @u_str, uconvert, unit, ustrip, NoUnits, #  Units
    GeoUnit, GEO_units, SI_units, NO_units,
    superscript, upreferred, GEO, SI, NONE, isDimensional,
    km, m, cm, mm, yr, s, MPa, Pa, Pas, K, C, kg, mol,
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
include("event_counts.jl")
include("surface_functions.jl")
include("movies_from_pics.jl")
include("sea_lvl.jl")
include("WaterFlow.jl")
include("import_topo.jl")
export import_topo

# Add optional routines (only activated when the packages are loaded)

# NCDatasets routines

"""
    write_ASAGI(fname::String, Data::CartData; fields::Union{Nothing, Tuple}=nothing, km_to_m::Bool=false)

Writes a CartData structure `Data` to the ASAGI file `fname * "_ASAGI.nc"`, which can be read by SeisSol or ExaHype, and
returns the file name. Only individual (scalar) array fields can be written, so vector or tensor fields need to be split first.
The element type of each field is kept.

Keywords:
- `fields`  : tuple with the names of the fields to write, for example `(:Rho, :Sxx)` (default: all fields)
- `km_to_m` : if true, the coordinates are multiplied by 1000 (default=false)

Requires loading `NCDatasets`.
"""
function write_ASAGI end
export write_ASAGI

"""
    data::CartData = read_ASAGI(fname_asagi::String)

Reads a 3D ASAGI NetCDF file, which is used as input for a number of codes such as SeisSol, and returns a `CartData` set.
`fname_asagi` must end in `.nc`. Requires loading `NCDatasets`.
"""
function read_ASAGI end
export read_ASAGI

"""
    data::GeoData = tomo_2_GeoData(filename::String; vel_type::String="vs")

Reads a seismic tomography dataset from a NetCDF file and returns a `GeoData` set. The file must contain variables for
depth, longitude, latitude, and the velocity type `vel_type` (default `"vs"`, the shear wave velocity).
Velocities that the file marks as missing (its fill value) are returned as `NaN`. Requires loading `NCDatasets`.
"""
function tomo_2_GeoData end
export tomo_2_GeoData

# GMT routines

"""
    data_GMT = import_GeoTIFF(fname::String; fieldname=:layer1, negative=false, iskm=true, NorthernHemisphere=true, constantDepth=false, removeNaN_z=false, removeNaN_field=false)

Imports a GeoTIFF dataset (usually containing a surface of some sort) using GMT and returns it as a `GeoData` or `UTMData` set.
The file should have either `UTM` or `longlat` coordinates. If it does not, use QGIS to convert it to `longlat` coordinates.

Keywords:
- `fieldname` : name of the field (default=:layer1)
- `negative`  : if true, the depth is multiplied by -1 (default=false)
- `iskm`      : if true, the depth is multiplied by 1e-3 (default=true)
- `NorthernHemisphere`: if true, the UTM zone is set to be in the northern hemisphere (default=true); only relevant if the data uses UTM projection
- `constantDepth`: if true, the surface is not warped by the z-values, but a constant value is used instead (default=false)
- `removeNaN_z`  : if true, NaN values are removed from the z-dataset (default=false)
- `removeNaN_field`: if true, NaN values are removed from the field (default=false)

Requires loading `GMT`.
"""
function import_GeoTIFF end
export import_GeoTIFF

# GLMakie routines

"""
    fig = visualise(DataSet; Topography=nothing, Topo_range=nothing)

Starts an interactive widget that allows you to explore the 3D data set `DataSet` (a `GeoData` set or an orthogonal `CartData` set)
and returns the `Figure`, which can be saved with `save("file.png", fig)`. All fields of the data set can be explored.
`CartData` sets that are not on an orthogonal grid need to be projected first with `project_CartData`.
The isosurface is off by default and can be toggled in the widget. Vector fields (tuples) are not offered.

Keywords:
- `Topography` : optional topography data set that is drawn on top of `DataSet` (default=nothing)
- `Topo_range` : optional range of the topography color scale (default: symmetric about 0, ± the largest absolute topography)

Requires loading `GLMakie`.
"""
function visualise end
export visualise

# JustPIC routines

"""
    pg = particle_grid(particles::Particles; length_unit)

Returns the coordinates of the `JustPIC` particles, in km, that the geometry routines (`add_box!`, `add_sphere!`, ...)
accept in place of a grid. `length_unit` is the length of one unit of the particle coordinates, for example `1m`, or
the characteristic length of a nondimensional model.

On the CPU, `pg` is a live view: it reads the particle coordinates and active flags when it is used, so it follows later
changes to `particles`. On a GPU, `pg` is a snapshot of the host copies of the coordinates made when `particle_grid` is
called, so it must be recreated after the particles move. The geometry routines process particle arrays on a GPU in host
copies and write the result back.

Requires loading `JustPIC`.
"""
function particle_grid end
export particle_grid

# Gmsh routines

"""
    fe_data, tag_names = import_Gmsh(fname::String)

Reads the 3D Gmsh file `fname` and returns a `FEData` object with the mesh vertices, the cell connectivity and the
region tag of each cell (`fe_data.cellfields.regions`), together with `tag_names`, the names of the regions of the mesh.
2D meshes are not supported. Requires loading `GridapGmsh`.
"""
function import_Gmsh end
export import_Gmsh

end
