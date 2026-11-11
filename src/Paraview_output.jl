# Write paraview (VTK/VTS) files for 3D volumes, cross-sections and lines
using WriteVTK

export write_paraview, movie_paraview


"""
    pvd = write_paraview(DataSet::ParaviewData, filename="test"; PointsData=false, pvd=nothing, time=nothing, directory=nothing, verbose=true)

Writes a structure with `Geodata` to a paraview (or VTK) file. If you have unstructured points (e.g., earthquake data), set `PointsData=true`.
In case you want to create a movie in Paraview, and this is a timestep of that movie you also have to pass `time` and `pvd`

# Example 1: Write a 3D volume 
```julia-repl
julia> Lon,Lat,Depth   =   lonlatdepth_grid(10:20,30:40,(-300:25:0)km);
julia> Data_set        =   GeoData(Lat,Lon,Depth,(Depthdata=Depth,LonData=Lon))  
julia> write_paraview(Data_set, "test_depth3D")
```

# Example 2: Horizontal slice @ given depth
```julia-repl
julia> Lon,Lat,Depth  =   lonlatdepth_grid(10:20,30:40,10km);
julia> Data_set       =   GeoData(Lat,Lon,Depth,(Topography=Depth,))  
julia> write_paraview(Data_set, "test")
```

# Example 3: Case with topography
```julia-repl
julia> Lon,Lat,Depth    =   lonlatdepth_grid(10:20,30:40,10km);
julia> Depth[2:4,2:4,1] .=  25km     
julia> Data_set         =   GeoData(Lat,Lon,Depth,(Topography=Depth,))  
julia> write_paraview(Data_set, "test2")
```

# Example 4: Profile
```julia-repl
julia> Lon,Lat,Depth  =   lonlatdepth_grid(10:20,35,(-300:25:0)km);
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth,Depth=Depth))  
julia> write_paraview(Data_set, "test")
```

# Example 5: Velocity vectors
```julia-repl
julia> Lon,Lat,Depth  =   lonlatdepth_grid(10:20,30:40,10km);
julia> Ve, Vn, Vz     =   ones(size(Depth)), ones(size(Depth))*0.5, zeros(size(Depth));
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth, Velocity=(Ve,Vn,Vz)))
GeoData 
  size  : (11, 11, 1)
  lon   ϵ [ 30.0 - 40.0]
  lat   ϵ [ 10.0 - 20.0]
  depth ϵ [ 10.0 km - 10.0 km]
  fields: (:DataSet, :Velocity)  
julia> write_paraview(Data_set, "test_Velocity")
```

# Example 6: Unconnected points (e.g., earthquake locations)
Note that these points should be 1D vectors.
```julia-repl
julia> Lon,Lat,Depth  =   lonlatdepth_grid(10:5:20,35:2:40,(-300:50:0)km);
julia> Lon=Lon[:]; Lat=Lat[:]; Depth=Depth[:];
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth[:],Depth=Depth*10));  
julia> write_paraview(Data_set, "test_Points", PointsData=true)
```

"""
write_paraview(DataSet::ParaviewData, filename = "test"; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true) =
    write_vtk_grid(filename, DataSet.x.val, DataSet.y.val, DataSet.z.val, DataSet.fields; PointsData, pvd, time, directory, verbose)

# Writes the coordinate arrays `x`, `y`, `z` and `fields` to a VTK file; see `write_paraview`.
function write_vtk_grid(filename, x, y, z, fields; PointsData, pvd, time, directory, verbose)

    # Error checking
    if !(ndims(x) == ndims(y) == ndims(z))
        error("The X/Y/Z or Lon/Lat/Depth arrays should be 3 dimensional")
    end

    # Create directory if required
    if !isnothing(directory)
        mkpath(directory)
        filename = joinpath(directory, filename)    # add directory name to pathname
    end

    # Create VT* file
    if PointsData
        # in case we write a dataset with unconnected points (e.g., GPS data, EQ locations etc.)
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:length(x)]
        vtkfile = vtk_grid(filename, vec(strip_units(x)), vec(strip_units(y)), vec(strip_units(z)), cells)
    else
        # for connected 3D grids, 2D planes or 1D lines
        vtkfile = vtk_grid(filename, strip_units(x), strip_units(y), strip_units(z))
    end

    # Add data fields to VT* file
    add_vtk_fields!(vtkfile, fields, PointsData ? VTKPointData() : nothing)
    outfiles = vtk_save(vtkfile)
    if verbose
        println("Saved file: $(outfiles[1])")
    end
    if !isnothing(pvd)
        # Write movie
        pvd[time] = vtkfile
    end

    return pvd
end

# The numbers of an array, without copying: arrays with units are reinterpreted, others returned as they are.
strip_units(A::AbstractArray{<:Quantity}) = ustrip(A)
strip_units(A) = A

# Adds `fields` to `vtkfile`, with the unit in each name. With `location = nothing` WriteVTK
# infers from the array size whether a field lives on points or cells; with `VTKPointData()`
# or `VTKCellData()` the field is written flattened at that location.
function add_vtk_fields!(vtkfile, fields, location = nothing)
    for (name, field) in pairs(fields)
        if field isa Tuple
            # vector fields are written without units
            unit(field[1][1]) != NoUnits && error("potential error as vector data fields have units; please save them with no units!")
            unit_name = ""
            Data = field
        else
            unit_name = unit(field[1])
            Data = strip_units(field)
        end
        name_with_units = "$name  [$unit_name]"
        if isnothing(location)
            vtkfile[name_with_units] = Data
        else
            vtkfile[name_with_units, location] = Data isa Tuple ? Data : vec(Data)
        end
    end
    return vtkfile
end

# Multiple dispatch such that we can also call the routine with GeoData input:
write_paraview(DataSet::GeoData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true) = write_paraview(convert(ParaviewData, DataSet), filename; PointsData, pvd, time, directory, verbose)

"""
    write_paraview(DataSet::UTMData, filename::Any; PointsData=false, pvd=nothing, time=nothing, directory=nothing, verbose=true) 

Writes a `UTMData` structure to paraview. Note that this data is *not* transformed into an Earth-like framework, but remains cartesian instead. 
"""
function write_paraview(DataSet::UTMData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true)

    return write_vtk_grid(filename, DataSet.EW.val, DataSet.NS.val, DataSet.depth.val, DataSet.fields; PointsData, pvd, time, directory, verbose)
end

"""
    write_paraview(DataSet::CartData, filename::Any; PointsData=false, pvd=nothing, time=nothing, directory=nothing, verbose=true) 

Writes a `CartData` structure to paraview. 
"""
function write_paraview(DataSet::CartData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true)

    return write_vtk_grid(filename, DataSet.x.val, DataSet.y.val, DataSet.z.val, DataSet.fields; PointsData, pvd, time, directory, verbose)
end


"""
    pvd = movie_paraview(; name="Movie", pvd=pvd, Finalize::Bool=false, Initialize::Bool=true)

If you want to make a movie of your data set, you can use this routine to initialize and to finalize the movie-file.
It will create a `*.pvd` file, which you can open in Paraview 

Individual timesteps are added to the movie by passing `pvd` and the time of the timestep to the `write_paraview` routine.

Example
=======

Usually this is used inside a `*.jl` script, as in this pseudo-example:
```julia
movie = movie_paraview(name="Movie", Initialize=true)
for itime=1:10
    name = "test"*string(itime)
    movie = write_paraview(Data, name, pvd=movie, time=itime)
end
movie_paraview(pvd=movie, Finalize=true)
```
"""
function movie_paraview(; name = "Movie", pvd = nothing, Finalize::Bool = false, Initialize::Bool = true)

    if (Initialize) & !(Finalize)
        pvd = paraview_collection(name)
    end
    if Finalize
        vtk_save(pvd)
        println("Saved PVD file")
    end

    return pvd
end


"""
    write_paraview(DataSet::Q1Data, filename="test"; directory=nothing, pvd=nothing, time=nothing, verbose=true)
Writes a `Q1Data` dataset to disk, which has cell and vertex field
"""
function write_paraview(DataSet::Q1Data, filename = "test"; directory = nothing, pvd = nothing, time = nothing, verbose = true)

    # Error checking
    if !(length(size(DataSet.x)) == length(size(DataSet.y)) == length(size(DataSet.z)))
        error("The X/Y/Z should be 3 dimensional")
    end

    # Create directory if required
    if !isnothing(directory)
        mkpath(directory)
        filename = joinpath(directory, filename)    # add directory name to pathname
    end

    # Create VT* file
    vtkfile = vtk_grid(filename, strip_units(DataSet.x.val), strip_units(DataSet.y.val), strip_units(DataSet.z.val))

    # Add vertex and cell data fields to VT* file
    add_vtk_fields!(vtkfile, DataSet.fields)
    add_vtk_fields!(vtkfile, DataSet.cellfields, VTKCellData())


    outfiles = vtk_save(vtkfile)
    if verbose
        println("Saved file: $(outfiles[1])")
    end
    if !isnothing(pvd)
        # Write movie
        pvd[time] = vtkfile
    end

    return pvd
end


"""
    write_paraview(DataSet::FEData, filename="test"; directory=nothing, pvd=nothing, time=nothing, verbose=true)
Writes a `FEData` dataset (general finite element) to disk, which has cell and vertex field
"""
function write_paraview(DataSet::FEData, filename = "test"; directory = nothing, pvd = nothing, time = nothing, verbose = true)

    # Create directory if required
    if !isnothing(directory)
        mkpath(directory)
        filename = joinpath(directory, filename)    # add directory name to pathname
    end

    connectivity = DataSet.connectivity
    if size(DataSet.connectivity, 1) == 4
        celltype = VTKCellTypes.VTK_TETRA

    elseif size(DataSet.connectivity, 1) == 8
        celltype = VTKCellTypes.VTK_HEXAHEDRON

        # we need to reorder this as pTatin uses a different ordering than VTK
        id_reorder = [1, 2, 4, 3, 5, 6, 8, 7]
        connectivity = connectivity[id_reorder, :]
    else
        error("This element is not yet implemented")
    end

    # Create VTU file
    points = DataSet.vertices
    cells = MeshCell[]
    for i in 1:size(connectivity, 2)
        push!(cells, MeshCell(celltype, connectivity[:, i]))
    end

    vtkfile = vtk_grid(filename, points, cells)

    # Add vertex and cell data fields to VT* file
    add_vtk_fields!(vtkfile, DataSet.fields)
    add_vtk_fields!(vtkfile, DataSet.cellfields, VTKCellData())

    outfiles = vtk_save(vtkfile)
    if verbose
        println("Saved file: $(outfiles[1])")
    end
    if !isnothing(pvd)
        # Write movie
        pvd[time] = vtkfile
    end

    return pvd
end
