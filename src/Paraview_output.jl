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
function write_paraview(DataSet::ParaviewData, filename = "test"; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true)

    # Error checking
    if !(length(size(DataSet.x)) == length(size(DataSet.y)) == length(size(DataSet.z)))
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
        npoints = length(DataSet.x)
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:npoints]
        x = ustrip.(DataSet.x.val);  x = x[:]
        y = ustrip.(DataSet.y.val);  y = y[:]
        z = ustrip.(DataSet.z.val);  z = z[:]

        vtkfile = vtk_grid(filename, x, y, z, cells)

    else
        # for connected 3D grids, 2D planes or 1D lines
        vtkfile = vtk_grid(filename, ustrip.(DataSet.x.val), ustrip.(DataSet.y.val), ustrip.(DataSet.z.val))
    end

    # Add data fields to VT* file
    names = String.(collect(keys(DataSet.fields)))  # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)

        if typeof(DataSet.fields[index]) <: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data = DataSet.fields[index]
            if unit(Data[1][1]) != NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.fields[index][1])
            Data = ustrip.(DataSet.fields[index])
        end

        name_with_units = join([name, "  [$(unit_name)]"])  # add units to the name of the field
        if PointsData
            vtkfile[name_with_units, VTKPointData()] = Data[:]
        else
            vtkfile[name_with_units] = Data
        end
    end
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

# Multiple dispatch such that we can also call the routine with GeoData input:
write_paraview(DataSet::GeoData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing) = write_paraview(convert(ParaviewData, DataSet), filename, PointsData = PointsData, pvd = pvd, time = time, directory = directory);

"""
    write_paraview(DataSet::UTMData, filename::Any; PointsData=false, pvd=nothing, time=nothing, directory=nothing, verbose=true) 

Writes a `UTMData` structure to paraview. Note that this data is *not* transformed into an Earth-like framework, but remains cartesian instead. 
"""
function write_paraview(DataSet::UTMData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true)

    PVData = ParaviewData(DataSet.EW, DataSet.NS, DataSet.depth.val, DataSet.fields)

    outfiles = write_paraview(PVData, filename, PointsData = PointsData, pvd = pvd, time = time, directory = directory, verbose = verbose)
    return outfiles
end

"""
    write_paraview(DataSet::CartData, filename::Any; PointsData=false, pvd=nothing, time=nothing, directory=nothing, verbose=true) 

Writes a `CartData` structure to paraview. 
"""
function write_paraview(DataSet::CartData, filename::Any; PointsData = false, pvd = nothing, time = nothing, directory = nothing, verbose = true)

    PVData = ParaviewData(DataSet.x.val, DataSet.y.val, DataSet.z.val, DataSet.fields)

    outfiles = write_paraview(PVData, filename, PointsData = PointsData, pvd = pvd, time = time, directory = directory, verbose = verbose)
    return outfiles
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
    vtkfile = vtk_grid(filename, ustrip.(DataSet.x.val), ustrip.(DataSet.y.val), ustrip.(DataSet.z.val))

    # Add vertex data fields to VT* file
    names = String.(collect(keys(DataSet.fields)))  # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)

        if typeof(DataSet.fields[index]) <: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data = DataSet.fields[index]
            if unit(Data[1][1]) != NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.fields[index][1])
            Data = ustrip.(DataSet.fields[index])
        end

        name_with_units = join([name, "  [$(unit_name)]"])  # add units to the name of the field
        #if PointsData
        #    vtkfile[name_with_units, VTKPointData()]    = Data[:];
        #else
        vtkfile[name_with_units] = Data
        #end
    end

    # Add cell data fields to VT* file
    names = String.(collect(keys(DataSet.cellfields)))  # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)

        if typeof(DataSet.cellfields[index]) <: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data = DataSet.cellfields[index]
            if unit(Data[1][1]) != NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.cellfields[index][1])
            Data = ustrip.(DataSet.cellfields[index])
        end

        name_with_units = join([name, "  [$(unit_name)]"])  # add units to the name of the field
        vtkfile[name_with_units, VTKCellData()] = Data[:]

    end


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

    # Add vertex data fields to VT* file
    names = String.(collect(keys(DataSet.fields)))  # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)

        if typeof(DataSet.fields[index]) <: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data = DataSet.fields[index]
            if unit(Data[1][1]) != NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.fields[index][1])
            Data = ustrip.(DataSet.fields[index])
        end

        name_with_units = join([name, "  [$(unit_name)]"])  # add units to the name of the field
        vtkfile[name_with_units] = Data
    end

    # Add cell data fields to VT* file
    names = String.(collect(keys(DataSet.cellfields)))  # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)

        if typeof(DataSet.cellfields[index]) <: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data = DataSet.cellfields[index]
            if unit(Data[1][1]) != NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.cellfields[index][1])
            Data = ustrip.(DataSet.cellfields[index])
        end

        name_with_units = join([name, "  [$(unit_name)]"])  # add units to the name of the field
        vtkfile[name_with_units, VTKCellData()] = Data[:]
    end

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
