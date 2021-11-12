# Write paraview (VTK/VTS) files for 3D volumes, cross-sections and lines
using WriteVTK

export Write_Paraview, Movie_Paraview


"""
    pvd = Write_Paraview(DataSet::ParaviewData, filename="test"; PointsData=false, pvd=nothing, time=nothing)

Writes a structure with `Geodata`` to a paraview (or VTK) file. If you have unstructured points (e.g., earthquake data), set `PointsData=true`.
In case you want to create a movie in Paraview, and this is a timestep of that movie you also have to pass `time` and `pvd`

# Example 1: Write a 3D volume 
```julia-repl
julia> Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
julia> Data_set        =   GeoData(Lat,Lon,Depth,(Depthdata=Depth,LonData=Lon))  
julia> Write_Paraview(Data_set, "test_depth3D")
```

# Example 2: Horizontal slice @ given depth
```julia-repl
julia> Lon,Lat,Depth  =   LonLatDepthGrid(10:20,30:40,10km);
julia> Data_set       =   GeoData(Lat,Lon,Depth,(Topography=Depth,))  
julia> Write_Paraview(Data_set, "test")
```

# Example 3: Case with topography
```julia-repl
julia> Lon,Lat,Depth    =   LonLatDepthGrid(10:20,30:40,10km);
julia> Depth[2:4,2:4,1] .=  25km     
julia> Data_set         =   GeoData(Lat,Lon,Depth,(Topography=Depth,))  
julia> Write_Paraview(Data_set, "test2")
```

# Example 4: Profile
```julia-repl
julia> Lon,Lat,Depth  =   LonLatDepthGrid(10:20,35,(-300:25:0)km);
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth,Depth=Depth))  
julia> Write_Paraview(Data_set, "test")
```

# Example 5: Velocity vectors
```julia-repl
julia> Lon,Lat,Depth  =   LonLatDepthGrid(10:20,30:40,10km);
julia> Ve, Vn, Vz     =   ones(size(Depth)), ones(size(Depth))*0.5, zeros(size(Depth));
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth, Velocity=(Ve,Vn,Vz)))
GeoData 
  size  : (11, 11, 1)
  lon   ϵ [ 30.0 - 40.0]
  lat   ϵ [ 10.0 - 20.0]
  depth ϵ [ 10.0 km - 10.0 km]
  fields: (:DataSet, :Velocity)  
julia> Write_Paraview(Data_set, "test_Velocity")
```

# Example 6: Unconnected points (e.g., earthquake locations)
Note that these points should be 1D vectors.
```julia-repl
julia> Lon,Lat,Depth  =   LonLatDepthGrid(10:5:20,35:2:40,(-300:50:0)km);
julia> Lon=Lon[:]; Lat=Lat[:]; Depth=Depth[:];
julia> Data_set       =   GeoData(Lat,Lon,Depth,(DataSet=Depth[:],Depth=Depth*10));  
julia> Write_Paraview(Data_set, "test_Points", PointsData=true)
```



"""
function Write_Paraview(DataSet::ParaviewData, filename="test"; PointsData=false, pvd=nothing, time=nothing)

    # Error checking
    if !(length(size(DataSet.x))==length(size(DataSet.y))==length(size(DataSet.z)))
        error("The X/Y/Z or Lon/Lat/Depth arrays should be 3 dimensional")
    end

    # Create VT* file 
    if PointsData    
        # in case we write a dataset with unconnected points (e.g., GPS data, EQ locations etc.)
        npoints =   length(DataSet.x)
        cells   =   [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]
        x       =   ustrip.(DataSet.x.val);  x = x[:];
        y       =   ustrip.(DataSet.y.val);  y = y[:];
        z       =   ustrip.(DataSet.z.val);  z = z[:]; 

        vtkfile =   vtk_grid(filename, x,y,z, cells)

    else
        # for connected 3D grids, 2D planes or 1D lines
        vtkfile =   vtk_grid(filename, ustrip.(DataSet.x.val), ustrip.(DataSet.y.val), ustrip.(DataSet.z.val)) 
    end

    # Add data fields to VT* file
    names       =   String.(collect(keys(DataSet.fields))); # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)
      
        if typeof(DataSet.fields[index])<: Tuple
            # if we do a tuple of velocities, it appears difficult to deal with units
            # This will require some more work
            unit_name = ""
            Data       =    DataSet.fields[index]  
            if unit(Data[1][1])!=NoUnits
                error("potential error as vector data fields have units; please save them with no units!")
            end
        else
            unit_name = unit(DataSet.fields[index][1])
            Data      = ustrip.(DataSet.fields[index])
        end
      
        name_with_units             = join([name,"  [$(unit_name)]"]); # add units to the name of the field
        if PointsData    
            vtkfile[name_with_units, VTKPointData()]    = Data[:];                        
        else
            vtkfile[name_with_units]    = Data;                        
        end
    end
    outfiles = vtk_save(vtkfile);
    println("Saved file: $(outfiles[1])")

    if !isnothing(pvd)
        # Write movie 
        pvd[time] = vtkfile
    end

    return pvd
end

# Multiple dispatch such that we can also call the routine with GeoData input:
Write_Paraview(DataSet::GeoData,  filename::Any; PointsData=false, pvd=nothing, time=nothing) = Write_Paraview(convert(ParaviewData,DataSet), filename, PointsData=PointsData, pvd=pvd, time=time);

"""
    Write_Paraview(DataSet::UTMData, filename::Any; PointsData=false) 

Writes a `UTMData` structure to paraview. Note that this data is *not* transformed into an Earth-like framework, but remains cartesian instead. 
"""
function Write_Paraview(DataSet::UTMData, filename::Any; PointsData=false, pvd=nothing, time=nothing) 
    
    PVData = ParaviewData(DataSet.EW, DataSet.NS, DataSet.depth.val, DataSet.fields)

    outfiles = Write_Paraview(PVData, filename, PointsData=PointsData, pvd=pvd, time=time);
    return outfiles
end

"""
    Write_Paraview(DataSet::CartData, filename::Any; PointsData=false) 

Writes a `CartData` structure to paraview. 
"""
function Write_Paraview(DataSet::CartData, filename::Any; PointsData=false, pvd=nothing, time=nothing) 
    
    PVData = ParaviewData(DataSet.x.val, DataSet.y.val, DataSet.z.val, DataSet.fields)

    outfiles = Write_Paraview(PVData, filename, PointsData=PointsData, pvd=pvd, time=time);
    return outfiles
end



"""
    pvd = Movie_Paraview(; name="Movie", pvd=pvd, Finalize::Bool=false, Initialize::Bool=true)

If you want to make a movie of your data set, you can use this routine to initialize and to finalize the movie-file.
It will create a `*.pvd` file, which you can open in Paraview 

Individual timesteps are added to the movie by passing `pvd` and the time of the timestep to the `Write_Paraview` routine.

Example
=======

Usually this is used inside a `*.jl` script, as in this pseudo-example:
```julia
movie = Movie_Paraview(name="Movie", Initialize=true)
for itime=1:10
    name = "test"*string(itime)
    movie = Write_Paraview(Data, name, pvd=movie, time=itime)
end
Movie_Paraview(pvd=movie, Finalize=true)
```
"""
function Movie_Paraview(; name="Movie", pvd=nothing, Finalize::Bool=false, Initialize::Bool=true)

    if (Initialize) & !(Finalize)
        pvd = paraview_collection(name) 
    end
    if Finalize
        vtk_save(pvd)
        println("Saved PVD file")
    end
    
    return pvd
end

    
   