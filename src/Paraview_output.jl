# Write paraview (VTK/VTS) files for 3D volumes, cross-sections and lines
using WriteVTK

export Write_Paraview


"""
    Write_Paraview(DataSet::CartData, filename="test")

Writes a structure with Geodata to a paraview (or VTK) file

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

"""
function Write_Paraview(DataSet::CartData, filename="test")

    # Error checking
    if !(length(size(DataSet.x))==length(size(DataSet.y))==length(size(DataSet.z)))
        error("The X/Y/Z or Lon/Lat/Depth arrays should be 3 dimensional")
    end

    # Create VT* file 
    vtkfile     =   vtk_grid(filename, ustrip(DataSet.x.val), ustrip(DataSet.y.val), ustrip(DataSet.z.val)) 

    # Add data fields to VT* file
    names       =   String.(collect(keys(DataSet.fields))); # this is how to retrieve the names of the data fields
    for (index, name) in enumerate(names)
        name_with_units             = join([name,"  [$(unit(DataSet.fields[index][1]))]"]); # add units to the name of the field
        vtkfile[name_with_units]    = ustrip(DataSet.fields[index]);                        
    end
    outfiles = vtk_save(vtkfile);

    return outfiles
end

# Multiple dispatch such that we can also call the routine with GeoData input:
Write_Paraview(DataSet::GeoData, filename::Any) = Write_Paraview(convert(CartData,DataSet), filename);

