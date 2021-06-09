# Write paraview (VTK/VTS) files for 3D volumes, cross-sections and lines
using WriteVTK

export Write_Paraview


"""
    Write_Paraview(DataSet::CartData, filename="test")

Writes a structure with Geodata to a paraview (or VTK) file

# Example 1: Write a 3D volume 
```julia-repl
julia> 
```

# Example 2: horizontal slice @ given depth

# Example 3: profile


"""
function Write_Paraview(DataSet::CartData, filename="test")

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

