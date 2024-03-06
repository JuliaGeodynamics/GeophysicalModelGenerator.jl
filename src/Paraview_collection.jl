using LightXML
using WriteVTK, Printf

export make_paraview_collection

"""
    make_paraview_collection(; dir=pwd(), pvd_name=nothing, files=nothing, file_extension = ".vts", time = nothing)

In case one has a list of `*.vtk` files, this routine creates a `*.pvd` file that can be opened in Paraview.
This is useful if you previously saved vtk files but didnt save it as a collection in the code itself.

Optional options
===
- `dir`:    directory where the `*.vtk` are stored.
- `pvd_name`:  filename of the resulting `*.pvd` file without extension; if not specified, `full_simulation` is used.
- `files`:  Vector of the `*.vtk` files without extension; if not specified, all `*.vtk` files in the directory are used.
- `file_extension`:  file extension of the vtk files. Default is `.vts` but all `vt*` work.
- `time`:  Vector of the timesteps; if not specified, pseudo time steps are assigned.
"""
function make_paraview_collection(; dir=pwd(), pvd_name=nothing, files=nothing, file_extension = ".vts", time = nothing)

    # if no files are given, use all vtm files in the directory
    curdir = pwd()
    cd(dir)
    if isnothing(files)
        files = filter(endswith(file_extension), readdir())
        files = joinpath.(dir, files)
    end
    # if no time is given, use pseudo time steps
    if isnothing(time)
        time = [string(i) for i in 1:length(files)]
    end
    # Check that the arrays have the same length
    @assert length(time) == length(files)
    @assert length(files) > 0

    # Create a new paraview collection
    if isnothing(pvd_name)
        pvd_name = "full_simulation"
    end

    pvd_name = joinpath(curdir, pvd_name)

    make_paraview_collection(pvd_name, files, time)
    cd(curdir)  # return to original directory
    return pvd_name
end

"""
    make_paraview_collection(pvd_name::String, files::Vector{String}, time::Vector{String})
This function makes a `*.pvd` file from a provided list of `*.vtk` files and time steps.

Inputs are:
- `pvd_name`:  filename of the resulting `*.pvd` file without extension.
- `files`: Vector of desired files to be included in the `*.pvd` file.
- `time`: Vector of time steps corresponding to the `files` vector.
"""

function make_paraview_collection(pvd_name::String, files::Vector{String}, time::Vector{String})
    # Check that the arrays have the same length
    @assert length(time) == length(files)

    # Create a new paraview collection
    pvd = paraview_collection(pvd_name)
    # save the collection
    vtk_save(pvd)

    # Parse the XML file
    xdoc = parse_file("$pvd_name.pvd")

    # Find the <Collection> element
    collection = find_element(root(xdoc), "Collection")

    # Loop over the time and data values
    for i in 1:length(time)
        # Create a new <DataSet> element
        new_dataset = new_child(collection, "DataSet")
        # Set the attributes of the new <DataSet> element
        set_attribute(new_dataset, "timestep", time[i])
        set_attribute(new_dataset, "part", "0")
        set_attribute(new_dataset, "file", files[i])
    end

    # Write the modified XML to a new file of your choice (or overwrite the original)
    open("$pvd_name.pvd", "w") do f
        print(f, xdoc)
    end
    # pretty print the xml file
    str = read("$pvd_name.pvd", String)
    str_out = replace(str, "/>" => "/>\n")
    open("$pvd_name.pvd", "w") do f
        print(f, str_out)
    end
    return pvd_name
end
