using Unitful
using NaturalSort
using GeophysicalModelGenerator
using LaMEM 
using Serialization




export  search_for_phase_properties, get_phase_number, search_for_model_constrains, get_phase, 
get_phase_bool,  get_data_timestep, split_at__to_type, track_point_over_time, deserialize_file


"""
    material_block = search_for_phase_properties(dat_file::String=dat_path, model_name::String = model_name, search_name::String="<MaterialStart>",  stop_name::String="<MaterialEnd>")

Parameters
====
- `dat_file` - path to the .dat file in the model folder
- `model_name`  - name of the model
- `search_name` -  searching string in dat file
- `stop_name` -  stopping string  in dat file


Examples
========

```julia

julia> dat_path = ("test/test_files/Subduction_VEP.dat")
julia> model_name = "VEP"
julia> material_block = Dict{String, Dict{String, Dict{String, String}}}()  # Dictionary to store material properties
julia> material_block = search_for_phase_properties(dat_path, model_name, "<MaterialStart>", "<MaterialEnd>")

Dict{String, Dict{String, Dict{String, String}}} with 1 entry:
  "VEP" => Dict("Phase5"=>Dict("ch"=>"20e6  # cohesion [Pa]", "Cp"=>"1.2e3 # heat capacity", "fr"=>"30    # friction angle [deg]", "k"=>"2.5", "alpha"=>"1e-5", "disl_prof"=>"Dry_Olivine_disl_creep…

 ```
"""



# Search for a phase in a file and store the properties

function search_for_phase_properties(dat_file::String, model_name::String, search_name::String, stop_name::String)
    number_lines = Int[]  # To store line numbers where search_name is found
    stop_lines = Int[]    # To store line numbers where stop_name is found


    # Open the file and read lines
    lines = String[]
    open(dat_file, "r") do file
        lines = readlines(file)  # Read all lines at once
    end

    # Search for 'search_name' and 'stop_name' in the file
    for i in eachindex(lines)
        if occursin(search_name, lines[i])
            push!(number_lines, i)
        end
        if occursin(stop_name, lines[i])
            push!(stop_lines, i)
        end
    end

    # Initialize material_block entry for the folder if it doesn't exist
    if !haskey(material_block, model_name)
        material_block[model_name] = Dict()
    end

    for k in eachindex(number_lines)
        # Create a phase entry within the folder
        phase_key = "Phase$(k-1)"  # e.g., "Phase0", "Phase1", etc.
        if !haskey(material_block[model_name], phase_key)
            material_block[model_name][phase_key] = Dict()  # Initialize phase in the folder
        end

        # Iterate through the lines between number_lines[k] and stop_lines[k]
        for i in number_lines[k]+1:stop_lines[k]-1
            line = strip(lines[i])  # Clean the line by removing extra spaces
            if occursin("=", line)  # Make sure the line contains an '=' sign
                parts = split(line, "=")  # Split at the '=' symbol
                if length(parts) == 2
                    property_name = strip(parts[1])  # Get the property name
                    property_value = strip(parts[2])  # Get the property value
                    # Store in dictionary, ensure we're adding properties correctly
                    material_block[model_name][phase_key][property_name] = property_value
                end
            end
        end
    end

    return material_block
end

###########################################################################################


"""
  search_value =   search_for_model_constrains(path::String, search_name::String)

Parameters
====
- `path` - path to the ascii file in the model folder
- `search_name` -  searching string in the ascii file


Examples
========

```julia

julia> dat_path = ("test/test_files/Subduction_VEP.dat")
julia> search_name = "surf_level" 
julia> surface_level = search_for_model_constrains(dat_path, search_name)

1-element Vector{Any}:
 0.0

 ```
"""

# search for values set in the julia file to create and are written to ascii files. 

function search_for_model_constrains(path::String, search_name::String)
    number_lines = Int[]  # To store line numbers where search_name is found
    model_constrain=[]

    # Open the file and read lines
    lines = open(path, "r") do file
        readlines(file)  # Read all lines at once
    end

    # Search for 'search_name' and 'stop_name' in the file
    for i in eachindex(lines)
        if occursin(search_name, lines[i])
            push!(number_lines, i)
            con_l = [m.match for m in eachmatch(r"(\d+\.?\d*)", lines[i])]
            for num in con_l
                push!(model_constrain, parse(Float64, num))
            end
            break
        end
    end

    return model_constrain
end

##############################################################################


"""
  phase_coord =   get_phase(path::String,FileName_pvtr::String,phaseIDs::Vector{Int},sep_ind=false)

Parameters
====
- `path` - path to pvtr file
- `FileName_pvtr` -  file name if the pvtr file
- `phaseIDs` -  number of phase which should be detected
- `sep_ind` - false: if each phase coordinates should be stored separately or true: if all coordinates in one vector 


Examples
========

```julia

julia> path = ("test/test_files/timestep/")
julia> FileName_pvtr = "output.pvtr"
julia> phase_coord =   get_phase(path,FileName_pvtr,[2],false)


343-element Vector{CartesianIndex{3}}:
 CartesianIndex(10, 10, 10)
 CartesianIndex(11, 10, 10)
 CartesianIndex(12, 10, 10)
 CartesianIndex(13, 10, 10)
 CartesianIndex(14, 10, 10)

 ```
"""



# get coordinates of one specific phase

function get_phase(path::String,FileName_pvtr::String,phaseID::Int)

    indices = []

    # processing folder
    proc_folder = replace(path,"\\" => "/")*"/"
    println("Processing folder:" * proc_folder)

    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr;fields="phase")
    
    # Get the indices of the phase in the current row
    ind = findall(x -> x == Float64(phaseID), data.fields.phase)
    indices = push!(indices, ind)

    return indices
end

###########################################################################################################################################

# get Coordinates of the searched phases. Return can either be in one list or separated by phase than as vector{vectors}
function get_phase(path::String,FileName_pvtr::String,phaseIDs::Vector{Int},sep_ind=false)

    indices = Int64[]
    idx = Int64[]

    if length(phaseIDs) == 1
        sep_ind=false
        ind = get_phase(path,FileName_pvtr,phaseIDs[1])
        indices = ind[1]
    else
        # Get the indices of the phase in the current row
        for i in eachindex(phaseIDs)

            ind = get_phase(path,FileName_pvtr,i)
            indices = push!(indices, ind)

            if sep_ind && i > 1
                idx = vcat(indices[i-1], indices[i])
            end

        end
    end

    if sep_ind
        indices = idx[1]
    end

    return indices
end

####################################################################################

"""
  matrix =   get_phase_bool(path::String,FileName_pvtr::String,ind)

Parameters
====
- `path` - path to pvtr file
- `FileName_pvtr` -  file name if the pvtr file
- `ind` -  CartesianIndex Points where the phase is detected


Examples
========

```julia

julia> path = ("test/test_files/timestep/")
julia> FileName_pvtr = "output.pvtr"
julia> phase_coord =   get_phase(path,FileName_pvtr,[2],false)

julia> matrix = get_phase_bool(path,FileName_pvtr,phase_coord)

33×33×33 Array{Int64, 3}:
[:, :, 1] =
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0 ....

 ```
"""

# boolse matrix, where coordinates of phase exist values are set to true
function get_phase_bool(path::String,FileName_pvtr::String,ind)

    idx = Int64[]
    proc_folder = replace(path,"\\" => "/")*"/"
    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr)
    if any(x -> isa(x, Vector), ind)
        for i in eachindex(ind)
            if i > 1
                idx = vcat(ind[i-1], ind[i])
            end
        end
        ind = idx[1]
    end

    matrix =zeros(size(data.x.val))
    @views matrix[ind[:]] .= 1
    matrix = Int64.(matrix)

    return matrix
end


##############################################################################


"""
  data =   get_data_timestep(model_path::String,timestep::String,FileName_pvtr::String,p_fields::Vector{String},surface_level::Vector{Any},print=false)

Parameters
====
- `model_path` - path to time steps 
- `timestep` -   time step folder name 
- `FileName_pvtr` -  file name of the pvtr file
- `p_fields` -  property fields which should be read
- `surface_level` - get surface level for depth correction
- `print` - if true prints current processing time step 


Examples
========

```julia

julia> model_path = "test/test_files/timestep/"
julia> dat_path = "test/test_files/Subduction_VEP.dat"
julia> timestep = "Timestep_00000000_0.00000000e+00"
julia> FileName_pvtr = "output.pvtr"
julia> p_fields = ["phase", "temperature"]
julia> surface_level = search_for_model_constrains(dat_path, "surf_level")

julia> data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level,print=false)

CartData 
    size    : (33, 33, 33)
    x       ϵ [ 0.0 : 1.0]
    y       ϵ [ 0.0 : 1.0]
    z       ϵ [ 0.0 : 1.0]
    fields  : (:phase, :temperature)

 ```
"""


# reads data of one time step and corrects the surface level

function get_data_timestep(model_path::String,timestep::String,FileName_pvtr::String,p_fields::Vector{String},surface_level::Vector{Any},print=false)
    # processing folder
    processing_folder = joinpath(model_path,timestep)
    proc_folder = replace(processing_folder,"\\" => "/")*"/"

    print && println("Processing folder:" * proc_folder)

    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr;fields=p_fields)

    # Correct surface level
    data.z.val .-= surface_level # surf_level read from dat file 

    return data
end


############################################################################################

"""
  split_entry =   split_at__to_type(string_to_split::Vector{String} = ["Timestep_00000000_0.00000000e+00"],i::Int64 = 2,type::String="Float")

Parameters
====
- `string_to_split` - string which will be split at _
- `i` -  position within the string which will be returned
- `type` - select type how it will be return: Int, Float, or String 


Examples
========

```julia

julia> string_to_split = ["Timestep_00000000_0.00000000e+00"]
julia> split_value = split_at__to_type(string_to_split,2,"Float")

1-element Vector{Float64}:
 0.0

 ```
"""

# extract time from time_files

function split_at__to_type(string_to_split::Vector{String},i::Int64,type::String)

    spl = if occursin("String", type)
       [split(entry, '_')[i] for entry in string_to_split]
    elseif occursin("Float", type)
       [parse(Float64, split(entry, '_')[i]) for entry in string_to_split]
    elseif occursin("Int", type)
       [parse(Int64, split(entry, '_')[i]) for entry in string_to_split]
    end
    return spl

end




####################################################################################################################

"""
  data =    track_point_over_time(Point_coord::CartesianIndex,fields::Vector{String},model_name::String,timefile_location::String,surface_level::Vector, name::String, output_dir::String,save = false)


Parameters
====
- `Point_coord` - tracked point in CartesianIndex
- `fields` -  property fields which should be tracked
- `timefile_location` -  path of the time step
- `surface_level` - get surface level for depth correction
- `name` -  saving file name
- `output_dir` - save location
- `save` - true: saved, false: not saved


Examples  
========

```julia

julia> Point_coord = CartesianIndex(10, 10, 10)
julia> model_name = "PEV"
julia> timefile_location = "test/test_files/timestep/"
julia> timestep = "Timestep_00000000_0.00000000e+00"
julia> p_fields = ["phase", "temperature"]
julia> surface_level = search_for_model_constrains(dat_path, "surf_level")
julia> name = "track"*string(CartesianIndex(200,1,200))
julia> output_dir = "test/test_files/timestep/"

julia> track_point = track_point_over_time(Point_coord,p_fields,model_name,timefile_location,surface_level, name::String, output_dir::String,save = false)

Dict{String, Dict{String, Float64}} with 2 entries:
  "Timestep_00000010_1.49079279e+01" => Dict("phase"=>1.53445, "temperature"=>0.0)
  "Timestep_00000000_0.00000000e+00" => Dict("phase"=>2.0, "temperature"=>0.0)
   ...
 ```
"""



# track one specific point over time for multiple fields

function track_point_over_time(Point_coord::CartesianIndex,fields::Vector{String},model_name::String,timefile_location::String,surface_level::Vector, name::String, output_dir::String,save = false)

    track_point = Dict{String,Dict{String,Float64}}()

    time_files = filter(f -> startswith(f, "Time"), readdir(timefile_location))
    save_dict = Dict{String,Dict{String,Dict{String,Float64}}}()

    for timestep in time_files

        data = get_data_timestep(timefile_location,timestep,FileName_pvtr,fields,surface_level,false)
        track_point[timestep] = Dict{String,Float64}()

        for field in fields

            field_name = Symbol(field)
            field_data = getfield(data.fields, field_name)

            track_point[timestep][field] = field_data[Point_coord.I[1],Point_coord.I[2],Point_coord.I[3]]

        end
    end


    if save
        save_dict[model_name] = Dict{String,Dict{String,Float64}}()
        save_dict[model_name] = track_point

        file_name = string(name)*".txt"
        output_name=joinpath(output_dir,file_name)

        # Serialize the vector of structures to a file
        open(output_name, "a") do file
            serialize(file, save_dict)
        end
    end


    return track_point

end


###################################################################################################################

"""
  file =   deserialize_file(output_dir::String,name::String)

Parameters
====
- `output_dir` - directory of the file
- `name` -  name of the file


Examples
========

```julia

julia> output_dir = "./output/
julia> name = "track"*string(CartesianIndex(10,10,10))
julia> file_info = deserialize_file(output_dir,name)

1-element Vector{Any}:
 Dict("VEP" => Dict("Timestep_00000010_1.49079279e+01" => Dict("phase" => 1.534447431564331, "temperature" => 0.0), 
        "Timestep_00000000_0.00000000e+00" => Dict("phase" => 2.0, "temperature" => 0.0)))

 ```
"""



# load and extract information of serialized files

function deserialize_file(output_dir::String,name::String)

    file_name = name * ".txt"
    output_name=joinpath(output_dir,file_name)

    det_info =[]

    open(output_name, "r") do file 
        while !eof(file)
            loaded_detach_instances = deserialize(file)
            push!(det_info, loaded_detach_instances)
        end
    end

    return det_info
end


