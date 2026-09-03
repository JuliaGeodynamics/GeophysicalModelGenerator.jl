using Unitful
using NaturalSort
using GeophysicalModelGenerator
using LaMEM 
using JLD2
using CairoMakie
using Printf
using Statistics
using Plots




export  search_for_phase_properties, search_for_model_constrains, search_for_all_model_contrains, get_phase, get_phase_bool, split_at__to_type, 
get_data_timestep, get_tracer_timestep,  get_surf_timestep,
post_plot, find_field_properties_grid, find_general_grid_prop, find_surf_evolution, find_tracer_info, 
track_point_over_time

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

julia> dat_path = ("test/input_files/Subduction_VEP.dat")
julia> search_name = "surf_level" 
julia> surface_level = search_for_model_constrains(dat_path, search_name)

1-element Vector{Any}:
 0.0

 ```
"""

# search for values set in the julia file to create and are written to ascii files. stops at first found
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
            line = first.(split.(lines[i], '#'))
            con_l = [m.match for m in eachmatch(r"(\d+\.?\d*)", line)]
            for num in con_l
                push!(model_constrain, parse(Float64, num))
            end
            break
        end
    end

    return model_constrain
end



"""
  search_value =   search_for_all_model_constrains(path::String, search_name::String)

Parameters
====
- `path` - path to the ascii file in the model folder
- `search_name` -  searching string in the ascii file


Examples
========

```julia

julia> dat_path = ("test/input_files/Subduction_VEP.dat")
julia> search_name = "rho" 
julia> surface_level = search_for_all_model_constrains(dat_path, search_name)

1-element Vector{Any}:
3-element Vector{Any}:
   60.0
 2700.0
 2750.0

    ⋮
 2800.0
 4000.0
 ```
"""

# search for values set in the julia file to create and are written to ascii files. stops at first found
function search_for_all_model_constrains(path::String, search_name::String)
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
            line = first.(split.(lines[i], '#'))
            con_l = [m.match for m in eachmatch(r"(\d+\.?\d*)", line)]
            for num in con_l
                push!(model_constrain, parse(Float64, num))
            end
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

    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr;fields="phase")
    
    # Get the indices of the phase in the current row
    ind = findall(x -> x == Float64(phaseID), data.fields.phase)
    indices = push!(indices, ind)

    return indices
end

###########################################################################################################################################

# get Coordinates of the searched phases. Return can either be in one list or separated by phase than as vector{vectors}
function get_phase(location::String,FileName_pvtr::String,phaseIDs::Vector{Int},sep_ind=false)

    indices = []    # preallocation
    idx = []    # preallocation

    if length(phaseIDs) == 1
        sep_ind=false
        ind = get_phase(location,FileName_pvtr,phaseIDs[1])         # get indices of the phase location in matrix
        indices = ind[1]  
    else
        for i in eachindex(phaseIDs)                                # Get indices for multiple phases
            ind = get_phase(location,FileName_pvtr,phaseIDs[i])     # get indices of the phase location in matrix
            indices = push!(indices, ind)       

            # sep ind --> indices for each phase separately listed. Make from Vector{} to Vector{Vector{}}
            if (sep_ind  == true && i > 1)
                idx = vcat(indices[i-1][1], indices[i][1])
            end

        end
    end

    if sep_ind  == true
        indices = idx
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
function get_phase_bool(location::String,FileName_pvtr::String,ind)

    proc_folder = replace(location,"\\" => "/")*"/"             # correct path for linux
    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr)      # read data from path

    # reduce matrix to Vector{Strings} for the positions
    if ind isa Vector{Any}
        ind = reduce(vcat, Iterators.flatten(ind))
    end

    matrix =zeros(size(data.x.val)) # preallocation
    matrix[ind[:]] .= 1 # set all indices which have the correct phase to 1 at the matrix position
    matrix = Int64.(matrix) # make the matrix to integer
    return matrix
end


##############################################################################


"""
  data =   get_data_timestep(model_path::String,timestep::String,FileName_pvtr::String,p_fields::Vector{String},surface_level::Vector{Any})

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
julia> dat_path = "test/input_files/Subduction_VEP.dat"
julia> timestep = "Timestep_00000000_0.00000000e+00"
julia> FileName_pvtr = "output.pvtr"
julia> p_fields = ["phase", "temperature"]
julia> surface_level = search_for_model_constrains(dat_path, "surf_level")

julia> data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level)

CartData 
    size    : (33, 33, 33)
    x       ϵ [ 0.0 : 1.0]
    y       ϵ [ 0.0 : 1.0]
    z       ϵ [ 0.0 : 1.0]
    fields  : (:phase, :temperature)

 ```
"""


# reads data of one time step and corrects the surface level

function get_data_timestep(model_path::String,timestep::String,FileName_pvtr::String,p_fields::Vector{String},surface_level::Vector{Any})
    # processing folder
    processing_folder = joinpath(model_path,timestep)
    proc_folder = replace(processing_folder,"\\" => "/")*"/"

    data = read_LaMEM_PVTR_file(proc_folder,FileName_pvtr;fields=p_fields)

    # Correct surface level
    data.z.val .-= surface_level # surf_level read from dat file 

    return data
end

###############################################################################################
# get tracer information

"""

tracer = get_tracer_timestep(model_path::String,timestep::String,FileName_pvtu::String)

Parameters
====
- `model_path` - path to time steps 
- `timestep` -   time step folder name 
- `FileName_pvtu` -  file name of the pvtr file
- `print` - if true prints current processing time step 



Examples
========

```julia

julia> model_path = "test/test_files/timestep/"
julia> dat_path = "test/input_files/Passive_tracer_ex2D.dat"
julia> timestep = "Timestep_00000000_0.00000000e+00"
julia> FileName_pvtu = "output.pvtu"
julia> p_fields = ["phase", "temperature"]
julia> surface_level = search_for_model_constrains(dat_path, "surf_level")

julia> tracer = get_tracer_timestep(model_path,timestep,FileName_pvtu)

"""

function get_tracer_timestep(model_path::String,timestep::String,FileName_pvtu::String)

    processing_folder = joinpath(model_path,timestep)    # path to processing folder
    proc_folder = replace(processing_folder,"\\" => "/")*"/"    # correct path for linux

    tracer = read_LaMEM_PVTU_file(proc_folder,FileName_pvtu)    # read tracer pvtu file

    return tracer
end
##############################################################################################################
# get surface information


"""

surf = get_surf_timestep(model_path::String,timestep::String,FileName_pvts::String,surface_level::Vector{Any})

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
julia> dat_path = "test/input_files/Passive_tracer_ex2D.dat"
julia> timestep = "Timestep_00000000_0.00000000e+00"
julia> FileName_pvts = "output.pvts"
julia> p_fields = ["phase", "temperature"]
julia> surface_level = search_for_model_constrains(dat_path, "surf_level")

julia> surf = get_surf_timestep(model_path,timestep,FileName_pvts,surface_level)

"""

function get_surf_timestep(model_path::String,timestep::String,FileName_pvts::String,surface_level::Vector{Any})

    processing_folder = joinpath(model_path,timestep)    # path to processing folder
    proc_folder = replace(processing_folder,"\\" => "/")*"/"    # correct path for linux

    surf = read_LaMEM_PVTS_file(proc_folder,FileName_pvts)  # read surface pvts file

    surf.fields.topography .= surf.fields.topography .- surface_level # Correct surface level

    return surf
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
  split_entry =   split_at__to_type(string_to_split::String = "Timestep_00000000_0.00000000e+00",i::Int64 = 2,type::String="Float")

Parameters
====
- `string_to_split` - string which will be split at _
- `i` -  position within the string which will be returned
- `type` - select type how it will be return: Int, Float, or String 


Examples
========

```julia

julia> string_to_split = "Timestep_00000000_0.00000000e+00"
julia> split_value = split_at__to_type(string_to_split,2,"Float")

Float64:
 0.0

 ```
"""

function split_at__to_type(timestep::String,i::Int64,type::String)
    
    # split time string on_ and convert to wished type
    if occursin("String", type)
        time = split(timestep, '_')[i]
    elseif occursin("Float", type)
        time = parse(Float64, split(timestep, '_')[i]) 
    elseif occursin("Int", type)
        time = parse(Int64, split(timestep, '_')[i])
    end
    return time

end


##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
# save data for postprocessing



"""

        post_plot(data::CartData,p_fields::Vector{String},timestep::String,number_phases::Int64,y_slices::Vector{Int64},output_folder::String,surface_level::Vector,dxdz::Vector{Float64}, textpos::Vector{Float64}, numb_ticks::Int64; phase_contour=false)


Parameters
====
- `data` - CartData structure of the model for one timestep,
- `p_fields` - Fields which should be plotted,
- `timestep` - Timestep name to plot on figure,
- `number_phases` - number of phase to generate the phase colorbar color consistent for all timesteps,
- `y_slices` -select which y slice should be plotted,
- `output_folder` - path to output folder,
- `surface_level` - Elevation of surface in Coordinates, 
- `dxdz` - Which part should be plotted in Coordinates, 
- `textpos` - position of where the timestep information should be printed on the figure,
- `numb_ticks` -number of ticks on the axis
- `phase_contour` - plot the phase contour on top of the filed, bool)

"""


# plot field properties to get first impression of detachment 2D plotting, y slice can be chosen
function post_plot(data::CartData,p_fields::Vector{String},timestep::String,number_phases::Int64,y_slices::Vector{Int64},output_folder::String,surface_level::Vector,dxdz::Vector{Float64}, textpos::Vector{Float64}, numb_ticks::Int64; phase_contour=false)

    # get size of plotting
    x_vec = data.x.val[:,1,1]
    z_vec = data.z.val[1,1,:]

    xposmin = argmin(abs.(x_vec .- dxdz[1]))
    xposmax = argmin(abs.(x_vec .- dxdz[2]))
    zposmin = argmin(abs.(z_vec .- (dxdz[3] .+ surface_level[1])))
    zposmax = argmin(abs.(z_vec .- (dxdz[4] .+ surface_level[1])))

    # get four ticks on each axis
    dx = (x_vec[xposmax] - x_vec[xposmin])/(numb_ticks -1)
    dz = (z_vec[zposmax] - z_vec[zposmin])/(numb_ticks -1)


    # determine required decimals
    digits_x = max(0, -floor(Int, log10(abs(dx))))
    digits_z = max(0, -floor(Int, log10(abs(dz))))

    xticks = round.(collect(x_vec[zposmin]:dx:x_vec[xposmax]),digits=digits_x)
    yticks = round.(collect(z_vec[zposmin]:dz:z_vec[zposmax]),digits=digits_z)


    # set figure size
    aspect = dx / dz
    base = 1200

    if aspect >= 1
        plot_width = base
        plot_height = round(Int, base / aspect)
    else
        plot_height = base
        plot_width = round(Int, base * aspect)
    end

    # prevent absurdly small figures
    plot_width = max(plot_width, 400)
    plot_height = max(plot_height, 300)
    cbar_width = 100
          
    # general figure properties
    size=(plot_width + cbar_width, plot_height)

    xtickfontsize=round(mean(size)*0.025)
    ytickfontsize=round(mean(size)*0.025)
    labelfontsize=round(mean(size)*0.025) +2
    xlabel="Distance [km]"
    ylabel="Depth [km]"

    # title 
    time = round(parse(Float64, string(split(timestep, "_")[3])),digits=3)
    title="t = "*  @sprintf("%.5f",time)* " Myr"
    titlefontsize =round(maximum(size)*0.03)

    # colorbar
    levels=number_phases
    colorbar_tickfontsize=round(mean(size)*0.025)
    thickness_scaling = 2
    cb_width = round(mean(size)*0.05)


    for i in eachindex(p_fields)

        # plot parameters reset for each field
        field_name = Symbol(p_fields[i])
        colorbar_scale= identity
        cb_ticks="automatic"
        ticks_pos = WilkinsonTicks(9)

        ###############################################################################
        # for matrices in each direction --> 9 directions
        if  p_fields[i] == "dev_stress" || p_fields[i] == "vel_gr_tensor"

            # get direction
            data_direction = ["_xx", "_xy", "_xz", "_yx","_yy","_yz", "_zx", "_zy", "_zz"]

            # load specific fields and set colormap
            for j in eachindex(data_direction)

                if "dev_stress" == p_fields[i]
                    cmap_name="deep"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Deviatoric stress tensor []"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                elseif "vel_gr_tensor" == p_fields[i]
                    cmap_name="deep"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Velocity gradient tensor []"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                end


                for s in y_slices
                    fig = CairoMakie.Figure(size = size)

                    ax = Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, aspect = DataAspect(),
                    xticks = xticks, yticks = yticks, xticklabelsize = xtickfontsize, yticklabelsize = ytickfontsize, 
                    xlabelsize = labelfontsize, ylabelsize = labelfontsize)

                    hm = CairoMakie.heatmap!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], field_data[:, s, :], colormap = cmap, colorscale = colorbar_scale)

                    CairoMakie.contour!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], data.fields.phase[xposmin:xposmax, s, zposmin:zposmax]; levels = levels, color = :white, linewidth = 0.5,)
                    CairoMakie.text!(ax,textpos[1], textpos[2]; text = title, align = (:left, :bottom), fontsize = titlefontsize, color = textcolor)
                    cb = CairoMakie.Colorbar(fig[1, 2], hm;label = cbar_label, width = cb_width, labelsize = colorbar_tickfontsize,ticklabelsize = colorbar_tickfontsize)
                    cb.height = ax.scene.viewport[].widths[2]
                    colgap!(fig.layout, 10)
                    ##########################################
                    # make directories and save figs

                    output_field = joinpath(output_folder,string(p_fields[i],data_direction[j]))
                    if !isdir(output_field)
                        mkdir(output_field)
                    end

                    fig_name = "Fig"*string(parse(Int,split(timestep,"_")[2]))*"_"*string(NumValue(data.y[1,s,1]))*".png"
                    output_name=joinpath(output_field,fig_name)
                    CairoMakie.save(output_name,fig)
                end
            end


        #################################################################################
        # for vector fields not in cell center --> three dircetions

        elseif p_fields[i] == "tot_displ" || p_fields[i] =="velocity" || p_fields[i] =="moment_res" || p_fields[i] == "EHmax" || p_fields[i] == "SHmax"

            # get direction
            data_direction = ["_x", "_y", "_z"]
            
            # load specific fields and set colormap
            for j in eachindex(data_direction)

                if "EHmax" == p_fields[i]
                    cmap_name="deep"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Maximum horizontal extension []"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                elseif "SHmax" == p_fields[i]
                    cmap_name="deep"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Maximum horizontal stress []"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                elseif "tot_displ"  == p_fields[i]
                    cmap_name="deep"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Total displacements [km]"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                elseif "velocity"  == p_fields[i]
                    cmap_name="speed"
                    field_data = getfield(data.fields, field_name)[j]
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    cbar_label="Velocity [cm/yr]"
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan

                elseif "moment_res"  == p_fields[i]
                    cmap_name="navia"
                    cbar_label="log₁₀(Momentum residual) [N³]"
                    field_data = log10.(abs.(getfield(data.fields, field_name)[j]))  # Access the field using the symbol
                    field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                    field_data = ifelse.(isfinite.(field_data), field_data, 0)
                    cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                    textcolor = :cyan
                end
                
                # plot field in specific direction
                for s in y_slices
                    fig = CairoMakie.Figure(size = size)

                    ax = Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, aspect = DataAspect(),
                    xticks = xticks, yticks = yticks, xticklabelsize = xtickfontsize, yticklabelsize = ytickfontsize, 
                    xlabelsize = labelfontsize, ylabelsize = labelfontsize)


                    hm = CairoMakie.heatmap!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], field_data[:, s, :], colormap = cmap, colorscale = colorbar_scale)


                    CairoMakie.contour!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], data.fields.phase[xposmin:xposmax, s, zposmin:zposmax]; levels = levels, color = :white, linewidth = 0.5,)
                    CairoMakie.text!(ax,textpos[1], textpos[2]; text = title, align = (:left, :bottom), fontsize = titlefontsize, color = textcolor)
                    cb = CairoMakie.Colorbar(fig[1, 2], hm;label = cbar_label,ticks =  ticks_pos, width = cb_width, labelsize = colorbar_tickfontsize,ticklabelsize = colorbar_tickfontsize)
                    cb.height = ax.scene.viewport[].widths[2]
                    colgap!(fig.layout, 10)
                    ##########################################
                    # make directories and save figs

                    output_field = joinpath(output_folder,string(p_fields[i],data_direction[j]))
                    if !isdir(output_field)
                        mkdir(output_field)
                    end

                    fig_name = "Fig"*string(parse(Int,split(timestep,"_")[2]))*"_"*string(NumValue(data.y[1,s,1]))*".png"
                    output_name=joinpath(output_field,fig_name)
                    CairoMakie.save(output_name,fig)
                end
            end
        else

        ####################################################################################
        # for information on cell center --> read specific fields and create colormap

            if "conductivity"== p_fields[i]
                cbar_label="Conductivity [W/m/K]"
                cmap_name="navia"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name),  30, categorical=true,rev=true)
                textcolor = :cyan

            elseif "cont_res" == p_fields[i]
                cbar_label="log₁₀(Continuity residual) [1/s]"
                cmap_name="navia"
                field_data = log10.(abs.(getfield(data.fields, field_name)))  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                cmap=cgrad(Symbol(cmap_name),  30, categorical=true,rev=true)
                textcolor = :cyan

            elseif "density"  == p_fields[i]
                cbar_label="Density [kg/m³]"
                cmap_name="dense"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_val = sort(unique(round.(Int, field_data./ 100) .* 100))
                cmap=cgrad(Symbol(cmap_name), 20, categorical= true, rev=true)
                textcolor = :magenta

            elseif "energ_res"  == p_fields[i]
                cbar_label="log₁₀(Energy residual)" *"[W³]"
                cmap_name="navia"
                field_data = log10.(abs.(getfield(data.fields, field_name)))  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                cmap=cgrad(Symbol(cmap_name),30,categorical=true, rev=true)
                cb_ticks = WilkinsonTicks(5)
                textcolor = :cyan

            elseif "fluid_density"  == p_fields[i]
                cbar_label="Fluid density [kg/m³]"
                cmap_name="matter"
                field_data = getfield(data.fields, field_name)./1e3  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "j2_dev_stress"  == p_fields[i]
                cbar_label="Deviatoric stress second invariant " *"[MPa]"#*L"$\tau_{II}$"*
                cmap_name="matter"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_val =unique(Int.(round.(unique(filter(!isnan, field_data)),digits=-1)))
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan
                
            elseif "j2_strain_rate"  == p_fields[i]
                cbar_label="log₁₀(Creep Deviatoric strain rate second invariant) s⁻¹" #L"\dot{\epsilon}" 
                cmap_name="default"
                field_data = log10.(abs.(getfield(data.fields, field_name)))  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true)
                textcolor = :cyan

            elseif "litho_press"  == p_fields[i]
                cbar_label="Lithospheric pressure [MPa]"#*L"$^3$"*
                cmap_name="matter"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "melt_fraction"  == p_fields[i]
                cbar_label="Melt fraction []"#*L"$^3$"*
                cmap_name="matter"
                field_data = getfield(data.fields, field_name)./1e3  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "over_press"  == p_fields[i]
                cbar_label="Overpressure [MPa]"#*L"$^3$"*
                cmap_name="matter"
                field_data = getfield(data.fields, field_name) # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "phase"  == p_fields[i]
                cbar_label="Phase"
                cmap_name="plasma"
                field_data = getfield(data.fields, field_name)# Access the field using the symbol
                field_data = Int.(round.(field_data[xposmin:xposmax,:,zposmin:zposmax]))
                field_val = sort(Int.(unique(filter(!isnan, field_data))))
                levels = length(field_val)
                col_map=cgrad(Symbol(cmap_name),number_phases,categorical = true)
                col = collect(col_map)
                col[1] = RGBA(1,1,1,1)
                if number_phases > maximum(field_val) +1
                col = col[field_val .+ 1]
                end
                cmap = PlotUtils.cgrad(col; categorical=true)
                cb_values = field_val #.+ 0.5 #collect(0.5:1:levels-0.5) # Custom tick positions
                cb_labels = string.(field_val)  # Custom labels
                ticks_pos = (cb_values ,cb_labels)
                textcolor = :cyan

            elseif "plast_dissip"  == p_fields[i]
                cbar_label="Plastic dissipation [W/m³]"
                cmap_name="thermal"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_val = sort(Int.(unique(round.(filter(!isnan, field_data)))))
                levels = length(field_val)
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true)
                textcolor = :cyan

            elseif "plast_strain"  == p_fields[i]
                cbar_label="log₁₀(Accumulated plastic strain) []"
                cmap_name="imola"
                field_data = log10.(abs.(getfield(data.fields, field_name)))  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                field_val = unique(round.(Int, field_data))
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "pore_press"  == p_fields[i]
                cbar_label="Pore pressure [MPa]"#*L"$^3$"*
                cmap_name="matter"
                field_data = getfield(data.fields, field_name) # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "pressure"  == p_fields[i]
                cbar_label="Pressure [MPa]"
                cmap_name="matter" # noch durch GPa teilen
                field_data = getfield(data.fields, field_name) # is already in MPa
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "rel_dif_rate"  == p_fields[i]
                cbar_label="Diffusion creep relative strain rate []"
                cmap_name="batlow"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :magenta

            elseif "rel_dis_rate"  == p_fields[i]
                cbar_label="Dislocation creep relative strain rate []"
                cmap_name="batlow"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :magenta

            elseif "rel_prl_rate"  == p_fields[i]
                cbar_label="Relative low-temperature-plasticity creep []"
                cmap_name="batlow"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :magenta

            elseif "rel_pl_rate"  == p_fields[i]
                cbar_label="Peierls creep relative strain rate []"
                cmap_name="batlow"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :magenta

            elseif "temperature"  == p_fields[i]
                cbar_label="Temperature [°C]"
                cmap_name="thermal" # noch durch GPa teilen
                field_data = getfield(data.fields, field_name)
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_val =unique(Int.(round.(unique(filter(!isnan, field_data)) ./50).*50))
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true)
                textcolor = :cyan

            elseif "total_pressure"  == p_fields[i]
                cbar_label="Total pressure [MPa]"
                cmap_name="matter" # noch durch GPa teilen
                field_data = getfield(data.fields, field_name) # is already in MPa
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "visc_creep"  == p_fields[i]
                cbar_label="log₁₀(Creep effective viscosity) [Pa s]"
                cmap_name="magma"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                field_val = unique(round.(field_data, digits=1))
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "visc_total"  == p_fields[i]
                cbar_label="log₁₀(Total effective viscosity) [Pa s]"
                cmap_name="magma"
                field_data = getfield(data.fields, field_name) # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                field_data = ifelse.(isfinite.(field_data), field_data, 0)
                field_val = unique(round.(field_data, digits=1))
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true, rev=true)
                textcolor = :cyan

            elseif "yield"  == p_fields[i]
                cbar_label="Yield stress [MPa]"
                cmap_name="matter"
                field_data = getfield(data.fields, field_name)  # Access the field using the symbol
                field_data = field_data[xposmin:xposmax,:,zposmin:zposmax]
                cmap=cgrad(Symbol(cmap_name), 30, categorical=true)
                textcolor = :cyan

            end
        
            # plot field 
            for s in y_slices
                    fig = CairoMakie.Figure(size = size)

                    ax = Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel, aspect = DataAspect(),
                    xticks = xticks, yticks = yticks, xticklabelsize = xtickfontsize, yticklabelsize = ytickfontsize, 
                    xlabelsize = labelfontsize, ylabelsize = labelfontsize)

                    hm = CairoMakie.heatmap!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], field_data[:, s, :], colormap = cmap, colorscale = colorbar_scale)

                    CairoMakie.contour!(ax, x_vec[xposmin:xposmax],  z_vec[zposmin:zposmax], data.fields.phase[xposmin:xposmax, s, zposmin:zposmax]; levels = levels, color = :white, linewidth = 0.5,)
                    CairoMakie.text!(ax.scene,textpos[1], textpos[2]; text = title, align = (:left, :bottom), fontsize = titlefontsize, color = textcolor)
                    cb = CairoMakie.Colorbar(fig[1, 2], hm;label = cbar_label,ticks =  ticks_pos, width = cb_width, labelsize = colorbar_tickfontsize,ticklabelsize = colorbar_tickfontsize)
                    cb.height = ax.scene.viewport[].widths[2]
                    colgap!(fig.layout, 10)
                ##########################################
                # make directories and save figs
                output_field = joinpath(output_folder,p_fields[i])
                if !isdir(output_field)
                    mkdir(output_field)
                end

                fig_name = "Fig"*string(parse(Int,split(timestep,"_")[2]))*"_"*string(NumValue(data.y[1,s,1]))*".png"
                output_name=joinpath(output_field,fig_name)
                CairoMakie.save(output_name,fig)
            end
        end
    end
end
 
########################################################################################################################################################
# save field properties of specific phases for each timestep

"""

        find_field_properties_grid(data::CartData,subsubfolder_path::String,timestep::String,phases::Vector{Int64},FileName_pvtr::String,fields::Vector{String},output_path::String,out_folder::String)


Parameters
====
- `data` - CartData structure of the model for one timestep,
- `subsubfolder_path` - path to folder where timesteps are stored
- `timestep` - timestep name to saved,
- `phases` - for which phase should be the information saved,
- `FileName_pvtr` - name of the pvtr file within the timesteps
- `fields` - Fields which should be saved,
- `output_path` - path to where data should stored,
- `out_folder` - name of folder where data gets stored in

"""


function find_field_properties_grid(data::CartData,subsubfolder_path::String,timestep::String,phases::Vector{Int64},FileName_pvtr::String,fields::Vector{String},output_path::String,out_folder::String)

    # get correct path
    processing_folder = joinpath(subsubfolder_path,timestep)
    location = replace(processing_folder,"\\" => "/")*"/"
    
    # get phase indices in matrix
    indices = get_phase(location,FileName_pvtr,phases,false)

    # create matrix from position of indices
    matrix = get_phase_bool(location,FileName_pvtr,indices)

    # write information to dictonary
    field_info = Dict()  # Dictionary to store material properties
    field_info["indices"]= indices
    field_info["matrix"]= matrix

    # write field information to dictonary for each direction
    for i in eachindex(fields)

        field_n = Symbol(fields[i])

        if  fields[i] == "dev_stress" || fields[i] == "vel_gr_tensor"

            data_direction = ["_xx", "_xy", "_xz", "_yx","_yy","_yz", "_zx", "_zy", "_zz"]

            for j in eachindex(data_direction)

                field_name_dir = string(fields[i]) * string(data_direction[j]) # get field direction name
                field_data = getfield(data.fields, field_n)[j] # get field data direction
                field_values = field_data[matrix .== 1] # get field of phase
                field_info[string(field_name_dir)]= field_values # write to dictonary
            end

       
        elseif p_fields[i] == "tot_displ" || p_fields[i] =="velocity" || p_fields[i] =="moment_res" || p_fields[i] == "EHmax" || p_fields[i] == "SHmax"

            data_direction = ["_x", "_y", "_z"]
            
            for j in eachindex(data_direction)

                field_name_dir = string(fields[i]) * string(data_direction[j])# get field direction name
                field_data = getfield(data.fields, field_n)[j] # get field data direction
                field_values = field_data[matrix .== 1] # get field of phase
                field_info[string(field_name_dir)]= field_values # write to dictonary

            end

        else

            field_data = getfield(data.fields, field_n)# get field data direction
            field_values = field_data[matrix .== 1]  # get field of phase
            field_info[string(fields[i])]= field_values # write to dictonary
        end

    end

    # create output directory if not existing
    output_path = joinpath(output_path,out_folder)
    if !isdir(output_path)
        mkdir(output_path)
    end

    # create output name and path and save it
    file_name = out_folder*string(parse(Int,split(split(location, "Timestep_")[2],"_")[1]))*".jld2"
    output_name=joinpath(output_path,file_name)

    # jld2 files
    JLD2.jldsave(output_name; field_info)

    return
end 

#####################################################################################################################################
# saves genreall information about the grid and timesteps


"""

        find_general_grid_prop(data::CartData,time_file::Vector{String},output_path::String,out_folder::String,material_block::Dict)


Parameters
====
- `data` - CartData structure of the model for one timestep,
- `time_file` - timestep where general information gets extracted
- `output_path` - path to where data should stored,
- `out_folder` - name of folder where data gets stored in, 
- `material_block` - Dict where material properties for each phase is stored ,

"""


function find_general_grid_prop(data::CartData,time_file::Vector{String},output_path::String,out_folder::String,material_block::Dict)

    field_info = Dict()  # dictionary to store material properties
    field_info["x"]= data.x.val # store x coorindates
    field_info["y"]= data.y.val # store y coorindates
    field_info["z"]= data.z.val # store z coorindates
    field_info["phase_init"] = data.fields.phase # store Initial phase setup
    field_info["timesteps"] = time_file # save timesteps
    field_info["material_block"] = material_block # Save material properties of all phases

    # create output path
    output_path = joinpath(output_path,out_folder)
    if !isdir(output_path)
        mkdir(output_path)
    end

    file_name = out_folder *".jld2"
    output_name=joinpath(output_path,file_name)

    # jld2 files
    JLD2.jldsave(output_name; field_info)

    return
end


#####################################################################################################################################
# saves surface development for each timestep

"""

        find_surf_evolution(subsubfolder_path::String,timestep::String,FileName_pvts::String,surface_level::Vector{Any},output_path::String,out_folder::String)


Parameters
====
- `subsubfolder_path` - Path to folder where timesteps are stored
- `timestep` - Timestep name to saved,
- `FileName_pvts` - name of the pvtr file within the timesteps
- `surface_level` - surface level in coordinates,
- `output_path` - path to where data should stored,
- `out_folder` - name of folder where data gets stored in

"""


function find_surf_evolution(subsubfolder_path::String,timestep::String,FileName_pvts::String,surface_level::Vector{Any},output_path::String,out_folder::String)

    # read surface information
    surf = get_surf_timestep(subsubfolder_path,timestep,FileName_pvts,surface_level)

    surf_info = Dict() # Dictionary to store material properties
    surf_info["velocity"]= surf.fields.velocity # Store surface velocity 
    surf_info["topography"]= surf.fields.topography # Store surface topography

    # create output directory
    output_path = joinpath(output_path,out_folder)
    if !isdir(output_path)
        mkdir(output_path)
    end

    file_name = out_folder*string(parse(Int,split(split(timestep, "Timestep_")[2],"_")[1]))*".jld2"
    output_name=joinpath(output_path,file_name)

    # jld2 files
    JLD2.jldsave(output_name; surf_info)

    return
end

######################################################################################################################################################
# saves tracer development for each timestep


"""

        find_tracer_info(subsubfolder_path::String,timestep::String,FileName_pvtu::String,surface_level::Vector{Any},phase_numb::Vector{Int},output_path::String,out_folder::String)


Parameters
====
- `subsubfolder_path` - Path to folder where timesteps are stored
- `timestep` - Timestep name to saved,
- `FileName_pvtu` - name of the pvtr file within the timesteps
- `surface_level` - surface level in coordinates,
- `phase_numb` - extract tracers only from specified phases,
- `output_path` - path to where data should stored,
- `out_folder` - name of folder where data gets stored in

"""



function find_tracer_info(subsubfolder_path::String,timestep::String,FileName_pvtu::String,surface_level::Vector{Any},phase_numb::Vector{Int},output_path::String,out_folder::String)

    tracer = get_tracer_timestep(subsubfolder_path,timestep,FileName_pvtu) # read all tracer
    tracer_in_slab = findall(x -> x in phase_numb, tracer.fields.Phase) # read all tracer in phase
    tracer_indices = tracer_in_slab .- 1 # correct position numbers to tracer numbers (start zero and one)

    tracer_info = Dict()  # Dictionary to store material properties
    tracer_info["x"]= try tracer.x.val[tracer_indices] catch; nothing end # save x position
    tracer_info["y"]= try tracer.y.val[tracer_indices] catch; nothing end # save y position
    tracer_info["z"]= try tracer.z.val[tracer_indices] .- surface_level catch; nothing end # save z position
    tracer_info["Phase"]= try tracer.fields.Phase[tracer_indices] catch; nothing end # save phase
    tracer_info["Temperature"]= try tracer.fields.Temperature[tracer_indices] catch; nothing end # save temperature
    tracer_info["Pressure"]= try tracer.fields.Pressure[tracer_indices] catch; nothing end # save pressure
    tracer_info["ID"]= try tracer.fields.ID[tracer_indices] catch; nothing end # save ID of tracer
    tracer_info["MF"]= try tracer.fields.Mf[tracer_indices] catch; nothing end # save melt fraction of tracer
    tracer_info["MF_Grid"]= try tracer.fields.Mf_Grid[tracer_indices] catch; nothing end # save grid of melt fraction of tracer
    tracer_info["Active"]= try tracer.fields.Active[tracer_indices] catch; nothing end  # save state 


    # create output directory
    output_folder = joinpath(output_path,out_folder)
    if !isdir(output_folder)
        mkdir(output_folder)
    end


    filename = out_folder * string(split_at__to_type([timestep],2,"Int64")[1])*".jld2"
    output_name=joinpath(output_folder,filename)

    #jld2 files
    JLD2.jldsave(output_name; tracer_info)

    
end


"""

Examples
========

 ```julia


julia> FileName = "output"
julia> FileName_pvtr=FileName*".pvtr"   # Name of pvtr file change if pvtr has a different name
julia> FileName_pvtu=FileName*".pvtu"   # Name of pvtr file change if pvtu has a different name
julia> FileName_pvts=FileName*".pvts"   # Name of pvtr file change if pvts has a different name

julia> folder = "test"  # folder where timesteps are saved in
julia> current_folder = "test" # location of where FileName folder is stored
julia> model_path = "test/input_files/timestep/"
julia> dat_path = "test/input_files/Passive_tracer_ex2D.dat"
julia> output_dir =joinpath(current_folder, "output")

julia> if !isdir(output_dir)
julia>    mkdir(output_dir)
julia> end

julia> p_fields = ["phase","temperature"] # field which you saved
julia> y_slice = [15]           # which y slice should be looked at, best is in the middle

julia> dxdz = [0.3, 0.8, 0.2, 0.7] # Maximum and minimum x and z values in Coordinates --> Float numbers
julia> phase_to_save = [2,3]       # phases to save the location an properties
julia> numb_ticks = 6              # phases to save the location an properties
julia> Savefilefolder = "fields"   # folder where field information are stored 
julia> Savegenfolder = "general"   # folder where general information are stored 
julia> Savetracerfolder = "tracer" # folder where tracer information are stored 
julia> Savesurffolder = "surf"     # folder where surface information are stored 

################################################################################################################################################################################
################################################################################################################################################################################

julia> surface_level = search_for_model_constrains(dat_path, "surf_level")

julia> material_block = Dict{String, Dict{String, Dict{String, String}}}()  # Dictionary to store material properties
julia> material_block = search_for_phase_properties(dat_path, folder, "<MaterialStart>", "<MaterialEnd>")
julia> number_phases = length(material_block[folder])  # extract total number of phases

julia> time_file = filter(f -> startswith(f, "Time"), readdir(model_path)) # Extract time information --> timestep and time

################################################################################################################################################################################
################################################################################################################################################################################

julia> for timestep in time_file       # go through all timesteps of a model

            data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level) # read model information

            post_plot(data,p_fields,timestep,number_phases,y_slice,output_dir, surface_level, dxdz,textpos, numb_ticks;phase_contour=true) # plot fields to see overall development

            if timestep == "Timestep_00000000_0.00000000e+00"   
                find_general_grid_prop(data,timestep,output_dir,Savegenfolder,material_block)    # get general information from first time step
            end

            find_field_properties_grid(data,model_path,timestep,phase_to_save,FileName_pvtr,p_fields,output_dir,Savefilefolder)    #save field properties

            #find_surf_evolution(model_path,timestep,FileName_pvts,surface_level,output_dir,Savesurffolder)   # save surface development

            #find_tracer_info(model_path,timestep,FileName_pvtu,surface_level,phase_to_save,output_dir,Savetracerfolder)  # save tracer development

julia> end
 ```

"""

#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################################################################################################

"""

    field_info = load_field_info(timestep::String,output_path::String,output_folder::String)

Parameters
==============

- `timestep` - Timestep Integer number to loaded,
- `output_path` - Path to folder where files are stored in,
- `output_folder` - name of files

"""

# Extract file information again --> FROM saved files

# load saved field data for one specific time step
function load_field_info(timestep::String,output_path::String,output_folder::String)

    # load saved field informations 
    output_slab = joinpath(output_path,output_folder)
    output_det = sort(readdir(output_slab), lt=natural) # sort files in output folder according to thte time steps

    file_name = filter(s -> occursin(timestep,s), output_det)[1] # read specific times step file

    output_file = joinpath(output_slab,file_name) # path to correct output file which should be read

    # preallocation
    phase_info = []

    jldopen(joinpath(output_file), "r") do f
        push!(phase_info,f[keys(f)[1]])
        end

    return phase_info[1]
end

############################################################################################################################
# load saved field data for all timesteps


"""


    field_info = load_field_info(output_path::String,output_folder::String)

Parameters
==============

- `output_path` - Path to folder where files are stored in,
- `output_folder` - Mame of files, loads all files in the folder with this name



Example
===============

julia> Savegenfolder = "general"   # folder where general information are stored 

julia> current_folder = "test" # location of where FileName folder is stored
julia> output_dir =joinpath(current_folder, "output")

julia> gen_info = load_field_info(output_dir,Savegenfolder)[1]

Dict{Any, Any} with 6 entries:
  "x"              => [0.0 0.0 … 0.0 0.0; 0.03125 0.03125 … 0.03125 0.03125; … ; 0.96875 0.96875 … 0.96875 0.96875; 1.0 1.0 … 1.0 1.0;;; 0.0 0.0 … 0.0 0.0; 0.03125 0.03125 … 0.03125 0.03125; … ; 0.96875 0.96875 ……
  "phase_init"     => Float32[0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.0 0.0 … 0.0 0.0;…
  "timesteps"      => "Timestep_00000000_0.00000000e+00"
  "z"              => [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.03125 0.03125 … 0.03125 0.03125; 0.03125 0.03125 … 0.03125 0.03125; … ; 0.03125 0.03125 … 0.03125 0.03125…
  "material_block" => Dict("test"=>Dict("Phase5"=>Dict("Cp"=>"1.2e3 # heat capacity", "ID"=>"5", "fr"=>"30    # friction angle [deg]", "k"=>"2.5", "alpha"=>"1e-5", "disl_prof"=>"Dry_Olivine_disl_creep-Hirth_Kohls…
  "y"              => [0.0 0.03125 … 0.96875 1.0; 0.0 0.03125 … 0.96875 1.0; … ; 0.0 0.03125 … 0.96875 1.0; 0.0 0.03125 … 0.96875 1.0;;; 0.0 0.03125 … 0.96875 1.0; 0.0 0.03125 … 0.96875 1.0; … ; 0.0 0.03125 … 0.9…


"""


function load_field_info(output_path::String,output_folder::String)

    output_slab = joinpath(output_path,output_folder) # get folder location
    phase_info = [] # preallocation

    if isdir(output_slab)
        output_det = sort(readdir(output_slab), lt=natural) # get all in folder

        # read each file in the folder and write it to the vector
        for t in eachindex(output_det) 
            output_file = joinpath(output_slab,output_det[t])
    
            jldopen(joinpath(output_file), "r") do f
                push!(phase_info,f[keys(f)[1]])
                end
        end

    elseif isfile(output_slab)
        # read each file in the folder and write it to the vector

        jldopen(joinpath(output_slab), "r") do f
            push!(phase_info,f[keys(f)[1]])
            end

    end

    return phase_info
end





####################################################################################################################

"""
  data =   track_point_over_time(Point_coord::CartesianIndex,folder::String,Savegenfolder::String,Savefieldfolder::String, name::String, output_path::String)


Parameters
====
- `Point_coord` - tracked point in CartesianIndex
- `folder` -  name of model
- `Savegenfolder` -  name of folder where general informations are stored
- `Savefieldfolder` -  name of folder where field informations are stored
- `surface_level` - get surface level for depth correction
- `name` -  file name to save information
- `output_path` - location where the output of the model is stored


Examples  
========

```julia

julia> Point_coord = CartesianIndex(10, 10, 10)
julia> folder = "PEV"
julia> Savegenfolder = "general"
julia> Savefieldfolder = "field"
julia> name = "track_output"
julia> output_path = "test/input_files/output/"

julia> track_point =  track_point_over_time(Point_coord::CartesianIndex,folder::String,Savegenfolder::String,Savefieldfolder::String, name::String, output_path::String)

Dict{String, Dict{String, Float64}} with 2 entries:
  "Timestep_00000010_1.49079279e+01" => Dict("phase"=>1.53445, "temperature"=>0.0)
  "Timestep_00000000_0.00000000e+00" => Dict("phase"=>2.0, "temperature"=>0.0)
   ...
 ```
"""



# track one specific point over time for multiple fields

function track_point_over_time(Point_coord::CartesianIndex,folder::String,Savegenfolder::String,Savefieldfolder::String, name::String, output_path::String)

    # read saved dictonaries for phases, generall info, timesteps, and indices and saved fields
    phase_info = load_field_info(output_path,Savefieldfolder)
    gen_info = load_field_info(output_path,Savegenfolder)[1]
    time_files = gen_info["timesteps"]
    fields = collect(keys(phase_info[1]))
    fields = filter(x -> x != "indices", fields)

    #get index of vector which contains only slab information
    ind = reduce(vcat, Iterators.flatten(phase_info[1]["indices"])) 
    index = findfirst(x -> x == Point_coord, ind)[1]

    # create empty dictory to save the time dependent locations of points 
    track_point = Dict{String,Dict{String,Float64}}()
    save_dict = Dict{String,Dict{String,Dict{String,Float64}}}()

    for t in eachindex(time_files)

        track_point[time_files[t]] = Dict{String,Float64}() # add timestep information to Dictonary

        for field in fields


            field_value = phase_info[t][field][index] # extract field information for indices at specific time step
            track_point[time_files[t]][field] = field_value# add field and index information to Dictonary

        end
    end

    # add model name to dictonary to save all models in one
    save_dict[folder] = Dict{String,Dict{String,Float64}}()
    save_dict[folder] = track_point

    # create output information
    file_name = string(name)*".jld2"
    output_name=joinpath(output_dir,file_name)

    jldopen(output_name, "w") do f
        f["save_dict"] = save_dict
        end
    
    return save_dict

end


