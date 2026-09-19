using GeophysicalModelGenerator, JLD2, Test, CairoMakie, Printf, Statistics

# load data file
dat_path = ("../test/test_files/Subduction2D_LaMEM.dat") # path to dat file
model_path = ("../test/test_files/")             # path to model timesteps
output_dir =("../test/test_files/output/")       # path to output_folder
model_name = "Subduction"                        # name of the model
timestep = "Timestep_00000000_0.00000000e+00"    # name of Timestep
FileName = "output"                              # name of model output files
FileName_pvtr=FileName*".pvtr"                   # name of pvtr file 
FileName_pvtu=FileName*"_passive_tracers.pvtu"   # name of pvtu file 
FileName_pvts=FileName*"_surf.pvts"              # name of pvts file 
p_fields = ["phase", "temperature"]              # field to save
Savefieldfolder = "fields"                       # folder to store field information 
Savegenfolder = "general"                        # folder to store general information
Savetracerfolder = "tracer"                      # folder to store tracer information
Savesurffolder = "surf"                          # folder to store surface information
y_slice = [1]                                    # slice in y-direction which should be looked at
dxdz = [-1000.0, 1000.0, -600.0, -50.0]          # pLot window, maximum and minimum x and z values in Coordinates --> Float numbers
textpos = [-500.0, -500.0]                       # position of the text on the field plot
numb_ticks = 6                                   # number of ticks on axis
phase_to_save = [2,3]                            # phases to save the fields
Point_coord = CartesianIndex(282, 1, 81)         # node coordinate to track over time
name = "track_point"                             # name to save the tracked point

# extract data information from data file
surface_level = search_for_model_constrains(dat_path, "surf_level")
gravity = search_for_model_constrains(dat_path, "gravity")
rhos = search_for_all_model_constrains(dat_path, "rho")

@test gravity == [0.0, 0.0, 9.81]
@test rhos == [3300.0, 3300.0, 3300.0, 2700.0, 3300.0, 50.0]
@test surface_level == Any[0.0]

# read output file
material_block = Dict{String, Dict{String, Dict{String, String}}}()  # dictionary to store material properties
material_block = search_for_phase_properties(dat_path, model_name, "<MaterialStart>", "<MaterialEnd>")

@test material_block[model_name]["Phase1"]["ID"] == "1     # Material phase ID"
@test material_block[model_name]["Phase0"]["rho"] =="3300.0     # Density [kg/m^3]"

# test time extraction
# get the time as a float number
time3 = split_at__to_type([timestep],3,"Float")
time2 = split_at__to_type([timestep],2,"Int")
@test time3 == [0.0]
@test time2 == [0]

# get information about where the phase is located
processing_folder = joinpath(model_path,timestep)
path = replace(processing_folder,"\\" => "/")*"/"
indices = get_phase(path,FileName_pvtr,[2],false)
matrix = get_phase_bool(path,FileName_pvtr,indices)
@test sum(matrix) == 9406
@test length(indices) == 9406
@test matrix[indices[1]] == 1
@test matrix[CartesianIndex(1,1,1)] == 0


# extract data information
pvt_path = joinpath(model_path,timestep)*"/"
data_pvtr = read_LaMEM_PVTR_file(pvt_path,FileName_pvtr;fields = p_fields)
surf_pvts = read_LaMEM_PVTS_file(pvt_path,FileName_pvts)
tracer_pvtu = read_LaMEM_PVTU_file(pvt_path,FileName_pvtu)
data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level)
surf = get_surf_timestep(model_path,timestep,FileName_pvts,surface_level)
tracer = get_tracer_timestep(model_path,timestep,FileName_pvtu)

@test (getindex(data.z[1,1,1])) == (getindex(data_pvtr.z[1,1,1]))
@test (getindex(surf.z[1,1,1])) == (getindex(surf_pvts.z[1,1,1]))
@test (getindex(tracer.z[1,1,1])) == (getindex(tracer_pvtu.z[1,1,1]))


number_phases = length(material_block[model_name])  # extract total number of phases
time_file = filter(f -> startswith(f, "Time"), readdir(model_path)) # extract time information --> timestep and time

post_plot(data,p_fields,time_file[1],number_phases,y_slice,output_dir, surface_level, dxdz,textpos, numb_ticks;phase_contour=true) # plot fields to see overall development

find_general_grid_prop(data,time_file,output_dir,Savegenfolder,material_block)    # get general information from first time step

find_field_properties_grid(data,model_path,time_file[1],phase_to_save,FileName_pvtr,p_fields,output_dir,Savefieldfolder)    #save field properties
find_surf_evolution(model_path,time_file[1],FileName_pvts,surface_level,output_dir,Savesurffolder)   # save surface development
find_tracer_info(model_path,time_file[1],FileName_pvtu,surface_level,phase_to_save,output_dir,Savetracerfolder)  # save tracer development

# load saved information
gen_info = load_field_info(output_dir,Savegenfolder)[1]
phase_info0 = load_field_info("0",output_dir,Savefieldfolder)
phase_info = load_field_info(output_dir,Savefieldfolder)
surf_info = load_field_info(output_dir,Savesurffolder)
tracer_info = load_field_info(output_dir,Savetracerfolder)[1]

@test collect(keys(gen_info)) == Any["x", "phase_init", "timesteps", "z", "material_block", "y"]
@test unique(Int.(round.(gen_info["phase_init"]))) == [0 , 1, 2, 3, 4, 5]

@test isfile(joinpath(output_dir,p_fields[1], "Fig0_-2.5.png"))

@test typeof(phase_info0) == Dict{Any,Any}
@test typeof(phase_info) == Vector{Any}
@test typeof(surf_info) == Vector{Any}
@test typeof(tracer_info) ==  Dict{Any, Any}


#track one point over time
track_point =  track_point_over_time(Point_coord,model_name,Savegenfolder,Savefieldfolder, name, output_dir)
tracker_info = load_field_info(output_dir,"track_point.jld2")[1]

@test first(keys(tracker_info[model_name])) == "Timestep_00000000_0.00000000e+00"
@test track_point[model_name]["Timestep_00000000_0.00000000e+00"]["phase"] == 2.0



