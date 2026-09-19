using GeophysicalModelGenerator, LaMEM, JLD2
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
track_name = "track_point"                             # name to save the tracked point


# extract data information from data file
surface_level = search_for_model_constrains(dat_path, "surf_level")

# read output file
material_block = Dict{String, Dict{String, Dict{String, String}}}()  # dictionary to store material properties
material_block = search_for_phase_properties(dat_path, model_name, "<MaterialStart>", "<MaterialEnd>")
number_phases = length(material_block[model_name])  # extract total number of phases

# get the time as a float number
time_file = filter(f -> startswith(f, "Time"), readdir(model_path)) # Extract time information --> timestep and time
time = split_at__to_type([timestep],3,"Float")

# extract data information from models
data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level)   # model fields
surf = get_surf_timestep(model_path,timestep,FileName_pvts,surface_level)            # surface development
tracer = get_tracer_timestep(model_path,timestep,FileName_pvtu)                      # tracer development

# get information about where specific phases are located
processing_folder = joinpath(model_path,timestep)
path = replace(processing_folder,"\\" => "/")*"/"
indices = get_phase(path,FileName_pvtr,[2],false)
matrix = get_phase_bool(path,FileName_pvtr,indices)

# save general grid properties
find_general_grid_prop(data,time_file,output_dir,Savegenfolder,material_block)  

# save data of fields, surface and tracer for each timestep
for timestep in time_file
    @show timestep
    data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level)    # load data from current timestep
    post_plot(data,p_fields,timestep,number_phases,y_slice,output_dir, surface_level, dxdz,textpos, numb_ticks;phase_contour=true) # plot fields to see overall development
    find_field_properties_grid(data,model_path,timestep,phase_to_save,FileName_pvtr,p_fields,output_dir,Savefieldfolder)    # save field properties
    find_surf_evolution(model_path,timestep,FileName_pvts,surface_level,output_dir,Savesurffolder)   # save surface development
    find_tracer_info(model_path,timestep,FileName_pvtu,surface_level,phase_to_save,output_dir,Savetracerfolder)  # save tracer development
end

# track a point over time
track_point =  track_point_over_time(Point_coord,model_name,Savegenfolder,Savefieldfolder, track_name, output_dir) # track one grid point over time

# load saved information
gen_info = load_field_info(output_dir,Savegenfolder)            # load general information
phase_info0 = load_field_info("0",output_dir,Savefieldfolder)   # load field info for one specific timestep
phase_info = load_field_info(output_dir,Savefieldfolder)        # load field info for all timestep
surf_info = load_field_info(output_dir,Savesurffolder)          # load surface information
tracer_info = load_field_info(output_dir,Savetracerfolder)      # load tracer information
tracker_info = load_field_info(output_dir,track_name*".jld2")   # load tracked point information




