using GeophysicalModelGenerator, LaMEM, Serialization

# load data file
dat_path = ("../test/test_files/Subduction_VEP.dat")
model_path = ("../test/test_files/timestep/")
model_name = "VEP"
timestep = "Timestep_00000000_0.00000000e+00"
FileName_pvtr = "output.pvtr"
p_fields = ["phase", "temperature"]
output_dir = model_path

# extract data information from data file
surface_level = search_for_model_constrains(dat_path, "surf_level")

# read output file
material_block = Dict{String, Dict{String, Dict{String, String}}}()  # Dictionary to store material properties
material_block = search_for_phase_properties(dat_path, model_name, "<MaterialStart>", "<MaterialEnd>")

# get the time as a float number
time = split_at__to_type([timestep],3,"Float")

# extract data information
data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level,false)

# get information about where the phase is located
processing_folder = joinpath(model_path,timestep)
path = replace(processing_folder,"\\" => "/")*"/"
indices = get_phase(path,FileName_pvtr,[2],false)
matrix = get_phase_bool(path,FileName_pvtr,indices)

# track one point over time
name = "track"*string(indices[1])
track_point = track_point_over_time(indices[1],p_fields,model_name,model_path,surface_level,name,output_dir,false)

tracked_point = deserialize_file(output_dir,name)





