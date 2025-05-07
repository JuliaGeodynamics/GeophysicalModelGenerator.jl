using GeophysicalModelGenerator, LaMEM, Serialization, Test


# load data file
dat_path = ("./test_files/Subduction_VEP.dat")
model_path = ("./test_files/timestep/")
model_name = "VEP"
timestep = "Timestep_00000000_0.00000000e+00"
FileName_pvtr = "output.pvtr"
p_fields = ["phase", "temperature"]
output_dir = model_path

# test extraction from ascii files
# extract data information from data file
surface_level = search_for_model_constrains(dat_path, "surf_level")
gravity = search_for_model_constrains(dat_path, "gravity")

@test surface_level == [0.0]
@test gravity == [0.0, 0.0, 10.0]

# read output file
material_block = Dict{String, Dict{String, Dict{String, String}}}()  # Dictionary to store material properties
material_block = search_for_phase_properties(dat_path, model_name, "<MaterialStart>", "<MaterialEnd>")

@test material_block[model_name]["Phase5"]["fr"] == "30    # friction angle [deg]"
@test material_block[model_name]["Phase0"]["rho"] == "100"


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
@test sum(matrix) == 343
@test length(indices) == 343
@test matrix[indices[1]] == 1
@test matrix[CartesianIndex(1,1,1)] == 0


# extract data information
pvtr_path = joinpath(model_path,timestep)*"/"
data_pvtr = read_LaMEM_PVTR_file(pvtr_path,FileName_pvtr;fields = p_fields)
data = get_data_timestep(model_path,timestep,FileName_pvtr,p_fields,surface_level,false)

@test (getindex(data.z[1,1,1])) == (getindex(data_pvtr.z[1,1,1]))


# track one point over time
name = "track"*string(indices[1])
track_point = track_point_over_time(indices[1],p_fields,model_name,model_path,surface_level,name,output_dir,false)

@test track_point["Timestep_00000010_1.49079279e+01"]["phase"] == 1.534447431564331
@test track_point["Timestep_00000000_0.00000000e+00"]["phase"] == 2.0


tracked_point = deserialize_file(output_dir,name)
@test first(keys(track_point)) == "Timestep_00000010_1.49079279e+01"