using Test

# Q1 data set
q1_data = Q1Data(xyz_grid(1:10, 1:10, 1:8))
q1_data = addfield(q1_data, (T = ones(Float64, size(q1_data) .+ 1),))
q1_data = addfield(q1_data, (region = zeros(Int64, size(q1_data)),), cellfield = true)

# convert to FEData
fe_data = convert2FEData(q1_data)

out = write_pTatin_mesh(fe_data)
@test isnothing(out)
@test isfile("md.bin")
@test isfile("region_cell.bin")
@test isfile("T_vertex.bin")
@test isfile("Z_vertex.bin")

rm("md.bin")
rm("region_cell.bin")
rm("T_vertex.bin")
rm("Z_vertex.bin")
