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

# the field files go next to `out_file`, not into the working directory
dir = mktempdir()
write_pTatin_mesh(fe_data; out_file = joinpath(dir, "md.bin"))
@test isfile(joinpath(dir, "md.bin"))
@test isfile(joinpath(dir, "region_cell.bin"))
@test isfile(joinpath(dir, "T_vertex.bin"))
@test !isfile("T_vertex.bin")

# a mesh with as many cells as vertices per cell is valid
write_pTatin_mesh(Q1Data(xyz_grid(1:3, 1:3, 1:3)); out_file = joinpath(dir, "small.bin"))
@test isfile(joinpath(dir, "small.bin"))

# swapping y and z keeps the hexahedra positively oriented
triple(a, b, c) = a[1] * (b[2] * c[3] - b[3] * c[2]) - a[2] * (b[1] * c[3] - b[3] * c[1]) + a[3] * (b[1] * c[2] - b[2] * c[1])
corner_volume(v, c) = triple(v[:, c[2]] - v[:, c[1]], v[:, c[3]] - v[:, c[1]], v[:, c[5]] - v[:, c[1]])
fe_small = convert2FEData(Q1Data(xyz_grid(1:3, 1:3, 1:3)))
fe_swap = swap_yz_dims(fe_small)
@test all(c -> corner_volume(fe_small.vertices, c) > 0, eachcol(fe_small.connectivity))
@test all(c -> corner_volume(fe_swap.vertices, c) > 0, eachcol(fe_swap.connectivity))
@test fe_swap.vertices[2, :] == fe_small.vertices[3, :]

# ditto for a tetrahedron
tet = FEData([0.0 1 0 0; 0 0 1 0; 0 0 0 1], reshape([1, 2, 3, 4], 4, 1))
tet_volume(v, c) = triple(v[:, c[2]] - v[:, c[1]], v[:, c[3]] - v[:, c[1]], v[:, c[4]] - v[:, c[1]])
@test tet_volume(tet.vertices, tet.connectivity[:, 1]) > 0
@test tet_volume(swap_yz_dims(tet).vertices, swap_yz_dims(tet).connectivity[:, 1]) > 0
