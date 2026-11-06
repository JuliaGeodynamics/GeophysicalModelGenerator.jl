using Test
using GeophysicalModelGenerator, GridapGmsh


# Read a Gmsh file
fname = "test_files/subduction_ptatin.msh"
fe_data, tag_names = import_Gmsh(fname)
@test sum(fe_data.cellfields.regions) == 830

# Swap y and z dimensions (as pTatin uses a different definition)
data_fe = swap_yz_dims(fe_data)

# the swap keeps the orientation of every tetrahedron
triple(a, b, c) = a[1] * (b[2] * c[3] - b[3] * c[2]) - a[2] * (b[1] * c[3] - b[3] * c[1]) + a[3] * (b[1] * c[2] - b[2] * c[1])
tet_volume(v, c) = triple(v[:, c[2]] - v[:, c[1]], v[:, c[3]] - v[:, c[1]], v[:, c[4]] - v[:, c[1]])
@test sign.(tet_volume.(Ref(fe_data.vertices), eachcol(fe_data.connectivity))) ==
    sign.(tet_volume.(Ref(data_fe.vertices), eachcol(data_fe.connectivity)))

# Define a CartData set with the same dimensions as the Gmsh file
bbox = extrema(data_fe);
nx, ny, nz = 100, 50, 80
data_cart = CartData(xyz_grid(range(bbox[1]..., length = nx), range(bbox[2]..., length = ny), range(bbox[3]..., length = nz)))

data_cart1 = project_FEData_CartData(data_cart, data_fe)
@test extrema(data_cart1.fields.regions) == (2, 11)

# 2D meshes are rejected
@test_throws "supports 3D meshes only" import_Gmsh("test_files/square_2D.msh")
