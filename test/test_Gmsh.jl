using Test
using GeophysicalModelGenerator, GridapGmsh


# Read a Gmsh file
fname="test_files/subduction_ptatin.msh"
fe_data, tag_names = import_Gmsh(fname)
@test sum(fe_data.cellfields.regions) == 830

# Swap y and z dimensions (as pTatin uses a different definition)
data_fe = swap_yz_dims(fe_data)

# Define a CartData set with the same dimensions as the Gmsh file
bbox = extrema(data_fe);
nx,ny,nz = 100,50,80
data_cart = CartData( xyz_grid(range(bbox[1]...,length=nx),range(bbox[2]...,length=ny),range(bbox[3]...,length=nz) ))

data_cart1 = project_FEData_CartData(data_cart, data_fe)
@test extrema(data_cart1.fields.regions) == (2, 11)
