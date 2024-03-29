using Test
using GeophysicalModelGenerator, GridapGmsh


# Read a Gmsh file
fname="test_files/subduction_ptatin.msh"
fe_data, tag_names = import_Gmsh(fname)
@test sum(fe_data.cellfields.regions) == 830


