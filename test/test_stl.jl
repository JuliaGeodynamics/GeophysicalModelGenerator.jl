# test STL routines
using Test, GeophysicalModelGenerator

# Load cat MESH
mesh = load("./test_files/cat.stl")
X, Y, Z = xyz_grid(150:180, -15:2:15, 10:5:60) # Create mesh

# Test isinside_closed_STL routine for individual points (note: bit slow)
Phase = zeros(size(X));
for i in eachindex(X)

    inside = isinside_closed_STL(mesh, [X[i], Y[i], Z[i]])
    if inside
        Phase[i] = 1
    end
end

@test Phase[14, 6, 2] == 1.0

#Data_Cat = ParaviewData(X,Y,Z, (Phase=Phase,))
#write_paraview(Data_Cat,"Data_Cat")
