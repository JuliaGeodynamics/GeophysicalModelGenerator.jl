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

# Closed tetrahedron with known inside/outside points
@testset "isinside_closed_STL tetrahedron" begin
    using GeometryBasics: Point3f, TriangleFace, Mesh
    pts = Point3f[(0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4)]
    faces = TriangleFace.([(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)])
    tet = Mesh(pts, faces)
    @test isinside_closed_STL(tet, [1.0, 1.0, 1.0])
    @test !isinside_closed_STL(tet, [3.0, 3.0, 3.0])
    @test !isinside_closed_STL(tet, [-1.0, 1.0, 1.0])
    @test !isinside_closed_STL(tet, [10.0, 0.5, 0.5])
    @test isinside_closed_STL(tet, [0.5, 0.5, 2.5])
    @test isinside_closed_STL(tet, [3.0, 0.5, 0.4])
end
