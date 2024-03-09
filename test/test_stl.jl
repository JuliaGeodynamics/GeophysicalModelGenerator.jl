# test STL routines
using Test, GeophysicalModelGenerator

@testset "STL" begin
    # Load cat MESH
    mesh    =   load("./test_files/cat.stl")
    X,Y,Z   =   XYZGrid(150:180, -15:2:15, 10:5:60) # Create mesh

    # Test IsInsideClosedSTL routine for individual points (note: bit slow)
    Phase = zeros(size(X));
    for i in eachindex(X)

        inside = IsInsideClosedSTL(mesh, [X[i], Y[i], Z[i]]) 
        if inside   
            Phase[i] = 1;
        end
    end

    @test Phase[14,6,2] == 1.0

    #Data_Cat = ParaviewData(X,Y,Z, (Phase=Phase,))
    #Write_Paraview(Data_Cat,"Data_Cat")
end