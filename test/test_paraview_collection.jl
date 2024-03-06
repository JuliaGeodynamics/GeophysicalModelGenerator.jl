using Test
using GeophysicalModelGenerator, WriteVTK
@testset "Paraview collection" begin

x, y, z = 0:10, 1:6, 2:0.1:3
times = range(0, 1; step = 1)

#generate `*.vti` files
for (n, time) ∈ enumerate(times)
    vtk_grid("./test_files/test_vti_$n", x, y, z) do vtk
        vtk["Pressure"] = rand(length(x), length(y), length(z))
    end
end

# Generate a 3D grid
Lon,Lat,Depth   =   LonLatDepthGrid(10:20,30:40,(-300:25:0)km);
Data            =   Depth*2; # some data
Data_set        =   GeoData(Lon,Lat,Depth,(Depthdata=Data,LonData=Lon))
Write_Paraview(Data_set, "./test_files/test_depth3D")

make_paraview_collection(;dir = "./test_files", pvd_name="test", file_extension=".vti")
@test isfile("test.pvd")
@test filesize("test.pvd") == 317

make_paraview_collection(;dir = "./test_files", file_extension=".vti")
@test isfile("full_simulation.pvd")
@test filesize("full_simulation.pvd") == 317

make_paraview_collection(;dir = "./test_files")
@test isfile("full_simulation.pvd")
@test filesize("full_simulation.pvd") == 251


files = ["test_files/test_vti_1.vti", "test_files/test_vti_2.vti"]
time  = ["1.0", "2.0"]
make_paraview_collection("test2", files, time)
@test isfile("test2.pvd")
@test filesize("test2.pvd") == 317

make_paraview_collection(; pvd_name="test3", files=files, time=time)
@test isfile("test3.pvd")
@test filesize("test3.pvd") == 317

rm("test.pvd")
rm("full_simulation.pvd")
rm("test_files/test_depth3D.vts")
rm("test_files/test_vti_1.vti")
rm("test_files/test_vti_2.vti")
rm("test2.pvd")
rm("test3.pvd")

end
