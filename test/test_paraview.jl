# this tests creating paraview output from the data
using Test
using GeophysicalModelGenerator

@testset "Paraview" begin

    # Generate a 3D grid
    Lon, Lat, Depth = lonlatdepth_grid(10:20, 30:40, (-300:25:0)km)
    Data = Depth * 2  # some data
    Data_set = GeoData(Lon, Lat, Depth, (Depthdata = Data, LonData = Lon))
    @test write_paraview(Data_set, "test_depth3D") == nothing

    # Horizontal profile @ 10 km height
    Lon, Lat, Depth = lonlatdepth_grid(10:20, 30:40, 10km)
    Depth[2:4, 2:4, 1] .= 25km     # add some fake topography

    Data_set2 = GeoData(Lon, Lat, Depth, (Topography = Depth,))
    @test write_paraview(Data_set2, "test2") == nothing

    # Cross sections
    Lon, Lat, Depth = lonlatdepth_grid(10:20, 35, (-300:25:0)km)
    Data_set3 = GeoData(Lon, Lat, Depth, (DataSet = Depth,))
    @test write_paraview(Data_set3, "test3") == nothing


    Lon, Lat, Depth = lonlatdepth_grid(15, 30:40, (-300:25:0)km)
    Data_set4 = GeoData(Lon, Lat, Depth, (DataSet = Depth,))
    @test write_paraview(Data_set4, "test4") == nothing

    Lon, Lat, Depth = lonlatdepth_grid(15, 35, (-300:25:0)km)
    Data_set5 = GeoData(Lon, Lat, Depth, (DataSet = Depth,))
    @test write_paraview(Data_set5, "test5") == nothing


    # Test saving vectors
    Lon, Lat, Depth = lonlatdepth_grid(10:20, 30:40, 50km)
    Ve = zeros(size(Depth)) .+ 1.0
    Vn = zeros(size(Depth))
    Vz = zeros(size(Depth))
    Velocity = (copy(Ve), copy(Vn), copy(Vz))              # tuple with 3 values, which
    Data_set_vel = GeoData(Lon, Lat, Depth, (Velocity = Velocity, Veast = Velocity[1] * cm / yr, Vnorth = Velocity[2] * cm / yr, Vup = Velocity[3] * cm / yr))
    @test  write_paraview(Data_set_vel, "test_Vel") == nothing

    # Test saving colors
    red = zeros(size(Lon))
    green = zeros(size(Lon))
    blue = zeros(size(Lon))
    Data_set_color = GeoData(Lon, Lat, Depth, (Velocity = Velocity, colors = (red, green, blue), color2 = (red, green, blue)))
    @test write_paraview(Data_set_color, "test_Color") == nothing

    # Manually test the in-place conversion from spherical -> cartesian (done automatically when converting GeoData->ParaviewData  )
    Vel_Cart = (copy(Ve), copy(Vn), copy(Vz))
    velocity_spherical_to_cartesian!(Data_set_vel, Vel_Cart)
    @test Vel_Cart[2][15] ≈ 0.9743700647852352
    @test Vel_Cart[1][15] ≈ -0.224951054343865
    @test Vel_Cart[3][15] ≈ 0.0

    # Test saving unstructured point data (EQ, or GPS points)
    Data_set_VelPoints = GeoData(Lon[:], Lat[:], ustrip.(Depth[:]), (Velocity = (copy(Ve[:]), copy(Vn[:]), copy(Vz[:])), Veast = Ve[:] * mm / yr, Vnorth = Vn[:] * cm / yr, Vup = Vz[:] * cm / yr))
    @test write_paraview(Data_set_VelPoints, "test_Vel_points", PointsData = true) == nothing

end

@testset "Paraview output types, directories and movies" begin
    dir = mktempdir()
    X, Y, Z = xyz_grid(1:4, 1:5, -3:0)
    V = (copy(X), copy(Y), copy(Z))

    # movie with output in a directory
    movie = movie_paraview(name = joinpath(dir, "Movie"), Initialize = true)
    for itime in 1:2
        movie = write_paraview(CartData(X, Y, Z * itime, (Z = Z,)), "cart$itime"; pvd = movie, time = itime, directory = dir, verbose = false)
    end
    movie_paraview(pvd = movie, Finalize = true)
    @test isfile(joinpath(dir, "Movie.pvd"))
    @test isfile(joinpath(dir, "cart2.vts"))

    @test_throws "vector data fields have units" write_paraview(ParaviewData(X, Y, Z, (V = (X * cm / yr, Y * cm / yr, Z * cm / yr),)), joinpath(dir, "units"), verbose = false)

    EW, NS, Depth = xyz_grid(422123.0:100:422523.0, 4.514137e6:100:4.514537e6, -500:250:0)
    write_paraview(UTMData(EW, NS, Depth, 33, true, (Depth = Depth,)), "utm"; directory = dir, verbose = false)
    @test isfile(joinpath(dir, "utm.vts"))

    # Q1Data with vertex & cell fields, scalar & vector
    Xc, Yc, Zc = xyz_grid(1.5:3.5, 1.5:4.5, -2.5:-0.5)
    q1 = Q1Data(X, Y, Z, (Z = Z * km, V = V), (Zc = Zc * km, Vc = (Xc, Yc, Zc)))
    movie = movie_paraview(name = joinpath(dir, "Movie_q1"))
    movie = write_paraview(q1, "q1"; directory = dir, pvd = movie, time = 1.0, verbose = false)
    movie_paraview(pvd = movie, Finalize = true)
    @test isfile(joinpath(dir, "q1.vts"))
    @test isfile(joinpath(dir, "Movie_q1.pvd"))
    q1_units = Q1Data(X, Y, Z, (Z = Z,), (Vc = (Xc * km, Yc * km, Zc * km),))
    @test_throws "vector data fields have units" write_paraview(q1_units, joinpath(dir, "q1_units"), verbose = false)
    q1_units = Q1Data(X, Y, Z, (V = (X * km, Y * km, Z * km),), NamedTuple())
    @test_throws "vector data fields have units" write_paraview(q1_units, joinpath(dir, "q1_units"), verbose = false)

    # FEData: hexahedra (from Q1Data) and tetrahedra
    fe = convert2FEData(q1)
    movie = movie_paraview(name = joinpath(dir, "Movie_fe"))
    movie = write_paraview(fe, "fe_hex"; directory = dir, pvd = movie, time = 1.0, verbose = false)
    movie_paraview(pvd = movie, Finalize = true)
    @test isfile(joinpath(dir, "fe_hex.vtu"))

    vertices = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    tet = FEData(vertices, reshape([1, 2, 3, 4], 4, 1), (Z = vertices[3, :] * km, V = (vertices[1, :], vertices[2, :], vertices[3, :])), (C = [1.0] * km, Vc = ([1.0], [0.0], [0.0])))
    write_paraview(tet, joinpath(dir, "fe_tet"), verbose = false)
    @test isfile(joinpath(dir, "fe_tet.vtu"))
    tet_units = FEData(vertices, reshape([1, 2, 3, 4], 4, 1), (V = (vertices[1, :] * km, vertices[2, :] * km, vertices[3, :] * km),), NamedTuple())
    @test_throws "vector data fields have units" write_paraview(tet_units, joinpath(dir, "fe_units"), verbose = false)
    tet_units = FEData(vertices, reshape([1, 2, 3, 4], 4, 1), NamedTuple(), (Vc = ([1.0] * km, [0.0] * km, [0.0] * km),))
    @test_throws "vector data fields have units" write_paraview(tet_units, joinpath(dir, "fe_units"), verbose = false)

    wedge = FEData(rand(3, 6), reshape(1:6, 6, 1), NamedTuple(), NamedTuple())
    @test_throws "This element is not yet implemented" write_paraview(wedge, joinpath(dir, "wedge"), verbose = false)
end
