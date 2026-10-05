using Test

# Collapse whitespace so the tests do not depend on column alignment
showstr(x) = replace(sprint(show, x), r"\s+" => " ")

@testset "GeoData" begin
    Lon, Lat, Depth = lonlatdepth_grid(10:20, 30:40, -10.0:-1)
    Depth[1] = NaN
    s = showstr(GeoData(Lon, Lat, Depth, (Vp = Depth,)))
    @test occursin("GeoData", s)
    @test occursin("has NaN's", s)
    @test occursin("fields : (:Vp,)", s)
    @test !occursin("attributes", s)

    Depth[1] = -10.0
    s = showstr(GeoData(Lon, Lat, Depth, (Vp = Depth,), Dict("note" => "test")))
    @test !occursin("has NaN's", s)
    @test occursin("attributes", s)
end

@testset "CartData / ParaviewData / Q1Data" begin
    X, Y, Z = xyz_grid(1:10, 1:10, -10:-1)
    Z[1] = NaN
    for T in (CartData, ParaviewData)
        s = showstr(T(X, Y, Z, (Z = Z,)))
        @test occursin(string(nameof(T)), s)
        @test occursin("has NaN's", s)
    end
    s = showstr(Q1Data(X, Y, Z, (Z = Z,), NamedTuple()))
    @test occursin("Q1Data", s)
    @test occursin("has NaN's", s)
    X, Y, Z = xyz_grid(1:10, 1:10, -10:-1)
    s = showstr(CartData(X, Y, Z, (Z = Z,), Dict("note" => "test")))
    @test !occursin("has NaN's", s)
    @test occursin("attributes", s)
    @test !occursin("has NaN's", showstr(ParaviewData(X, Y, Z, (Z = Z,))))
    @test !occursin("has NaN's", showstr(Q1Data(X, Y, Z, (Z = Z,), NamedTuple())))
end

@testset "UTMData" begin
    EW, NS, Depth = xyz_grid(422123.0:100:433623.0, 4.514137e6:100:4.523637e6, -5400:250:600)
    Depth[1] = NaN
    s = showstr(UTMData(EW, NS, Depth, 33, true, (Depth = Depth,)))
    @test occursin("UTM zone : 33-33 North", s)
    @test occursin("has NaNs", s)

    Depth[1] = 0.0
    s = showstr(UTMData(EW, NS, Depth, 33, false, (Depth = Depth,), Dict("note" => "test")))
    @test occursin("UTM zone : 33-33 South", s)
    @test !occursin("has NaNs", s)
    @test occursin("attributes", s)
end

@testset "FEData" begin
    fe_data = convert2FEData(Q1Data(xyz_grid(1:3, 1:3, 1:3)))
    s = showstr(fe_data)
    @test occursin("FEData{3,8}", s)
    @test occursin("elements : 8", s)
    @test occursin("vertices : 27", s)
end

@testset "CartGrid" begin
    @test occursin("x ∈ [0.0, 10.0]", showstr(create_CartGrid(size = 10, x = (0.0, 10))))
    @test occursin("x ∈ [0.0, 10.0], z ∈ [2.0, 10.0]", showstr(create_CartGrid(size = (10, 30), x = (0.0, 10), z = (2.0, 10))))
    s = showstr(create_CartGrid(size = (10, 20, 30), x = (0.0, 10), y = (0.0, 10), z = (2.0, 10)))
    @test occursin("CartGrid{Float64, 3}", s)
    @test occursin("x ∈ [0.0, 10.0], y ∈ [0.0, 10.0], z ∈ [2.0, 10.0]", s)
end

@testset "LaMEM" begin
    s = showstr(read_LaMEM_inputfile("test_files/SaltModels.dat"))
    @test occursin("LaMEM Grid:", s)
    @test occursin("nel :", s)

    P = setup_model_domain([-35.0, 26.0], [-24.0, 34.0], [-6.4, 6.4], 128, 64, 32, 8)
    s = showstr(P)
    @test occursin("LaMEM Partitioning info:", s)
    @test occursin("nProcX : 4", s)
end

@testset "GMG_Dataset" begin
    @test showstr(GMG_Dataset("Moho", "Surface", "moho.jld2", true)) == "GMG Surface Dataset (active) : Moho @ moho"
    @test showstr(GMG_Dataset("Moho", "Surface", "moho", true)) == "GMG Surface Dataset (active) : Moho @ moho"
    @test occursin("(inactive)", showstr(GMG_Dataset("Moho", "Surface", "moho", false)))
end

@testset "Trench" begin
    s = showstr(Trench(direction = -1.0, WeakzoneThickness = 10, WeakzonePhase = 9))
    @test occursin("left to right", s)
    @test occursin("Weakzone phase : 9", s)
    s = showstr(Trench(direction = 1.0))
    @test occursin("right to left", s)
    @test !occursin("Weakzone phase", s)
end
