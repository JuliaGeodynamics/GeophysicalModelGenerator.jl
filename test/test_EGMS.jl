using Test, Dates
using GeophysicalModelGenerator
import GeophysicalModelGenerator: laea_to_lonlat, lonlat_to_laea

@testset "EGMS" begin
    # worked example of EPSG guidance note 7-2 for ETRS89-LAEA (EPSG:3035)
    E, N = lonlat_to_laea(5.0, 50.0)
    @test E ≈ 3962799.45 atol = 0.01
    @test N ≈ 2999718.85 atol = 0.01
    @test all(isapprox.(laea_to_lonlat(E, N), (5.0, 50.0), atol = 1.0e-7))
    @test laea_to_lonlat(4321000.0, 3210000.0) == (10.0, 52.0)

    # Synthetic files in the layout of the EGMS products: points on a 100 m grid that move
    # with a velocity that increases eastwards, plus a seasonal signal
    dates = collect(Date(2019, 1, 5):Day(6):Date(2023, 12, 28))
    t = [Dates.value(d - dates[1]) / 365.25 for d in dates]
    vel_true(i) = -10.0 + i
    dir = mktempdir()

    function write_file(name, level; i_range = 1:20)
        open(joinpath(dir, name), "w") do io
            # L3 with the column names of the releases from 2020_2024 on, L2 with the older ones
            hdr = level == 3 ? "pid,easting,northing,height_ortho,rmse_ts,mean_velocity,mean_velocity_std,acceleration,acceleration_std,seasonality,seasonality_std,gnss_velocity_n,gnss_velocity_e,gnss_velocity_u" :
                "pid,mp_type,latitude,longitude,easting,northing,height,height_wgs84,line,pixel,rmse,temporal_coherence,amplitude_dispersion,incidence_angle,track_angle,los_east,los_north,los_up,mean_velocity,mean_velocity_std,acceleration,acceleration_std,seasonality,seasonality_std"
            println(io, hdr, ",", join(Dates.format.(dates, "yyyymmdd"), ","))
            for i in i_range, j in 1:10
                E, N = 4.4e6 + 100 * i, 2.6e6 + 100 * j
                local lon, lat = laea_to_lonlat(E, N)
                disp = vel_true(i) .* t .+ 2.0 .* sin.(2π .* t)
                pre = level == 3 ? "p$(i)_$j,$E,$N,150.0,0.5,$(vel_true(i)),0.1,0.0,0.1,2.0,0.1,2.4,0.5,-0.6" :
                    "p$(i)_$j,0,$lat,$lon,$E,$N,150.0,195.0,10,20,0.5,0.9,0.1,38.0,-10.0,-0.6,-0.1,0.78,$(vel_true(i)),0.1,0.0,0.1,2.0,0.1"
                println(io, pre, ",", join(round.(disp, digits = 1), ","))
            end
        end
        return joinpath(dir, name)
    end

    f3 = write_file("EGMS_L3_E44N26_100km_U_2019_2023_1.csv", 3)
    f2 = write_file("EGMS_L2b_044_0297_IW2_VV_2019_2023_1.csv", 2)

    d3 = import_EGMS(f3)
    d2 = import_EGMS(f2)
    @test length(d3) == length(d2) == 200
    @test d3.dates == dates
    @test size(d3.displacement) == (length(dates), 200)
    @test d3.lon ≈ d2.lon           # L3 from easting/northing, L2 from its longitude column
    @test d3.lat ≈ d2.lat
    @test all(d3.height .== 150.0)           # from `height_ortho`
    @test all(d2.height .== 150.0)           # from `height`
    @test haskey(d3.fields, :gnss_velocity_u)
    @test !haskey(d3.fields, :height_ortho)
    @test d3.fields.mean_velocity[1] == vel_true(1)
    @test !haskey(d3.fields, :pid)
    @test haskey(d2.fields, :temporal_coherence)
    @test occursin("EGMSData", sprint(show, d3))

    # velocities: full period, and a window of whole years
    v = EGMS_velocity(d3)
    @test v[1] ≈ vel_true(1) atol = 0.1
    @test v[end] ≈ vel_true(20) atol = 0.1
    v = EGMS_velocity(d3, from = Date(2021, 1, 1), to = Date(2022, 12, 31))
    @test v[end] ≈ vel_true(20) atol = 0.3
    @test_throws ErrorException EGMS_velocity(d3, from = Date(2030, 1, 1))

    @test EGMS_velocity(d3, seasonal = false)[end] ≈ vel_true(20) atol = 0.3
    @test_throws ErrorException EGMS_velocity(d3, from = Date(2021, 1, 1), to = Date(2021, 6, 1), seasonal = true)
    @test !isnan(EGMS_velocity(d3, from = Date(2021, 1, 1), to = Date(2021, 6, 1))[1])

    # missing values are skipped
    d3.displacement[5, 1] = NaN32
    @test EGMS_velocity(d3)[1] ≈ vel_true(1) atol = 0.1

    # limits, several files and a directory
    lims = [minimum(d3.lon), sum(extrema(d3.lon)) / 2, -90, 90]
    dl = import_EGMS(f3, limits = lims)
    @test 0 < length(dl) < 200
    @test all(lims[1] .<= dl.lon .<= lims[2])
    fb = write_file("EGMS_L3_E45N26_100km_U_2019_2023_1.csv", 3; i_range = 21:30)
    @test length(import_EGMS([f3, fb])) == 300
    @test_throws ErrorException import_EGMS([f3, f2])     # different columns
    @test_throws ErrorException import_EGMS(joinpath(dir, "nofile.csv"))
    dir3 = mktempdir()
    cp(f3, joinpath(dir3, basename(f3))); cp(fb, joinpath(dir3, basename(fb)))
    @test length(import_EGMS(dir3)) == 300

    # points
    P = EGMS_points(d3)
    @test P isa GeoData
    @test P.depth.val[1] ≈ 0.15
    @test P.fields.velocity[end] ≈ vel_true(20) atol = 0.1

    # on a grid: coarse cells average the points in them, and cells without points are NaN
    lon = range(minimum(d3.lon) - 0.01, maximum(d3.lon) + 0.01, length = 12)
    lat = range(minimum(d3.lat) - 0.01, maximum(d3.lat) + 0.01, length = 8)
    G = EGMS_on_grid(d3, lon, lat)
    @test size(G.fields.velocity) == (12, 8, 1)
    @test isnan(G.fields.velocity[1, 1, 1])
    @test G.fields.npoints[1, 1, 1] == 0
    inside = findall(!isnan, G.fields.velocity)
    @test !isempty(inside)
    @test all(vel_true(1) - 0.1 .<= G.fields.velocity[inside] .<= vel_true(20) + 0.1)
    @test sum(G.fields.npoints) > 0

    # smooth: fills more of the grid, varies less from node to node, and stays within the data
    Gs = EGMS_on_grid(d3, range(extrema(d3.lon)..., length = 60), range(extrema(d3.lat)..., length = 8), radius = 0.6, smooth = true)
    Gp = EGMS_on_grid(d3, range(extrema(d3.lon)..., length = 60), range(extrema(d3.lat)..., length = 8), radius = 0.6)
    vs, vp = Gs.fields.velocity[:, 4, 1], Gp.fields.velocity[:, 4, 1]
    @test !any(isnan, vs)
    @test minimum(EGMS_velocity(d3)) <= minimum(vs) && maximum(vs) <= maximum(EGMS_velocity(d3))
    @test vs[1] < vp[1] && vs[end] > vp[end]      # less smeared out at the edges than the plain average
    @test issorted(vs)

    # a step in the velocity (a fault) stays sharp where data is dense, however large the radius
    step = import_EGMS(f3)
    step.displacement .= Float32.([d > Date(2019, 1, 5) ? (x > sum(extrema(step.lon)) / 2 ? 10.0 : -10.0) * Dates.value(d - Date(2019, 1, 5)) / 365.25 : 0.0 for d in step.dates, x in step.lon])
    glon, glat = range(extrema(step.lon)..., length = 60), range(extrema(step.lat)..., length = 8)
    sharp = EGMS_on_grid(step, glon, glat, radius = 1.0, smooth = true).fields.velocity[:, 4, 1]
    blurred = EGMS_on_grid(step, glon, glat, radius = 1.0, smooth = true, minpoints = Inf).fields.velocity[:, 4, 1]
    @test count(v -> abs(v) < 9, sharp) < count(v -> abs(v) < 9, blurred)
    @test count(v -> abs(v) < 9, sharp) <= 6

    # vector output: (east, north, up); without an east component only the vertical one is filled
    @test G.fields.Velocity[3] isa Array && isequal(G.fields.Velocity[3], G.fields.velocity)
    @test all(G.fields.Velocity[1] .== 0) && all(G.fields.Velocity[2] .== 0)
    fe = write_file("EGMS_L3_E44N26_100km_E_2019_2023_1.csv", 3)
    Gv = EGMS_on_grid(d3, lon, lat, east = import_EGMS(fe))
    @test isequal(Gv.fields.Velocity[1], Gv.fields.velocity_east)
    @test isapprox(Gv.fields.velocity_east, Gv.fields.velocity, atol = 0.1, nans = true)      # the test file has the same values
    before = copy(Gv.fields.velocity)
    write_paraview(Gv, joinpath(dir, "EGMS_vectors"))
    @test isequal(Gv.fields.velocity, before)      # scalar fields untouched by the rotation on writing

    # on an existing surface: keeps its coordinates and fields
    Lon, Lat, Depth = lonlatdepth_grid(lon, lat, 0.0)
    Topo = GeoData(Lon, Lat, fill(0.15, size(Lon)), (Topography = fill(0.15, size(Lon)),))
    S = EGMS_on_grid(d3, Topo, radius = 0.5)
    @test haskey(S.fields, :Topography)
    @test S.depth.val == Topo.depth.val
    @test count(!isnan, S.fields.velocity) >= length(inside)
end
