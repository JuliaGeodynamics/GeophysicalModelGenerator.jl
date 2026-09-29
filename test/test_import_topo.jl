# Tests `import_topo`, which downloads from the GMT data server without needing GMT.
#
# These hit the network, so they are written to be cheap: small regions, and the tiles are
# cached between runs.
using Test
using GeophysicalModelGenerator

@testset "import_topo" begin

    # the tile addressing, which is where the awkward cases live
    @testset "tile naming" begin
        G = GeophysicalModelGenerator
        # named after the south-west corner, zero padded
        @test G.tile_filename(28, -18, "earth_relief", "03s", "g") == "N28W018.SRTMGL3.jp2"
        @test G.tile_filename(0, -30, "earth_relief", "01m", "p") == "N00W030.earth_relief_01m_p.jp2"
        @test G.tile_filename(-30, 60, "earth_relief", "15s", "p") == "S30E060.earth_relief_15s_p.jp2"

        # pixel registration samples cell centres, gridline the corners, so one more of them
        @test G.tile_size_px("01m", 30, "p") == 1800
        @test G.tile_size_px("01m", 30, "g") == 1801
        @test G.tile_size_px("03s", 1, "g") == 1201

        # the sets that exist in only one registration
        @test G.native_registration("earth_relief", "03s", "p") == "g"
        @test G.native_registration("earth_relief", "15s", "g") == "p"
    end

    # `"@earth_relief_01m"`, as this function used to be called
    @testset "file= spelling" begin
        G = GeophysicalModelGenerator
        @test G.parse_topo_file("@earth_relief_01m") == ("earth_relief", "01m", "g")
        @test G.parse_topo_file("earth_relief_15s_p") == ("earth_relief", "15s", "p")
        @test G.parse_topo_file("@earth_gebco_30s_g") == ("earth_gebco", "30s", "g")
    end

    # La Palma: small, has both land and deep water, and a summit to check against
    limits = [-18.7, -17.1, 28.0, 29.2]

    @testset "downloads and georeferencing" begin
        Topo = import_topo(limits, res = "01m")
        @test Topo isa GeoData

        z = ustrip.(Topo.fields.Topography)          # km
        @test size(z, 3) == 1
        @test -10 < minimum(z) < 0                   # deep water off the island
        @test 2 < maximum(z) < 3                     # Roque de los Muchachos, 2.4 km

        # the region we asked for, and no more
        @test minimum(Topo.lon.val) ≥ limits[1] - 0.1
        @test maximum(Topo.lon.val) ≤ limits[2] + 0.1
        @test minimum(Topo.lat.val) ≥ limits[3] - 0.1
        @test maximum(Topo.lat.val) ≤ limits[4] + 0.1
    end

    @testset "summit is in the right place" begin
        # 3 arcseconds resolves the summit; it should land on the real one to within a
        # pixel, which is what catches a transposed or mis-assembled grid
        Topo = import_topo(limits, res = "03s")
        z = ustrip.(Topo.fields.Topography)
        i = argmax(z)
        @test isapprox(Topo.lon.val[i], -17.885, atol = 0.002)
        @test isapprox(Topo.lat.val[i], 28.754, atol = 0.002)
        @test isapprox(maximum(z) * 1000, 2414, atol = 2)      # metres
    end

    @testset "bathymetry" begin
        # 3 arcseconds is SRTM, which is land only -- the sea comes from the coarser set the
        # server names for it, so it must not be flat
        Topo = import_topo(limits, res = "03s")
        z = ustrip.(Topo.fields.Topography)
        @test minimum(z) < -3                                   # km, real depth
        @test count(<(-1), z) > 0.5 * length(z)                 # mostly ocean here
    end

    @testset "other ways of asking" begin
        a = import_topo(limits, res = "30s")
        b = import_topo(lon = [limits[1], limits[2]], lat = [limits[3], limits[4]], res = "30s")
        @test ustrip.(a.fields.Topography) == ustrip.(b.fields.Topography)

        c = import_topo(limits, file = "@earth_relief_30s")
        @test size(c.fields.Topography) == size(a.fields.Topography)
    end

    @testset "other datasets" begin
        gebco = import_topo(limits, dataset = "earth_gebco", res = "15s")
        @test minimum(ustrip.(gebco.fields.Topography)) < -3     # km
    end
end
