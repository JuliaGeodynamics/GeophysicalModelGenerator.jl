# Tests `import_topo`, which downloads from the GMT data server without needing GMT.
#
# Everything these ask for is downloaded up front by topo_prefetch.jl (keep its list in step
# with the calls below), so they are served from the tile cache and touch the network only
# when the prefetch could not.
using Test
using GeophysicalModelGenerator

# The prefetch has already tried the GMT data server; if it could not be reached there is no
# point in asking again here. On a push/PR run the download tests are skipped, so a hiccup of
# that server does not block unrelated work; the scheduled run insists, so a real breakage of
# the download is still caught.
topo_server_ok = get(ENV, "GMG_TOPO_PREFETCH_OK", "true") == "true" ||
    get(ENV, "GITHUB_EVENT_NAME", "") == "schedule"

@testset "import_topo" begin

    # the tile addressing, which is where the awkward cases live
    @testset "tile naming" begin
        G = GeophysicalModelGenerator
        # the SRTM sets live under earth_relief whichever dataset asks for them
        @test G.tile_url("N28W018.SRTMGL3.jp2", "earth_gebco", "03s", "g") == "earth_relief/earth_relief_03s_g/N28W018.SRTMGL3.jp2"
        @test G.tile_url("N00W030.earth_gebco_01m_p.jp2", "earth_gebco", "01m", "p") == "earth_gebco/earth_gebco_01m_p/N00W030.earth_gebco_01m_p.jp2"
        # limits outside the globe or without extent are rejected before any download
        @test_throws "must satisfy" import_topo([180, 185, 40, 41])
        @test_throws "must satisfy" import_topo([350, 355, 40, 41])
        @test_throws "must satisfy" import_topo([10, 10.5, 60, 60])
        @test_throws "antimeridian" import_topo([170, -170, -20, -10])
        # named after the south-west corner, zero padded
        @test G.tile_filename(28, -18, "earth_relief", "03s", "g") == "N28W018.SRTMGL3.jp2"
        @test G.tile_filename(0, -30, "earth_relief", "01m", "p") == "N00W030.earth_relief_01m_p.jp2"
        @test G.tile_filename(-30, 60, "earth_relief", "15s", "p") == "S30E060.earth_relief_15s_p.jp2"

        # pixel registration samples cell centres, gridline the corners, so one more of them
        @test G.tile_size_px("01m", 30, "p") == 1800
        @test G.tile_size_px("01m", 30, "g") == 1801
        @test G.tile_size_px("03s", 1, "g") == 1201

        # the sets that exist in only one registration (read from the server index, or from
        # the copy bundled with the package when the server is out of reach)
        @test G.native_registration("earth_relief", "03s", "p") == "g"
        @test G.native_registration("earth_relief", "15s", "g") == "p"
    end

    # the `"@earth_relief_01m"` spelling of the data set
    @testset "file= spelling" begin
        G = GeophysicalModelGenerator
        @test G.parse_topo_file("@earth_relief_01m") == ("earth_relief", "01m", "g")
        @test G.parse_topo_file("earth_relief_15s_p") == ("earth_relief", "15s", "p")
        @test G.parse_topo_file("@earth_gebco_30s_g") == ("earth_gebco", "30s", "g")
    end

    # the bundled index must know the tiled sets, or a missing server index would send a
    # tiled resolution down the single-file path
    @testset "bundled server index" begin
        G = GeophysicalModelGenerator
        @test isfile(G.BUNDLED_SERVER_INDEX)
        @test G.tile_size_deg("earth_relief", "01m", "g") == 30
        @test G.tile_size_deg("earth_relief", "15s", "p") == 10
        @test G.tile_size_deg("earth_relief", "03s", "g") == 1
        @test isnothing(G.tile_size_deg("earth_relief", "10m", "g"))
    end

    # what enters the cache, checked with local files rather than the server
    @testset "download validation" begin
        G = GeophysicalModelGenerator
        dir = mktempdir()
        html = joinpath(dir, "error.html")
        write(html, "<html><body>Service unavailable</body></html>")
        jp2 = joinpath(dir, "source.jp2")
        write(jp2, b"\0\0\0\x0cjP  \r\n\x87\n", zeros(UInt8, 20))
        fileurl(p) = "file://" * (Sys.iswindows() ? "/" * replace(p, '\\' => '/') : p)

        # an error page is not cached, and leaves no temporary file behind
        dest = joinpath(dir, "N00E000.earth_relief_01m_g.jp2")
        @test_throws "did not deliver a valid" G.download_validated(fileurl(html), dest; timeout = 10)
        @test readdir(dir) == ["error.html", "source.jp2"]
        G.download_validated(fileurl(jp2), dest; timeout = 10)
        @test read(dest) == read(jp2)
        @test length(readdir(dir)) == 3

        index = joinpath(dir, "gmt_data_server.txt")
        @test_throws "did not deliver a valid" G.download_validated(fileurl(html), index; timeout = 10)
        @test G.is_valid_download(G.BUNDLED_SERVER_INDEX)

        # a tile that passes the check but does not decode reports what the decoder said
        @test_throws "opj_decompress could not decode" G.decode_jp2(dest, 10)
        @test !isfile(dest * ".raw")
    end

    # La Palma: small, has both land and deep water, and a summit to check against
    limits = [-18.7, -17.1, 28.0, 29.2]

    if !topo_server_ok
        @test_skip "GMT data server unreachable: the download tests are skipped"
    else
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
            # pixel, which is what catches a transposed or wrongly assembled grid
            Topo = import_topo(limits, res = "03s")
            z = ustrip.(Topo.fields.Topography)
            i = argmax(z)
            @test isapprox(Topo.lon.val[i], -17.885, atol = 0.002)
            @test isapprox(Topo.lat.val[i], 28.754, atol = 0.002)
            @test isapprox(maximum(z) * 1000, 2414, atol = 2)      # metres
        end

        @testset "bathymetry" begin
            # 3 arcseconds is SRTM, which is land only -- the sea comes from the coarser set
            # the server names for it, so it must not be flat
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
end

@testset "import_topo offline behavior" begin
    G = GeophysicalModelGenerator
    dir = mktempdir()
    dest = joinpath(dir, "N00E000.SRTMGL3.jp2")
    touch(dest * ".absent")          # a remembered 404
    @test G.download_tile("earth_relief/earth_relief_03s_g/N00E000.SRTMGL3.jp2", dest) == false
    @test_throws "res must be one of" import_topo([10, 11, 40, 41], res = "07m")
    @test_throws "dataset \"earth_releif\" is not on the GMT data server" import_topo([10, 11, 40, 41], dataset = "earth_releif")
end
