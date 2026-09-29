using Test
using GeophysicalModelGenerator, GMT

# The topography comes from the tile cache that topo_prefetch.jl filled. If the prefetch could
# not reach the GMT data server, asking again here is pointless: skip on a push/PR run, insist
# on the scheduled run so a real breakage is caught.
topo_server_ok = get(ENV, "GMG_TOPO_PREFETCH_OK", "true") == "true" ||
    get(ENV, "GITHUB_EVENT_NAME", "") == "schedule"

if !topo_server_ok
    @test_skip "GMT data server unreachable: topography tests skipped"
else
    Topo = import_topo(lon = [8, 9], lat = [50, 51])
    @test sum(Topo.depth.val) ≈ 1076.7045 rtol = 1.0e-1

    Topo = import_topo([8, 9, 50, 51])
    @test sum(Topo.depth.val) ≈ 1076.7045 rtol = 1.0e-1
end

test_fwd = import_GeoTIFF("test_files/length_fwd.tif", fieldname = :forward)
@test  maximum(test_fwd.fields.forward) ≈ 33.17775km

test2 = import_GeoTIFF("test_files/UTM2GTIF.TIF")
@test   test2.fields.layer1[20, 20] == 233.0
