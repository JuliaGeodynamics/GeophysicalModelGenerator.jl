using Test, GMT

# The topography comes from the tile cache that topo_prefetch.jl filled. If the prefetch could
# not reach the GMT data server, asking again here is pointless: skip on a push/PR run, insist
# on the scheduled run so a real breakage is caught.
topo_server_ok = get(ENV, "GMG_TOPO_PREFETCH_OK", "true") == "true" ||
    get(ENV, "GITHUB_EVENT_NAME", "") == "schedule"

if !topo_server_ok
    @test_skip "GMT data server unreachable: water flow tests skipped"
else
    # Some topographic data
    Topo = import_topo([6.5, 7.3, 50.2, 50.6], file = "@earth_relief_03s")

    # Flow the water through the area:
    Topo_water, sinks, pits, bnds = waterflows(Topo)

    @test maximum(Topo_water.fields.area) ≈ 9.309204547276944e8
    @test sum(Topo_water.fields.c) == 834501044
    @test sum(Topo_water.fields.nin) == 459361
    @test sum(Topo_water.fields.dir) == 2412566

    # With rain in m3/s per cell
    rainfall = ones(size(Topo.lon.val[:, :, 1])) * 1.0e-3 # 2D array with rainfall per cell area
    Topo_water1, sinks, pits, bnds = waterflows(Topo, rainfall = rainfall)

    @test maximum(Topo_water1.fields.area) ≈ 169.79800000000208
end
