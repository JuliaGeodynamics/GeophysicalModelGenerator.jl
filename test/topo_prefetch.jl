# Warm the topography cache before any test runs.
#
# This is the one place in the test suite that downloads topography. Every region and
# resolution that a test or tutorial asks for is listed below, so once this has run the
# tests find their tiles in the cache and never touch the network. That matters because the
# GMT data server throttles when several jobs pull tiles at once, which is exactly what the
# parallel test workers (and the CI matrix) would otherwise do.
#
# The download is still tested: it happens right here. On CI the tile cache is restored for
# push/PR runs, so this normally finds everything in place; the scheduled run clears the
# cache first, so the download is exercised for real every two weeks.
#
# If the server cannot be reached, that is recorded in ENV["GMG_TOPO_PREFETCH_OK"], which
# the parallel workers inherit: the tests that need topography then skip themselves rather
# than fail (except on the scheduled run, where they fail, so a real breakage is caught).
# `import_topo` also remembers an unreachable server for a while, so anything that does ask
# again fails in seconds instead of waiting out the full retry budget.
#
# NOTE: .github/workflows/CI.yml caches the scratch space that `import_topo` keeps its tiles
# in, keyed on the hash of this file, so changing the list below invalidates that cache
# automatically. The previous tiles are still restored; only the new ones are downloaded.

using GeophysicalModelGenerator

# (region, keyword arguments) pairs, matching the calls made in the tests and tutorials
const TOPO_PREFETCH = [
    ([8.0, 9.0, 50.0, 51.0], (file = "@earth_relief_01m",)),                # test_GMT
    ([6.5, 7.3, 50.2, 50.6], (file = "@earth_relief_03s",)),                # test_WaterFlow
    ([-18.2, -17.5, 28.4, 29.0], (file = "@earth_relief_15s.grd",)),        # LaPalma tutorial
    ([-18.7, -17.1, 28.0, 29.2], (res = "01m",)),                           # test_import_topo ...
    ([-18.7, -17.1, 28.0, 29.2], (res = "03s",)),                           # (also fetches the 15s filler)
    ([-18.7, -17.1, 28.0, 29.2], (res = "30s",)),
    ([-18.7, -17.1, 28.0, 29.2], (dataset = "earth_gebco", res = "15s")),
]

function prefetch_topography()
    cache = GeophysicalModelGenerator.topo_cache_dir()
    count_tiles() = isdir(cache) ? count(f -> endswith(f, ".jp2") || endswith(f, ".grd"), readdir(cache)) : 0
    n_before = count_tiles()

    # a stale "server unreachable" marker must not stop the one attempt this run makes
    GeophysicalModelGenerator.reset_topo_server()

    ok = true
    for (limits, kwargs) in TOPO_PREFETCH
        t = @elapsed fetched = try
            import_topo(copy(limits); kwargs..., maxattempts = 3)
            true
        catch e
            @warn "Could not prefetch topography $kwargs for $limits" exception = (e, catch_backtrace())
            false
        end
        ok &= fetched
        fetched && @info "Prefetched topography $kwargs for $limits in $(round(t, digits = 1)) s"
    end

    n_after = count_tiles()
    @info "Topography cache: $n_after tiles in $cache ($(n_after - n_before) newly downloaded)"
    ok || @warn "The GMT data server could not be reached; tests that need topography will be skipped"
    ENV["GMG_TOPO_PREFETCH_OK"] = string(ok)
    return ok
end

prefetch_topography()
