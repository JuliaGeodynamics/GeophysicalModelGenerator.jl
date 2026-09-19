# Warm the GMT topography cache before the parallel test workers start.
#
# Several tests (test_GMT, test_WaterFlow, the LaPalma tutorial) download
# topography tiles through GMT. Under ParallelTestRunner they run concurrently,
# so they all hit the GMT data server at the same time, and the transfers then
# regularly time out in CI ("Libcurl Error: Timeout was reached").
#
# Downloading the tiles here, one at a time and with generous retries, means the
# workers later find them in GMT's cache (~/.gmt/server) and never touch the
# network. The download itself is still exercised: it happens right here, and
# on CI the cache is bypassed for the scheduled run so it is tested regularly.
#
# NOTE: .github/workflows/CI.yml caches ~/.gmt keyed on the hash of this file,
# so changing the list below invalidates that cache automatically.

using GeophysicalModelGenerator, GMT

# (region, file) pairs, matching the calls made in the tests/tutorials
const TOPO_PREFETCH = [
    ([8.0, 9.0, 50.0, 51.0], "@earth_relief_01m"),               # test_GMT
    ([6.5, 7.3, 50.2, 50.6], "@earth_relief_03s"),               # test_WaterFlow
    ([-18.2, -17.5, 28.4, 29.0], "@earth_relief_15s.grd"),       # LaPalma tutorial
]

function prefetch_topography()
    gmt_cache = joinpath(get(ENV, "GMT_USERDIR", joinpath(homedir(), ".gmt")), "server")
    n_before = isdir(gmt_cache) ? count(f -> endswith(f, ".nc"), collect(Iterators.flatten(map(t -> t[3], walkdir(gmt_cache))))) : 0

    for (limits, file) in TOPO_PREFETCH
        t = @elapsed try
            # copy: import_topo may modify `limits` in place (negative longitudes)
            import_topo(copy(limits); file = file, maxattempts = 8)
        catch e
            # Do not abort the test run: the test that needs this tile will try
            # again itself and report a proper error if the server really is down.
            @warn "Could not prefetch topography $file for $limits; the test will retry" exception = (e, catch_backtrace())
        end
        @info "Prefetched topography $file for $limits in $(round(t, digits = 1)) s"
    end

    n_after = isdir(gmt_cache) ? count(f -> endswith(f, ".nc"), collect(Iterators.flatten(map(t -> t[3], walkdir(gmt_cache))))) : 0
    @info "GMT topography cache: $n_after tiles in $gmt_cache ($(n_after - n_before) newly downloaded)"
    return nothing
end

prefetch_topography()
