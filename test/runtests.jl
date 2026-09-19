using GeophysicalModelGenerator
using ParallelTestRunner

# Download the topography tiles the tests need *before* spawning the parallel
# workers, so they are served from GMT's cache instead of being downloaded
# concurrently (which times out regularly in CI). See topo_prefetch.jl.
include("topo_prefetch.jl")

testsuite = find_tests(@__DIR__)
# Not a test file, only a helper
delete!(testsuite, "topo_prefetch")

# Add `using GeophysicalModelGenerator` to each test
for (name, expr) in testsuite
    if name != "test_tutorials" # Tutorials are standalone
        testsuite[name] = quote
            using GeophysicalModelGenerator
            $expr
        end
    end
end

try
    ParallelTestRunner.runtests(GeophysicalModelGenerator, ARGS; testsuite)
finally
    # Cleanup
    foreach(f -> rm(joinpath(@__DIR__, f)), filter(endswith(".vts"), readdir(@__DIR__)))
    foreach(f -> rm(joinpath(@__DIR__, f)), filter(endswith(".vtu"), readdir(@__DIR__)))
    if isdir(joinpath(@__DIR__, "markers"))
        rm(joinpath(@__DIR__, "markers"), recursive = true)
    end
end
