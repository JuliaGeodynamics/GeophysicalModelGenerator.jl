using GeophysicalModelGenerator
using ParallelTestRunner

testsuite = find_tests(@__DIR__)

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
