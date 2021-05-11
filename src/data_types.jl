# This is data_types.jl
# contains type definitions to be used in GeophysicalModelGenerator


# 1D data structure
struct ScatteredPoints
    name::String
    unit::String
    values::Vector{Float64}
end