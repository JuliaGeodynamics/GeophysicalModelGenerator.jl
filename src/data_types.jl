# This is data_types.jl
# contains type definitions to be used in GeophysicalModelGenerator


# data structure for a list of values
mutable struct ValueList
    name::String
    unit::String
    values::Vector{Float64}
end

# data structure for lan/lot/depth data together with assigned values
mutable struct GeoData
    lat::ValueList
    lon::ValueList
    depth::ValueList
    values::tuple
end

mutable struct CartData
    x::ValueList
    y::ValueList
    z::ValueList
    values::tuple
end