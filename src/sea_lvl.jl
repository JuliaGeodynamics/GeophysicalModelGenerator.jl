export sea_level_files, SeaLevel, load_sea_level, curve_name

pkg_dir = pkgdir(GeophysicalModelGenerator)

const sea_level_path = joinpath(pkg_dir, joinpath("src", "sea_level_data"))

"Maps each built-in sea-level curve name to its data file."
const sea_level_files = Dict(
    :Spratt_800ka => "Spratt2016-800ka.txt",
    :Spratt_450ka => "Spratt2016-450ka.txt",
    :Bintanja_3Ma => "Bintanja-3Ma.txt",
    :Grant_153ka => "Grant2012_153ka.txt",
    :Rohl_Bint_3Ma => "Rohl-Bint-3Ma.txt",
    :Rohling2009_516ka => "Rohling2009-516ka.txt",
    :Siddall2003_379ka => "Siddall2003-379ka.txt",
    :Waelbroeck2002 => "Waelbroeck2002.txt",
)

"""
    SeaLevel(name::Symbol; flip_elevation=false, flip_age=false)

Sea-level curve `name` (a key of `sea_level_files`) with fields `elevation` and `age` (ka).
`flip_elevation` and `flip_age` must be given together: each reverses its vector, and
reversing only one would mispair ages and elevations.
"""
struct SeaLevel{T}
    elevation::Vector{T}
    age::Vector{T}
    name::Symbol

    function SeaLevel(name::Symbol; flip_elevation = false, flip_age = false)
        elevation, age = load_sea_level(name; flip_age, flip_elevation)
        return new{eltype(age)}(elevation, age, name)
    end
end

Base.getindex(x::SeaLevel, i::Int64) = x.elevation[i]
Base.getindex(x::SeaLevel, i::Int64, j::Int64) = x.elevation[i], x.age[j]
Base.getindex(x::SeaLevel, i::Tuple{Int64}) = x.elevation[i...], x.age[i...]
Base.size(x::SeaLevel) = size(x.elevation)
Base.eachindex(x::SeaLevel) = eachindex(x.elevation)
Base.axes(x::SeaLevel) = axes(x.elevation)
Base.length(x::SeaLevel) = length(x.elevation)
"Return the name of the sea-level curve `x`."
curve_name(x::SeaLevel) = x.name

"Load the built-in sea-level curve `name` (a key of `sea_level_files`) as a `SeaLevel`."
function load_sea_level(name::Symbol; flip_elevation = false, flip_age = false)
    flip_age == flip_elevation || error("flip_age and flip_elevation must be equal: flipping only one mispairs age and elevation.")
    fname = sea_level_files[name]
    data = readdlm(joinpath(sea_level_path, fname))
    age, h = data[:, 1], data[:, 2]
    flip_elevation && reverse!(h)
    flip_age && reverse!(age)
    return h, age
end
