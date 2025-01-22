export sea_level_files, SeaLevel, load_sea_level, curve_name

const sea_level_path = "sea_level_data"

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

struct SeaLevel{T}
    elevation::Vector{T}
    age::Vector{T}
    name::Symbol

    function SeaLevel(name::Symbol; flip_elevation = false, flip_age = false)
        age, elevation = load_sea_level(
            name;
            flip_age = flip_age,
            flip_elevation = flip_elevation
        )
        return new{eltype(age)}(age, elevation, name)
    end
end

Base.getindex(x::SeaLevel, i::Int64) = x.elevation[i]
Base.getindex(x::SeaLevel, i::Int64, j::Int64) = x.elevation[i], x.age[j]
Base.getindex(x::SeaLevel, i::Tuple{Int64}) = x.elevation[i...], x.age[i...]
Base.size(x::SeaLevel) = size(x.elevation)
Base.eachindex(x::SeaLevel) = eachindex(x.elevation)
Base.axes(x::SeaLevel) = axes(x.elevation)
Base.length(x::SeaLevel) = length(x.elevation)
curve_name(x::SeaLevel) = x.name

function load_sea_level(name::Symbol; flip_elevation = false, flip_age = false)
    fname = sea_level_files[name]
    data = readdlm(joinpath(sea_level_path, fname))
    h = data[:, 1]
    age = data[:, 2]
    flip_elevation && reverse!(h)
    flip_age && reverse!(age)
    return h, age
end
