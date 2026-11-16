# This provides functions to load and save GeoData structs & friends to file
using JLD2, Downloads

export save_GMG, load_GMG, download_data

"""
    save_GMG(filename::String, data::Union{GeoData, CartDat, UTMData}; dir=pwd())
Saves the dataset `data` to a JLD2 `file` (name without extension) in the directory `dir`

Example
====
```julia
julia> Lon3D,Lat3D,Depth3D = lonlatdepth_grid(1.0:3:10.0, 11.0:4:20.0, (-20:5:-10)*km);
julia> Data_set    =   GeophysicalModelGenerator.GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Depth3D,))   
julia> save_GMG("test",Data_set)
```
"""
function save_GMG(filename::String, data::Union{GeoData, CartData, UTMData}; dir = pwd())
    file_ext = joinpath(dir, filename * ".jld2")
    jldsave(file_ext; data)

    return nothing
end

"""
    load_GMG(filename::String, dir=pwd(); maxattempts=5)

Loads a `GeoData`/`CartData`/`UTMData` data set from jld2 file `filename`
Note: the `filename` can also be a remote `url`, in which case the file is downloaded to a temporary file, loaded and the temporary file removed again.
We make `maxattempts` attempts to download it before giving up.

Example 1 - Load local file
====
```julia
julia> data = load_GMG("test")
GeoData 
  size      : (4, 3, 3)
  lon       ϵ [ 1.0 : 10.0]
  lat       ϵ [ 11.0 : 19.0]
  depth     ϵ [ -20.0 : -10.0]
  fields    : (:DataFieldName,)
  attributes: ["note"]
```

Example 2 - remote download
====
```julia
julia> url  = "https://seafile.rlp.net/f/10f867e410bb4d95b3fe/?dl=1"
julia> load_GMG(url)
GeoData 
  size      : (149, 242, 1)
  lon       ϵ [ -24.875 : 35.375]
  lat       ϵ [ 34.375 : 71.375]
  depth     ϵ [ -11.76 : -34.7]
  fields    : (:MohoDepth,)
  attributes: ["author", "year"]
```

"""
function load_GMG(filename::String, dir = pwd(); maxattempts = 5)
    if is_url(filename)
        tmp = tempname() * ".jld2"
        try
            download_data(filename, tmp; maxattempts)
            return load_object(tmp)
        finally
            rm(tmp; force = true)
        end
    end
    return load_object(joinpath(dir, filename * ".jld2"))
end

is_url(name::AbstractString) = startswith(name, r"(https?|file)://")

# Rejects empty files and the HTML error pages that servers deliver with status 200. Bytes outside
# ASCII, which cannot start an HTML page, are replaced so that binary files can be compared as text.
function is_data_file(file::AbstractString)
    head = lowercase(String(map(b -> b < 0x80 ? b : UInt8('?'), open(io -> read(io, 15), file))))
    return !isempty(head) && !startswith(head, "<!doctype html") && !startswith(head, "<html")
end


"""
    download_data(url::String, local_filename="temp.dat"; dir=pwd(), maxattempts=5 )

Downloads a remote dataset with name `url` from a remote location and saves it to the current directory.
The destination is replaced only by a download that is not empty and not an HTML page.
If download fails, we make `maxattempts` attempts before giving up.

Example
====
```julia
julia> url  = "https://seafile.rlp.net/f/10f867e410bb4d95b3fe/?dl=1";
julia> download_data(url)
"/Users/kausb/.julia/dev/GeophysicalModelGenerator/temp.dat"
```

"""
function download_data(url::String, local_filename = "temp.dat"; dir = pwd(), maxattempts = 5)

    if !is_url(url)
        @warn "the url does not start with http(s):// or file://"
    end

    dest = joinpath(dir, local_filename)
    for attempt in 1:maxattempts
        try
            return download_validated(url, dest; valid = is_data_file)
        catch e
            e isa InterruptException && rethrow()
            @warn "Failed downloading $url on attempt $attempt/$maxattempts: $e"
            if attempt < maxattempts
                sleep(5 * attempt)  # back off progressively
            end
        end
    end
    error("Could not download data from $url after $maxattempts attempts")
end
