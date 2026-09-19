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
Note: the `filename` can also be a remote `url`, in which case we first download that file to a temporary directory before opening it.
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

    local_filename = "download_GMG_temp_$(getpid()).jld2"
    if contains(filename, "http")
        file_ext = download_data(filename, local_filename, dir = dir, maxattempts = maxattempts)
    else
        # local file
        file_ext = joinpath(dir, filename * ".jld2")
    end

    # load data:
    data = load_object(file_ext)

    # remove local temporary file
    if contains(filename, "http")
        rm(joinpath(dir, local_filename), force = true)
    end

    return data
end


"""
    download_data(url::String, local_filename="temp.dat"; dir=pwd(), maxattempts=5 )

Downloads a remote dataset with name `url` from a remote location and saves it to the current directory.
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

    if !contains(url, "http")
        @warn "the url does not contain http; please double check that it worked"
    end

    #download remote file to a local temporary directory
    file_ext = []
    for attempt in 1:maxattempts
        try
            file_ext = Downloads.download(url, joinpath(dir, local_filename))
            # A truncated transfer can still return a path, so only accept a
            # download that actually produced a non-empty file.
            if filesize(file_ext) == 0
                error("downloaded file is empty")
            end
            break
        catch e
            file_ext = []
            @warn "Failed downloading $url on attempt $attempt/$maxattempts: $e"
            if attempt < maxattempts
                sleep(5 * attempt)  # back off progressively
            end
        end
    end
    if isempty(file_ext)
        error("Could not download data from $url after $maxattempts attempts")
    end

    return file_ext
end
