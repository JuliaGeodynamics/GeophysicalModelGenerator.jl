# NOTE: these are useful routines that are only made available when the GMT package is already loaded in the REPL
module GMT_utils

import GeophysicalModelGenerator: import_GeoTIFF

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
using GMT

using GeophysicalModelGenerator: lonlatdepth_grid, GeoData, UTMData, km, remove_NaN_surface!

println("Loading GMT routines within GMG")


"""
  data_GMT = import_GeoTIFF(fname::String; fieldname=:layer1, negative=false, iskm=true, NorthernHemisphere=true, constantDepth=false, removeNaN_z=false, removeNaN_field=false)

This imports a GeoTIFF dataset (usually containing a surface of some sort) using GMT.
The file should either have `UTM` coordinates of `longlat` coordinates. If it doesn't, you can 
use QGIS to convert it to `longlat` coordinates.

Optional keywords:
- `fieldname` : name of the field (default=:layer1)
- `negative`  : if true, the depth is multiplied by -1 (default=false)
- `iskm`      : if true, the depth is multiplied by 1e-3 (default=true)
- `NorthernHemisphere`: if true, the UTM zone is set to be in the northern hemisphere (default=true); only relevant if the data uses UTM projection
- `constantDepth`: if true we will not warp the surface by z-values, but use a constant value instead
- `removeNaN_z`  : if true, we will remove NaN values from the z-dataset
"""
function import_GeoTIFF(fname::String; fieldname = :layer1, negative = false, iskm = true, NorthernHemisphere = true, constantDepth = false, removeNaN_z = false, removeNaN_field = false)
    G = gmtread(fname)

    # Transfer to GeoData
    nx, ny = length(G.x) - 1, length(G.y) - 1
    Lon, Lat, Depth = lonlatdepth_grid(G.x[1:nx], G.y[1:ny], 0)
    if hasfield(typeof(G), :z)
        Depth[:, :, 1] = G.z'
        if negative
            Depth[:, :, 1] = -G.z'
        end
        if iskm
            Depth *= 1.0e-3 * km
        end
    end

    # Create GeoData structure
    data = zero(Lon)
    if hasfield(typeof(G), :z)
        data = Depth

    elseif hasfield(typeof(G), :image)
        if length(size(G.image)) == 3
            data = permutedims(G.image, [2, 1, 3])
        elseif length(size(G.image)) == 2
            if size(G.image)==(nx,ny)
              data[:,:,1] = G.image
            elseif size(G.image)==(ny,nx)
              data[:,:,1] = G.image'
            else
              error("unknown size; ")
            end
        end

    end

    if removeNaN_z
        remove_NaN_surface!(Depth, Lon, Lat)
    end
    if removeNaN_field
        remove_NaN_surface!(data, Lon, Lat)
    end
    data_field = NamedTuple{(fieldname,)}((data,))

    if constantDepth
        Depth = zero(Lon)
    end

    if contains(G.proj4, "utm")
        zone = parse(Int64, split.(split(G.proj4, "zone=")[2], " ")[1])  # retrieve UTM zone
        data_GMT = UTMData(Lon, Lat, Depth, zone, NorthernHemisphere, data_field)

    elseif contains(G.proj4, "longlat")
        data_GMT = GeoData(Lon, Lat, Depth, data_field)

    else
        error("I'm sorry, I don't know how to handle this projection yet: $(G.proj4)\n
           We recommend that you transfer your GeoTIFF to longlat by using QGIS \n
           Open the GeoTIFF there and Export -> Save As , while selecting \"EPSG:4326 - WGS 84\" projection.")
    end

    return data_GMT
end


end
