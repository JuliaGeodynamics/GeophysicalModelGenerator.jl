module NCDatasets_utils

using NCDatasets

using NCDatasets: nc_create, NC_NETCDF4, NC_CLOBBER, NC_NOWRITE, nc_def_dim, nc_def_compound, nc_insert_compound, nc_def_var, nc_put_var, nc_close, NC_INT, nc_unsafe_put_var, libnetcdf, check, ncType, nc_open, nc_inq_vartype, nc_inq_compound_nfields, nc_inq_compound_size, nc_inq_compound_name, nc_inq_compound_fieldoffset, nc_inq_compound_fieldndims, nc_inq_compound_fielddim_sizes, nc_inq_compound_fieldname, nc_inq_compound_fieldindex, nc_inq_compound_fieldtype, nc_inq_compound, nc_inq_varid, nc_get_var!, nc_insert_array_compound, nc_def_vlen

import GeophysicalModelGenerator: write_ASAGI, read_ASAGI, tomo_2_GeoData, read_topo_netcdf
using GeophysicalModelGenerator: CartData, GeoData, lonlatdepth_grid, xyz_grid

function write_ASAGI(
        fname::String, Data::CartData;
        fields::Union{Nothing, Tuple} = nothing,
        km_to_m::Bool = false
    )

    nx, ny, nz = size(Data.x)
    x = Data.x.val[:, 1, 1]
    y = Data.y.val[1, :, 1]
    z = Data.z.val[1, 1, :]
    if km_to_m == true
        println("convert to meters")
        x = x .* 1000
        y = y .* 1000
        z = z .* 1000
    end

    # Transfer data to a single array with NamedTuple entries
    material = fields_to_namedtuple(Data.fields, fields)

    fname_asagi = fname * "_ASAGI.nc"

    ds = NCDataset(fname_asagi, "c", format = :netcdf4)

    # Write dimensions
    x_dimid = nc_def_dim(ds.ncid, "x", nx)
    y_dimid = nc_def_dim(ds.ncid, "y", ny)
    z_dimid = nc_def_dim(ds.ncid, "z", nz)

    v_x = defVar(ds, "x", eltype(x), ("x",))
    v_y = defVar(ds, "y", eltype(x), ("y",))
    v_z = defVar(ds, "z", eltype(x), ("z",))
    v_x[:] = x
    v_y[:] = y
    v_z[:] = z

    # add Tuple with data to file
    dimids = [x_dimid, y_dimid, z_dimid]
    T = eltype(material)
    typeid = nc_def_compound(ds.ncid, sizeof(T), "material")

    for i in 1:fieldcount(T)
        offset = fieldoffset(T, i)
        nctype = ncType[fieldtype(T, i)]
        nc_insert_compound(
            ds.ncid, typeid, fieldname(T, i),
            offset, nctype
        )
    end

    varid = nc_def_var(ds.ncid, "data", typeid, reverse(dimids))
    nc_put_var(ds.ncid, varid, material)

    # close file
    close(ds)

    return fname_asagi
end


# Transfer fields to a single array with NamedTuple entries
function fields_to_namedtuple(fields::NamedTuple, selected_fields)
    names = keys(fields)
    if !isnothing(selected_fields)
        names = selected_fields
    end
    for n in names
        haskey(fields, n) || error("Field $n is not in the data set; available fields: $(keys(fields))")
        isa(fields[n], Array) ||
            error("Field $n is not an Array but instead a $(typeof(fields[n])); only Arrays are supported")
    end

    arrays = map(n -> fields[n], Tuple(names))
    allequal(map(size, arrays)) || error("The fields $(Tuple(names)) must have the same size, got $(map(size, arrays))")

    # element type of each field is kept (e.g. Int32 phases stay Int32)
    return records(NamedTuple{Tuple(names), Tuple{map(eltype, arrays)...}}, arrays...)
end

# The inverse of `records`: one array per field of the NamedTuple records in `data`
split_records(data::AbstractArray{NamedTuple{N, T}}) where {N, T} =
    NamedTuple{N}(ntuple(i -> map(r -> getfield(r, i), data), Val(length(N))))

# One record of type `T` per point; `T` is a type parameter so that the loop is compiled for it
records(::Type{T}, arrays...) where {T} = map((x...) -> T(x), arrays...)


function read_ASAGI(fname_asagi::String)

    @assert fname_asagi[(end - 2):end] == ".nc"

    ds = NCDataset(fname_asagi, "r")

    x = ds["x"][:]
    y = ds["y"][:]
    z = ds["z"][:]

    nx, ny, nz = length(x), length(y), length(z)

    data_set_names = keys(ds)
    id = findall(data_set_names .!= "x" .&& data_set_names .!= "y" .&& data_set_names .!= "z")
    data_name = data_set_names[id]

    varid = nc_inq_varid(ds.ncid, data_name[1])
    xtype = nc_inq_vartype(ds.ncid, varid)

    # retrieve names of the fields
    numfields = nc_inq_compound_nfields(ds.ncid, xtype)
    cnames = Symbol.(nc_inq_compound_fieldname.(ds.ncid, xtype, 0:(numfields - 1)))

    types = []
    for fieldid in 0:(numfields - 1)
        local dim_sizes
        fT = NCDatasets.jlType[nc_inq_compound_fieldtype(ds.ncid, xtype, fieldid)]

        fieldndims = nc_inq_compound_fieldndims(ds.ncid, xtype, fieldid)

        if fieldndims == 0
            push!(types, fT)
        else
            dim_sizes = nc_inq_compound_fielddim_sizes(ds.ncid, xtype, fieldid)
            fT2 = NTuple{Int(dim_sizes[1]), fT}
            push!(types, fT2)
        end
    end

    # Create a single NamedTuple with correct type and names
    data_element = ()
    for ifield in 1:numfields
        data_element = (data_element..., types[ifield](0.0))
    end

    T2 = typeof(NamedTuple{(cnames...,)}(data_element))
    data = Array{T2, 3}(undef, nx, ny, nz)
    nc_get_var!(ds.ncid, varid, data)

    # data is an array of NamedTuples with the names and types of the fields; split it into one array per field
    read_fields = split_records(data)

    close(ds)

    return CartData(xyz_grid(x, y, z)..., read_fields)
end


function tomo_2_GeoData(filename::String; vel_type::String = "vs")
    depth, lon, lat, vel = NCDataset(filename) do data
        vel = data[vel_type][:, :, :]
        T = float(nonmissingtype(eltype(vel)))
        (
            nomissing(data["depth"][:]), nomissing(data["longitude"][:]), nomissing(data["latitude"][:]),
            map(v -> ismissing(v) ? T(NaN) : T(v), vel),
        )
    end

    # create lon, lat, depth grid
    Lon, Lat, Depth = lonlatdepth_grid(lon, lat, .- depth)

    return GeoData(Lon, Lat, Depth, (vel = vel,))
end


"""
    read_topo_netcdf(file)

Read one of the global `.grd` grids that the GMT data server publishes for its coarse
resolutions, and return `(lon, lat, z)`. The tiled resolutions do not come through here --
they are JPEG2000 and are decoded without NCDatasets.
"""
function read_topo_netcdf(file::AbstractString)
    return NCDataset(file) do ds
        lon = Array(ds["lon"][:])
        lat = Array(ds["lat"][:])
        z = ds["z"][:, :]
        any(ismissing, z) && error("$file has missing elevations")
        # the file stores z as (lon, lat); the rest of the code wants (lat, lon)
        (lon, lat, permutedims(nomissing(z)))
    end
end


end
