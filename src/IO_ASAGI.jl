# This allows creating ASAGI files from GMG structures, which can be used for example 
# as input for codes such as SeisSol or ExaHype
#
# see https://github.com/TUM-I5/ASAGI
#
#

using NCDatasets
using NCDatasets: nc_create, NC_NETCDF4, NC_CLOBBER, NC_NOWRITE, nc_def_dim, nc_def_compound, nc_insert_compound, nc_def_var, nc_put_var, nc_close, NC_INT, nc_unsafe_put_var, libnetcdf, check, ncType, nc_open, nc_inq_vartype, nc_inq_compound_nfields, nc_inq_compound_size, nc_inq_compound_name, nc_inq_compound_fieldoffset,nc_inq_compound_fieldndims,nc_inq_compound_fielddim_sizes, nc_inq_compound_fieldname, nc_inq_compound_fieldindex, nc_inq_compound_fieldtype, nc_inq_compound, nc_inq_varid, nc_get_var!, nc_insert_array_compound, nc_def_vlen

export write_ASAGI


"""
Writes a CartData structure to an ASAGI file
"""
function write_ASAGI(fname::String, Data::CartData)
    
    nx,ny,nz = size(Data.x)
    x = Data.x.val[:,1,1]
    y = Data.y.val[1,:,1]
    z = Data.z.val[1,1,:]

    # Transfer data to a single array with NamedTuple entries
    material = fields_to_namedtuple(Data.fields)

    fname_asagi = fname*"_ASAGI1.nc"
    #ncid = nc_create(fname_asagi, NC_NETCDF4|NC_CLOBBER)
    ds = NCDataset(fname_asagi,"c")

    # Write dimensions
    x_dimid = nc_def_dim(ds.ncid, "nx", nx)
    y_dimid = nc_def_dim(ds.ncid, "ny", ny)
    z_dimid = nc_def_dim(ds.ncid, "nz", nz)

    v_x = defVar(ds,"x",Float32,("nx",))
    v_y = defVar(ds,"y",Float32,("ny",))
    v_z = defVar(ds,"z",Float32,("nz",))
    v_x[:] = x
    v_y[:] = y
    v_z[:] = z

    # add Tuple with data to file
    dimids = [x_dimid, y_dimid, z_dimid]
    T = eltype(material)
    typeid = nc_def_compound(ds.ncid, sizeof(T), "sample_compound_type")

    for i = 1:fieldcount(T)
        local dim_sizes
        offset = fieldoffset(T,i)
        nctype = ncType[fieldtype(T,i)]
        nc_insert_compound(
                ds.ncid, typeid, fieldname(T,i),
                offset, nctype)
    end

    varid = nc_def_var(ds.ncid, "material", typeid, reverse(dimids))
    nc_put_var(ds.ncid, varid, material)
    
    # close file
    close(ds)

end



# Transfer fields to a single array with NamedTuple entries
function fields_to_namedtuple(fields::NamedTuple)
    names   =   keys(Data.fields)
    nfield  =   length(fields)
    ndim    =   length(size(fields[1]))

    s2      =   NamedTuple{names}(zeros(nfield))
    
    material = Array{typeof(s2),ndim}(undef, size(Data.fields[1]))
    for I in eachindex(material)
        data_local = []
        for ifield in 1:nfield
           push!(data_local,Data.fields[ifield][I])
        end

        local_tup = NamedTuple{names}(data_local)

        material[I] = local_tup
    end

    return material
end