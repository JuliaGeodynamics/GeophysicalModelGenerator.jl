# This allows creating ASAGI files from GMG structures, which can be used for example 
# as input for codes such as SeisSol or ExaHype
#
# see https://github.com/TUM-I5/ASAGI
#
#

using NCDatasets
using NCDatasets: nc_create, NC_NETCDF4, NC_CLOBBER, NC_NOWRITE, nc_def_dim, nc_def_compound, nc_insert_compound, nc_def_var, nc_put_var, nc_close, NC_INT, nc_unsafe_put_var, libnetcdf, check, ncType, nc_open, nc_inq_vartype, nc_inq_compound_nfields, nc_inq_compound_size, nc_inq_compound_name, nc_inq_compound_fieldoffset,nc_inq_compound_fieldndims,nc_inq_compound_fielddim_sizes, nc_inq_compound_fieldname, nc_inq_compound_fieldindex, nc_inq_compound_fieldtype, nc_inq_compound, nc_inq_varid, nc_get_var!, nc_insert_array_compound, nc_def_vlen

export write_ASAGI, read_ASAGI


"""
    write_ASAGI(fname::String, Data::CartData)
Writes a CartData structure to an ASAGI file, which can be read by SeisSol or ExaHype.
"""
function write_ASAGI(fname::String, Data::CartData)
    
    nx,ny,nz = size(Data.x)
    x = Data.x.val[:,1,1]
    y = Data.y.val[1,:,1]
    z = Data.z.val[1,1,:]

    # Transfer data to a single array with NamedTuple entries
    material = fields_to_namedtuple(Data.fields)

    fname_asagi = fname*"_ASAGI.nc"
    
    #ncid = nc_create(fname_asagi, NC_NETCDF4|NC_CLOBBER)
    ds = NCDataset(fname_asagi,"c", format=:netcdf4)

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
    typeid = nc_def_compound(ds.ncid, sizeof(T), "material")

    for i = 1:fieldcount(T)
        #local dim_sizes
        offset = fieldoffset(T,i)
        nctype = ncType[fieldtype(T,i)]
        nc_insert_compound(
                ds.ncid, typeid, fieldname(T,i),
                offset, nctype)
    end

    varid = nc_def_var(ds.ncid, "data", typeid, reverse(dimids))
    nc_put_var(ds.ncid, varid, material)
    
    # close file
    close(ds)

    return fname_asagi
end



# Transfer fields to a single array with NamedTuple entries
function fields_to_namedtuple(fields::NamedTuple)
    names   =   keys(fields)
    nfield  =   length(fields)
    ndim    =   length(size(fields[1]))

    s2      =   NamedTuple{names}(zeros(nfield))
    
    material = Array{typeof(s2),ndim}(undef, size(fields[1]))
    for I in eachindex(material)
        data_local = []
        for ifield in 1:nfield
           push!(data_local,fields[ifield][I])
        end

        local_tup = NamedTuple{names}(data_local)

        material[I] = local_tup
    end

    return material
end

module NCReconstructedTypes end


""" 
    data::CartData = read_ASAGI(fname_asagi::String)

This reads an ASAGI NetCDF file, which is used as input for a number of codes
"""
function read_ASAGI(fname_asagi::String)

    @assert fname_asagi[end-2:end]==".nc"

    ds = NCDataset(fname_asagi,"r")

    x = ds["x"][:];
    y = ds["y"][:];
    z = ds["z"][:];

    nx,ny,nz = length(x),length(y),length(z)
    
    data_set_names = keys(ds)
    id = findall(data_set_names.!="x" .&& data_set_names.!="y" .&& data_set_names.!="z")
    data_name = data_set_names[id]

    varid = nc_inq_varid(ds.ncid,data_name[1])
    xtype = nc_inq_vartype(ds.ncid,varid)
    
    # retrieve names of the fields
    numfields = nc_inq_compound_nfields(ds.ncid,xtype)
    cnames = Symbol.(nc_inq_compound_fieldname.(ds.ncid,xtype,0:(numfields-1)))
    
    types = []
    for fieldid = 0:(numfields-1)
        local dim_sizes
        fT = NCDatasets.jlType[nc_inq_compound_fieldtype(ds.ncid,xtype,fieldid)]

        fieldndims = nc_inq_compound_fieldndims(ds.ncid,xtype,fieldid)

        if fieldndims == 0
            push!(types,fT)
        else
            dim_sizes = nc_inq_compound_fielddim_sizes(ds.ncid,xtype,fieldid)
            fT2 = NTuple{Int(dim_sizes[1]),fT}
            push!(types,fT2)
        end
    end

    reconname = Symbol(string(nc_inq_compound_name(ds.ncid,xtype)))

    # would be better if this would be able to directly read the NamdTuple,
    # instead of creating a bogus struct
    Core.eval(
        NCReconstructedTypes,
        Expr(:struct, false, reconname,
             Expr(:block,
                  Any[ Expr(Symbol("::"), cnames[i], types[i]) for i = 1:length(types) ]...,
                  # suppress default constructors, plus a bogus `new()` call to make sure
                  # initialized is zero.
                  Expr(:if, false, Expr(:call, :new))
                  )))

    T2   = getfield(NCReconstructedTypes, reconname)
    data = Array{T2,3}(undef,nx,ny,nz)
    nc_get_var!(ds.ncid, varid, data)
    
    # at this stage we have an array with Main.NCReconstructedTypes
    # with the correct names and types
    #
    # Now we need to split them into different fields.

    read_fields_data = ()
    for ifield=1:numfields
        data_1 = zeros(types[ifield],nx,ny,nz)
        for I in CartesianIndices(data)
            loc = data[I]
            data_1[I] = getproperty(loc, cnames[ifield])
        end
        
        read_fields_data = (read_fields_data..., data_1)
    end
    read_fields = NamedTuple{(cnames...,)}(read_fields_data)

    close(ds)

    return CartData(xyz_grid(x,y,z)..., read_fields)
end