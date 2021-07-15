using Base: Int64, Float64, NamedTuple
using Printf

# LaMEM I/O
# 
# These are routines that help to create a LaMEM model setup from a CartData structure

export ParseValue_LaMEM_InputFile, LaMEM_grid, ReadLaMEM_InputFile
export Save_LaMEMMarkersParallel, GetProcessorPartitioning

"""
    Structure that holds information about the LaMEM grid (usually read from an input file)
"""
struct LaMEM_grid
    nmark_x :: Int64
    nmark_y :: Int64
    nmark_z :: Int64
    nump_x  :: Int64 
    nump_y  :: Int64
    nump_z  :: Int64
    
    nel_x   :: Int64
    nel_y   :: Int64
    nel_z   :: Int64

    W       ::  Float64
    L       ::  Float64
    H       ::  Float64

    coord_x 
    coord_y 
    coord_z

    x1D_c
    y1D_c
    z1D_c

    X 
    Y 
    Z

   
end

""" 
    CartData(Grid::LaMEM_grid, fields::NamedTuple)

Creates a `CartData` struct from a LaMEM grid and from fields stored on that grid. Note that one needs to have a field `Phases` and optionally a field `Temp` to create LaMEM marker files.
"""
CartData(Grid::LaMEM_grid, fields::NamedTuple) = CartData(Grid.X, Grid.Y, Grid.Z, fields)


"""
    value = ParseValue_LaMEM_InputFile(file,keyword,type)

Extracts a certain `keyword` from a LaMEM input `file` and convert it to a certain type 

#Example
```julia
julia> nmark_z = ParseValue_LaMEM_InputFile("SaltModels.dat","nmark_z",Int64)
```
"""
function ParseValue_LaMEM_InputFile(file,keyword,type)
    value = nothing
    for line in eachline(file)
        line_strip = lstrip(line)       # strip leading tabs/spaces

        # Strip comments
        ind        = findfirst("#", line)
        if isnothing(ind)
            # no comments
        else
            line_strip = line_strip[1:ind[1]-2];
        end
        line_strip = rstrip(line_strip)       # strip last tabs/spaces

        if startswith(line_strip, keyword)
            ind = findfirst("=", line_strip)
            if type==String
                value = split(line_strip)[3:end]
            else
                value = parse.(type,split(line_strip)[3:end])

                if length(value)==1
                    value=value[1];
                end
            end
        end
        
    end

    return value
end

"""
    Grid::LaMEM_grid = ReadLaMEM_InputFile(file) 

Parses a LaMEM input file and stores grid information in the `Grid` structure

#Example
```julia
julia> Grid = ReadLaMEM_InputFile("SaltModels.dat") 
LaMEM Grid: 
nel         : (32, 32, 32)
marker/cell : (3, 3, 3)
markers     : (96, 96, 96)
x           ϵ [-3.0 : 3.0]
y           ϵ [-2.0 : 2.0]
z           ϵ [-2.0 : 0.0]
"""
function ReadLaMEM_InputFile(file)

    nmark_x = ParseValue_LaMEM_InputFile(file,"nmark_x",Int64)
    nmark_y = ParseValue_LaMEM_InputFile(file,"nmark_y",Int64)
    nmark_z = ParseValue_LaMEM_InputFile(file,"nmark_z",Int64)

    nel_x   = ParseValue_LaMEM_InputFile(file,"nel_x",Int64)
    nel_y   = ParseValue_LaMEM_InputFile(file,"nel_y",Int64)
    nel_z   = ParseValue_LaMEM_InputFile(file,"nel_z",Int64)

    coord_x = ParseValue_LaMEM_InputFile(file,"coord_x",Float64)
    coord_y = ParseValue_LaMEM_InputFile(file,"coord_y",Float64)
    coord_z = ParseValue_LaMEM_InputFile(file,"coord_z",Float64)

    if (length(coord_x)>2) || (length(coord_y)>2) || (length(coord_z)>2)
        error("Routine currently not working for variable grid spacing")
    end

    W       = coord_x[end]-coord_x[1];
    L       = coord_y[end]-coord_y[1];
    H       = coord_z[end]-coord_z[1];

    nump_x  = nel_x*nmark_x;
    nump_y  = nel_y*nmark_y;
    nump_z  = nel_z*nmark_z;

    dx      =   W/nump_x;
    dy      =   L/nump_y;
    dz      =   H/nump_z;

    # these lines should be replaced with a separate routine for variable spacing   
    x       =   coord_x[1]+dx/2: dx : coord_x[end]-dx/2;
    y       =   coord_y[1]+dy/2: dy : coord_y[end]-dy/2;
    z       =   coord_z[1]+dz/2: dz : coord_z[end]-dz/2;

    X,Y,Z   =   LonLatDepthGrid(x,y,z); # create 3D grid using regular spacng
    
    Grid    =  LaMEM_grid(  nmark_x,    nmark_y,    nmark_z,
                            nump_x,     nump_y,     nump_z,
                            nel_x,      nel_y,      nel_z,    
                            W,          L,          H,
                            coord_x,    coord_y,    coord_z,
                            x,          y,          z,
                            X,          Y,          Z);

    return Grid
end

# Print an overview of the LaMEM Grid struct:
function Base.show(io::IO, d::LaMEM_grid)
    println(io,"LaMEM Grid: ")
    println(io,"  nel         : ($(d.nel_x), $(d.nel_y), $(d.nel_z))")
    println(io,"  marker/cell : ($(d.nmark_x), $(d.nmark_y), $(d.nmark_z))")
    println(io,"  markers     : ($(d.nump_x), $(d.nump_x), $(d.nump_x))")
    println(io,"  x           ϵ [$(d.coord_x[1]) : $(d.coord_x[2])]")
    println(io,"  y           ϵ [$(d.coord_y[1]) : $(d.coord_y[2])]")
    println(io,"  z           ϵ [$(d.coord_z[1]) : $(d.coord_z[2])]")
end

"""
    Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers")

Saves a LaMEM marker file from the CartData structure `Grid`. It must have a field called `Phases`, holding phase information (as integers) and optionally a field `Temp` with temperature info. 
It is possible to provide a LaMEM partitioning file `PartitioningFile`. If not, output is assumed to be for one processor.

The size of `Grid` should be consistent with what is provided in the LaMEM input file. In practice, the size of the mesh can be retrieved from a LaMEM input file using `ReadLaMEM_InputFile`.

# Example 

```
julia> Grid    = ReadLaMEM_InputFile("LaMEM_input_file.dat")
julia> Phases  = zeros(Int32,size(Grid.X));
julia> Temp    = ones(Float64,size(Grid.X));
julia> Model3D = CartData(Grid, (Phases=Phases,Temp=Temp))
julia> Save_LaMEMMarkersParallel(Model3D)
Writing LaMEM marker file -> ./markers/mdb.00000000.dat
```    
If you want to create a LaMEM input file for multiple processors:
```
julia> Save_LaMEMMarkersParallel(Model3D, PartitioningFile="ProcessorPartitioning_4cpu_1.2.2.bin")
Writing LaMEM marker file -> ./markers/mdb.00000000.dat
Writing LaMEM marker file -> ./markers/mdb.00000001.dat
Writing LaMEM marker file -> ./markers/mdb.00000002.dat
Writing LaMEM marker file -> ./markers/mdb.00000003.dat
```

"""
function Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers")

    x = ustrip(Grid.x.val[:,1,1]);
    y = ustrip(Grid.y.val[1,:,1]);
    z = ustrip(Grid.z.val[1,1,:]);
    
    if haskey(Grid.fields,:Phases)
        Phases = Grid.fields[:Phases];
    else
        error("You must provide the field :Phases in the structure")
    end
    
    if haskey(Grid.fields,:Temp)
        Temp = Grid.fields[:Temp];
    else
        println("Field :Temp is not provided; setting it to zero")
        Temp = zeros(size(Phases));
    end
    
    if PartitioningFile==empty
        # in case we run this on 1 processor only
        Nprocx  =   1;
        Nprocy  =   1;
        Nprocz  =   1;
        xc,yc,zc = x,y,z;
    else
        Nprocx,Nprocy,Nprocz, 
        xc,yc,zc, 
        nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(PartitioningFile)
    end

    Nproc                       =   Nprocx*Nprocy*Nprocz;
    num, num_i, num_j, num_k    =   get_numscheme(Nprocx, Nprocy, Nprocz);

    xi,ix_start,ix_end          =   get_ind(x,xc,Nprocx);
    yi,iy_start,iy_end          =   get_ind(y,yc,Nprocy);
    zi,iz_start,iz_end          =   get_ind(z,zc,Nprocz);

    x_start                     =   ix_start[num_i[:]];
    y_start                     =   iy_start[num_j[:]];
    z_start                     =   iz_start[num_k[:]];
    x_end                       =   ix_end[num_i[:]];
    y_end                       =   iy_end[num_j[:]];
    z_end                       =   iz_end[num_k[:]];

    # Loop over all processors partition
    for n=1:Nproc
        # Extract coordinates for current processor
        
        part_x   = Grid.x.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]];
        part_y   = Grid.y.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]];
        part_z   = Grid.z.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]];
        part_phs = Phases[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]];
        part_T   =   Temp[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]];
        num_particles = size(part_x,1)* size(part_x,2) * size(part_x,3);

        # Information vector per processor
        num_prop        =   5;      # number of properties we save [x/y/z/phase/T]
        lvec_info       =   num_particles;
    
        lvec_prtcls     =   zeros(Float64,num_prop*num_particles);
    
        lvec_prtcls[1:num_prop:end] = part_x[:];
        lvec_prtcls[2:num_prop:end] = part_y[:];
        lvec_prtcls[3:num_prop:end] = part_z[:];
        lvec_prtcls[4:num_prop:end] = part_phs[:];
        lvec_prtcls[5:num_prop:end] = part_T[:];

        # Write output files
        if ~isdir(directory); mkdir(directory); end         # Create dir if not existent
        fname = @sprintf "%s/mdb.%1.8d.dat"  directory (n-1);   # Name
        println("Writing LaMEM marker file -> $fname")                   # print info
        lvec_output    = [lvec_info; lvec_prtcls];          # one vec with info about length

        PetscBinaryWrite_Vec(fname, lvec_output)            # Write PETSc vector as binary file

    end
end


# Internal routine to retrieve indices of local portion of the grid
function get_ind(x,xc,Nprocx)
    if Nprocx == 1
        xi       = length(x);
        ix_start = [1];
        ix_end   = [length(x)];
    else

        xi = zeros(Int64,Nprocx)
        for k= 1:Nprocx
            if k==1
                xi[k] = length(x[ (x .>=xc[k]) .& (x .<=xc[k+1]) ]);
            else
                xi[k] = length(x[ (x.>xc[k]) .& (x.<=xc[k+1])]);
            end
        end
        ix_start = cumsum( [0; xi[1:end-1]] ) .+ 1;
        ix_end   = cumsum(xi[1:end]);
    end


    return xi,ix_start,ix_end
end

# Internal routine
function get_numscheme(Nprocx,Nprocy,Nprocz)
    n   = zeros(Int64, Nprocx*Nprocy*Nprocz)
    nix = zeros(Int64, Nprocx*Nprocy*Nprocz)
    njy = zeros(Int64, Nprocx*Nprocy*Nprocz)
    nkz = zeros(Int64, Nprocx*Nprocy*Nprocz)
    
    num=0;
    for k=1:Nprocz
        for j=1:Nprocy
            for i=1:Nprocx
                num=num+1;
                n[num]   = num;
                nix[num]= i;
                njy[num]= j;
                nkz[num]= k;
            end
        end
    end
    
    return n,nix,njy,nkz
end


# Internal routine, to write a PETSc vector (as Float64)
"""
    PetscBinaryWrite_Vec(filename, A)

Writes a vector `A` to disk, such that it can be read with `PetscBinaryRead` (which assumes a Big Endian type)

"""
function PetscBinaryWrite_Vec(filename, A)

    # Note: use "hton" to transfer to Big Endian type, which is what PETScBinaryRead expects
    open(filename,"w+") do f
        n               =   length(A);
        nummark         =   A[1];           # number of markers

        write(f,hton(Float64(1211214)));    # header (not actually used)
        write(f,hton(Float64(nummark)));    # info about # of markers written
        
        for i=2:n
             write(f,hton(Float64(A[i])));  # Write data itself    
        end

    end


end

"""
    nProcX,nProcY,nProcZ, xc,yc,zc, nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(filename)

Reads a LaMEM processor partitioning file, used to create marker files, and returns the parallel layout 

"""
function GetProcessorPartitioning(filename)

    io = open(filename, "r")
    
    nProcX = ntoh(read(io,Int32))
    nProcY = ntoh(read(io,Int32))
    nProcZ = ntoh(read(io,Int32))

    nNodeX = ntoh(read(io,Int32))
    nNodeY = ntoh(read(io,Int32))
    nNodeZ = ntoh(read(io,Int32))

    iX = [ntoh(read(io,Int32)) for i=1:nProcX+1];
    iY = [ntoh(read(io,Int32)) for i=1:nProcY+1];
    iZ = [ntoh(read(io,Int32)) for i=1:nProcZ+1];

    CharLength = ntoh(read(io,Float64))
    xcoor = [ntoh(read(io,Float64)) for i=1:nNodeX].*CharLength;
    ycoor = [ntoh(read(io,Float64)) for i=1:nNodeY].*CharLength;
    zcoor = [ntoh(read(io,Float64)) for i=1:nNodeZ].*CharLength;
    
    xc = xcoor[iX .+ 1]
    yc = ycoor[iY .+ 1]
    zc = zcoor[iZ .+ 1]

    close(io)

    return  nProcX,nProcY,nProcZ, 
            xc,yc,zc, 
            nNodeX,nNodeY,nNodeZ
           
end