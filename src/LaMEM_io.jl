using Base: Int64, Float64, NamedTuple
using Printf
using Glob
using Interpolations

# LaMEM I/O
#
# These are routines that help to create a LaMEM marker files from a ParaviewData structure, which can be used to perform geodynamic simulations
# We also include routines with which we can read LaMEM *.pvtr files into julia

export LaMEM_grid, ReadLaMEM_InputFile
export Save_LaMEMMarkersParallel, Save_LaMEMTopography
export GetProcessorPartitioning, ReadData_VTR, ReadData_PVTR, CreatePartitioningFile

"""
Structure that holds information about the LaMEM grid (usually read from an input file).
"""
struct LaMEM_grid <: AbstractGeneralGrid
    # number of markers per element
    nmark_x :: Int64
    nmark_y :: Int64
    nmark_z :: Int64
    # total number of markers
    nump_x  :: Int64
    nump_y  :: Int64
    nump_z  :: Int64
    # total number of elements in grid
    nel_x   :: Int64
    nel_y   :: Int64
    nel_z   :: Int64
    # extent of the grid
    W       ::  Float64
    L       ::  Float64
    H       ::  Float64
    # start and end coordinates of grid segments
    coord_x
    coord_y
    coord_z
    # 1D vectors with marker coordinates
    x_vec
    y_vec
    z_vec
    # grid with marker coordinates
    X
    Y
    Z
    # 1D vectors with node coordinates
    xn_vec
    yn_vec
    zn_vec
    # grid with node coordinates
    Xn
    Yn
    Zn
end

"""
    ParaviewData(Grid::LaMEM_grid, fields::NamedTuple)

Creates a `ParaviewData` struct from a LaMEM grid and from fields stored on that grid. Note that one needs to have a field `Phases` and optionally a field `Temp` to create LaMEM marker files.
"""
ParaviewData(Grid::LaMEM_grid, fields::NamedTuple) = ParaviewData(Grid.X, Grid.Y, Grid.Z, fields)

"""
    CartData(Grid::LaMEM_grid, fields::NamedTuple)

Creates a `CartData` struct from a LaMEM grid and from fields stored on that grid. Note that one needs to have a field `Phases` and optionally a field `Temp` to create LaMEM marker files.
"""
CartData(Grid::LaMEM_grid, fields::NamedTuple) = CartData(Grid.X, Grid.Y, Grid.Z, fields)

"""
    Below = belowSurface(Data_LaMEM::LaMEM_grid, DataSurface_Cart::CartData)

Determines if points within the 3D `LaMEM_grid` structure are below the Cartesian surface DataSurface_Cart
"""
function belowSurface(Grid::LaMEM_grid, DataSurface_Cart::CartData)
    return aboveSurface(CartData(Grid,(Z=Grid.Z,)), DataSurface_Cart; above=false)
end

"""
    Above = aboveSurface(Data_LaMEM::LaMEM_grid, DataSurface_Cart::CartData)

Determines if points within the 3D `LaMEM_grid` structure are above the Cartesian surface DataSurface_Cart
"""
function aboveSurface(Grid::LaMEM_grid, DataSurface_Cart::CartData)
    return aboveSurface(CartData(Grid,(Z=Grid.Z,)), DataSurface_Cart; above=true)
end


"""
    value = ParseValue_LaMEM_InputFile(file,keyword,type; args::String=nothing)

Extracts a certain `keyword` from a LaMEM input `file` and convert it to a certain type.
Optionally, you can also pass command-line arguments which will override the value read from the input file.

# Example 1:
```julia
julia> nmark_z = ParseValue_LaMEM_InputFile("SaltModels.dat","nmark_z",Int64)
```

# Example 2:
```julia
julia> nmark_z = ParseValue_LaMEM_InputFile("SaltModels.dat","nmark_z",Int64, args="-nel_x 128 -coord_x -4,4")
```

"""
function ParseValue_LaMEM_InputFile(file,keyword,type; args::Union{String,Nothing}=nothing)
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

    # parse command-line arguments
    if !isnothing(args)
        value = ParseValue_CommandLineArgs(args, keyword, type, value)
    end

    return value
end


"""
    This parses a LaMEM command line argument string and checks if the keyword exists there
"""
function ParseValue_CommandLineArgs(args,keyword,type, value)
    args_vec = split(args,"-"*keyword)

    if length(args_vec)==2
        # we found the keyword
        args_vec_keyword = split(args_vec[2])
        str = args_vec_keyword[1]               # first block after keyword is what we want
        str_strip = replace(str, "," => " ")    # in case we have an array of values
        value = parse.(type, split(str_strip))  # puts an array of values in a vector

        if length(value)==1
            value=value[1];
        end
    end

    return value
end


"""
    Grid::LaMEM_grid = ReadLaMEM_InputFile(file, args::Union{String,Nothing}=nothing)

Parses a LaMEM input file and stores grid information in the `Grid` structure.
Optionally, you can pass LaMEM command-line arguments as well.

# Example 1
```julia
julia> Grid = ReadLaMEM_InputFile("SaltModels.dat")
LaMEM Grid:
nel         : (32, 32, 32)
marker/cell : (3, 3, 3)
markers     : (96, 96, 96)
x           ϵ [-3.0 : 3.0]
y           ϵ [-2.0 : 2.0]
z           ϵ [-2.0 : 0.0]
```

# Example 2 (with command-line arguments)
```julia
julia> Grid = ReadLaMEM_InputFile("SaltModels.dat", args="-nel_x 64 -coord_x -4,4")
LaMEM Grid:
  nel         : (64, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (192, 96, 96)
  x           ϵ [-4.0 : 4.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
```

"""
function ReadLaMEM_InputFile(file; args::Union{String,Nothing}=nothing )

    # read information from file
    nmark_x   = ParseValue_LaMEM_InputFile(file,"nmark_x",Int64, args=args);
    nmark_y   = ParseValue_LaMEM_InputFile(file,"nmark_y",Int64, args=args);
    nmark_z   = ParseValue_LaMEM_InputFile(file,"nmark_z",Int64, args=args);

    nel_x     = ParseValue_LaMEM_InputFile(file,"nel_x",Int64, args=args);
    nel_y     = ParseValue_LaMEM_InputFile(file,"nel_y",Int64, args=args);
    nel_z     = ParseValue_LaMEM_InputFile(file,"nel_z",Int64, args=args);

    coord_x   = ParseValue_LaMEM_InputFile(file,"coord_x",Float64, args=args);
    coord_y   = ParseValue_LaMEM_InputFile(file,"coord_y",Float64, args=args);
    coord_z   = ParseValue_LaMEM_InputFile(file,"coord_z",Float64, args=args);

    nseg_x   = ParseValue_LaMEM_InputFile(file,"nseg_x",Int64, args=args);
    nseg_y   = ParseValue_LaMEM_InputFile(file,"nseg_y",Int64, args=args);
    nseg_z   = ParseValue_LaMEM_InputFile(file,"nseg_z",Int64, args=args);

    bias_x   = ParseValue_LaMEM_InputFile(file,"bias_x",Float64, args=args);
    bias_y   = ParseValue_LaMEM_InputFile(file,"bias_y",Float64, args=args);
    bias_z   = ParseValue_LaMEM_InputFile(file,"bias_z",Float64, args=args);

    # compute information from file
    W         = coord_x[end]-coord_x[1];
    L         = coord_y[end]-coord_y[1];
    H         = coord_z[end]-coord_z[1];

    nel_x_tot = sum(nel_x);
    nel_y_tot = sum(nel_y);
    nel_z_tot = sum(nel_z);

    nump_x    = nel_x_tot*nmark_x;
    nump_y    = nel_y_tot*nmark_y;
    nump_z    = nel_z_tot*nmark_z;

    # Create 1D coordinate vectors (either regular or refined)
    xn, x = Create1D_grid_vector(coord_x, nel_x, nmark_x, nseg_x, bias_x)
    yn, y = Create1D_grid_vector(coord_y, nel_y, nmark_y, nseg_y, bias_y)
    zn, z = Create1D_grid_vector(coord_z, nel_z, nmark_z, nseg_z, bias_z)

    # node grid
    Xn,Yn,Zn = XYZGrid(xn, yn, zn);

    # marker grid
    X,Y,Z    = XYZGrid(x, y, z);

    # finish Grid
    Grid    =  LaMEM_grid(  nmark_x,    nmark_y,    nmark_z,
    nump_x,     nump_y,     nump_z,
    nel_x_tot,  nel_y_tot,  nel_z_tot,
    W,          L,          H,
    coord_x,    coord_y,    coord_z,
    x,          y,          z,
    X,          Y,          Z,
    xn,         yn,         zn,
    Xn,         Yn,         Zn);

    return Grid
end

"""
Returns 1D coordinate vectors of grid points and of marker locations for a regular spacing
"""
function Create1D_grid_vector(coord::Vector{Float64}, nel::Int64, nmark::Int64, nseg::Union{Nothing, Int64}, bias::Union{Nothing, Float64})
    W  = coord[end] - coord[1]
    Δ  = W / nel;
    xn = range(coord[1], coord[end], length=nel+1);   # coordinates of the normals to the cells

    nump = nmark*nel
    Δ_m = W / nump;
    x  = range(coord[1]+ Δ_m/2, coord[end] - Δ_m/2, length=nump);
    return xn, x
end

"""
Returns 1D coordinate vectors of grid points and of marker locations for a regular spacing
"""
function Create1D_grid_vector(coord::Vector{T}, nel::Vector{I}, nmark::I, nseg::I, bias::Union{Nothing, T, Vector{T}}) where {T<:Float64, I<:Int64}
    if isnothing(bias)
        bias = ones(length(nel))
    end

    xn  = make1DCoords(nseg, nel, coord, bias);
    x   = make1DMarkerCoords(xn, nmark);

    return xn, x
end

function make1DMarkerCoords(xn::Array{Float64, 1}, nmark::Int64)
    # preallocate
    nel  = length(xn) - 1
    nump = nel * nmark;
    x    = zeros(Float64, nump);

    # compute coordinates
    for i = 1 : nel
        # start of cell
        x0 = xn[i];
        # markers spacing inside cell
        dx = (xn[i+1] - x0) / nmark;

        # compute position
        for j = 1 : nmark
            x[nmark*i-(nmark-j)] = x0 + dx/2 + (j-1)*dx;
        end
    end

    return x
end

function make1DCoords(nseg::Int64, nel, coord::Array{Float64, 1}, bias)
    # preallocate
    nel_tot = sum(nel);
    x       = zeros(Float64, nel_tot+1);

    for i = 1 : nseg
        # indices of this segment in the coordinate vector
        if i == 1
            indE = nel[1] + 1
        else
            indE = sum(nel[1:i]) + 1;
        end
        indS = indE - nel[i];

        # compute coordinates
        x[indS:indE] = makeCoordSegment(coord[i], coord[i+1], nel[i], bias[i]);
    end

    return x
end

function makeCoordSegment(xStart::Float64, xEnd::Float64, numCells::Int64, bias::Float64)
    # average cell size
    avgSize = (xEnd - xStart) / numCells;

    # uniform case
    if bias == 1.0
        x = Array(xStart : avgSize : xEnd);
    # non-uniform case
    else
        x = zeros(Float64, numCells+1)
        # cell size limits
        begSize = 2.0 * avgSize / (1.0 + bias);
        endSize = bias * begSize;

        # cell size increment (negative for bias < 1)
        dx      = (endSize - begSize) / (numCells - 1);

        # generate coordinates
        x[1]    = xStart;
        for i = 2 : numCells + 1
            x[i] = x[i-1] + begSize + (i-2)*dx;
        end

        # overwrite last coordinate
        x[end] = xEnd;
    end

    return x
end

# Print an overview of the LaMEM Grid struct:
function Base.show(io::IO, d::LaMEM_grid)
    println(io,"LaMEM Grid: ")
    println(io,"  nel         : ($(d.nel_x), $(d.nel_y), $(d.nel_z))")
    println(io,"  marker/cell : ($(d.nmark_x), $(d.nmark_y), $(d.nmark_z))")
    println(io,"  markers     : ($(d.nump_x), $(d.nump_y), $(d.nump_z))")
    println(io,"  x           ϵ [$(d.coord_x[1]) : $(d.coord_x[end])]")
    println(io,"  y           ϵ [$(d.coord_y[1]) : $(d.coord_y[end])]")
    println(io,"  z           ϵ [$(d.coord_z[1]) : $(d.coord_z[end])]")
end

"""
    Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers", verbose=true, is64bit=false)

Saves a LaMEM marker file from the `CartData` structure `Grid`. It must have a field called `Phases`, holding phase information (as integers) and optionally a field `Temp` with temperature info.
It is possible to provide a LaMEM partitioning file `PartitioningFile`. If not, output is assumed to be for one processor. By default it is assumed that the partitioning file was generated on a 32bit PETSc installation. If `Int64` was used instead, set the flag.

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
function Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers", verbose=true, is64bit=false)

    x = ustrip.(Grid.x.val[:,1,1]);
    y = ustrip.(Grid.y.val[1,:,1]);
    z = ustrip.(Grid.z.val[1,1,:]);

    if haskey(Grid.fields,:Phases)
        Phases = Grid.fields[:Phases];
    else
        error("You must provide the field :Phases in the structure")
    end

    if haskey(Grid.fields,:Temp)
        Temp = Grid.fields[:Temp];
    else
        if verbose
            println("Field :Temp is not provided; setting it to zero")
        end
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
        nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(PartitioningFile, is64bit=is64bit)
        if verbose
            @show  Nprocx,Nprocy,Nprocz, xc,yc,zc, nNodeX,nNodeY,nNodeZ
        end
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

        part_x   = ustrip.(Grid.x.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]]);
        part_y   = ustrip.(Grid.y.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]]);
        part_z   = ustrip.(Grid.z.val[x_start[n]:x_end[n],y_start[n]:y_end[n],z_start[n]:z_end[n]]);
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
        if verbose
            println("Writing LaMEM marker file -> $fname")                   # print info
        end
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
    nProcX,nProcY,nProcZ, xc,yc,zc, nNodeX,nNodeY,nNodeZ = GetProcessorPartitioning(filename; is64bit=false)

Reads a LaMEM processor partitioning file, used to create marker files, and returns the parallel layout.
By default this is done for a 32bit PETSc installation, which will fail if you actually use a 64bit version.

"""
function GetProcessorPartitioning(filename; is64bit=false)

    if is64bit
        typ=Int64
    else
        typ=Int32
    end
    io = open(filename, "r")


    nProcX = ntoh(read(io,typ))
    nProcY = ntoh(read(io,typ))
    nProcZ = ntoh(read(io,typ))

    nNodeX = ntoh(read(io,typ))
    nNodeY = ntoh(read(io,typ))
    nNodeZ = ntoh(read(io,typ))

    iX = [ntoh(read(io,typ)) for i=1:nProcX+1];
    iY = [ntoh(read(io,typ)) for i=1:nProcY+1];
    iZ = [ntoh(read(io,typ)) for i=1:nProcZ+1];

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




"""
    coord, Data_3D_Arrays, Name_Vec = ReadData_VTR(fname)

Reads a VTR (structured grid) VTK file `fname` and extracts the coordinates, data arrays and names of the data.
In general, this only contains a piece of the data, and one should open a `*.pvtr` file to retrieve the full data
"""
function ReadData_VTR(fname, FullSize)
    file = open(fname, "r")

    header = true
    num = 1;
    CoordOffset = zeros(Int64,3);
    Offset_Vec  =   [];
    Name_Vec    =   [];     Type_Vec = [];
    NumComp_Vec =   [];     PieceExtent=[]; WholeExtent=[];
    while header==true

        line        = readline(file)
        line_strip  = lstrip(line)
        if startswith(line_strip, "<RectilinearGrid WholeExtent")
            id_start    = findfirst("\"", line_strip)[1]+1
            id_end      = findlast("\"", line_strip)[1]-1
            WholeExtent = parse.(Int64,split(line_strip[id_start:id_end]))
        end
        if startswith(line_strip, "<Piece Extent=")
            id_start    = findfirst("\"", line_strip)[1]+1
            id_end      =  findlast("\"", line_strip)[1]-1
            PieceExtent =  parse.(Int64,split(line_strip[id_start:id_end]))


        end
        if startswith(line_strip, "<Coordinates>")
            # Read info where the coordinates are stored
            Type, Name, NumberOfComponents, CoordOffset[1]  = Parse_VTR_Line(readline(file)); num += 1
            Type, Name, NumberOfComponents, CoordOffset[2]  = Parse_VTR_Line(readline(file)); num += 1
            Type, Name, NumberOfComponents, CoordOffset[3]  = Parse_VTR_Line(readline(file)); num += 1
        end

        if startswith(line_strip, "<PointData>")
            line_strip  = lstrip(readline(file))
            while ~startswith(line_strip, "</PointData>")
                Type, Name, NumberOfComponents, Offset  = Parse_VTR_Line(line_strip);  num += 1

                Offset_Vec  = [Offset_Vec;  Offset];
                Name_Vec    = [Name_Vec;    Name];
                Type_Vec    = [Type_Vec;   Type];
                NumComp_Vec = [NumComp_Vec; NumberOfComponents];
                line_strip  = lstrip(readline(file))
            end
        end

        if startswith(line_strip, "<CellData>")
            line_strip  = lstrip(readline(file))
            while ~startswith(line_strip, "</CellData>")
                Type, Name, NumberOfComponents, Offset  = Parse_VTR_Line(line_strip);  num += 1

                Offset_Vec  = [Offset_Vec;  Offset];
                Name_Vec    = [Name_Vec;    Name];
                Type_Vec    = [Type_Vec;   Type];
                NumComp_Vec = [NumComp_Vec; NumberOfComponents];
                line_strip  = lstrip(readline(file))
                 # if we have cell Data, for some reason we need to increment this by one.
                PieceExtent[1:2:end] .+= 1
            end

        end

        if startswith(line_strip, "<AppendedData ")
            header=false
        end

        num += 1
    end

    # Skip to beginning of raw data (linebreak)
    skip(file, 5)
    start_bin = position(file);     # start of binary data

    # Determine the end of the raw data
    seekend(file);
    skip(file, -29)

    end_bin = position(file);

    # Start with reading the coordinate arrays:
    coord_x     =   ReadBinaryData(file, start_bin, CoordOffset[1],   (PieceExtent[2]-PieceExtent[1]+1)*sizeof(Float32))
    coord_y     =   ReadBinaryData(file, start_bin, CoordOffset[2],   (PieceExtent[4]-PieceExtent[3]+1)*sizeof(Float32))
    coord_z     =   ReadBinaryData(file, start_bin, CoordOffset[3],   (PieceExtent[6]-PieceExtent[5]+1)*sizeof(Float32))


    # Read data arrays:
    Data_3D_Arrays  = [];
    ix = PieceExtent[1]:PieceExtent[2];
    iy = PieceExtent[3]:PieceExtent[4];
    iz = PieceExtent[5]:PieceExtent[6];
    numPoints       = length(ix)*length(iy)*length(iz);

    coord_x_full = zeros(Float64, FullSize[1]);
    coord_y_full = zeros(Float64, FullSize[2]);
    coord_z_full = zeros(Float64, FullSize[3]);

    coord_x_full[ix] = coord_x[1:length(ix)];
    coord_y_full[iy] = coord_y[1:length(iy)];
    coord_z_full[iz] = coord_z[1:length(iz)];

    for i=1:length(Name_Vec)-1

        data3D      =   ReadBinaryData(file, start_bin, Offset_Vec[i],    numPoints*NumComp_Vec[i]*sizeof(Float32) )
        data3D      =   getArray(data3D, PieceExtent, NumComp_Vec[i]);

        data3D_full =   zeros(Float64,NumComp_Vec[i],FullSize[1],FullSize[2],FullSize[3])  # Generate full data

        # ugly hack to make it work with parallel files
        ix_left = ix; ix_right = 1:length(ix_left);
        iy_left = iy; iy_right = 1:length(iy_left);
        iz_left = iz; iz_right = 1:length(iz_left);
        if ix_left[1]>1; ix_left = ix_left[2:end]; ix_right=ix_right[2:end];    end
        if iy_left[1]>1; iy_left = iy_left[2:end]; iy_right=iy_right[2:end];    end
        if iz_left[1]>1; iz_left = iz_left[2:end]; iz_right=iz_right[2:end];    end

        data3D_full[1:NumComp_Vec[i], ix_left, iy_left, iz_left]        = data3D[1:NumComp_Vec[i],ix_right, iy_right, iz_right];
        #data3D_full[1:NumComp_Vec[i], ix, iy, iz]        = data3D;

        Data_3D_Arrays = [Data_3D_Arrays; data3D_full]
    end
    i=length(Name_Vec);

    if Type_Vec[i]=="UInt8"
        data3D   =   ReadBinaryData(file, start_bin, Offset_Vec[i],    numPoints*NumComp_Vec[i]*sizeof(UInt8), DataType=UInt8)
    else
        data3D   =   ReadBinaryData(file, start_bin, Offset_Vec[i],    numPoints*NumComp_Vec[i]*sizeof(Float32) )
    end

    data3D   =   getArray(data3D, PieceExtent, NumComp_Vec[i]);
    data3D_full =   zeros(Float64,NumComp_Vec[i],FullSize[1],FullSize[2],FullSize[3])  # Generate full d
    data3D_full[1:NumComp_Vec[i], ix, iy, iz]        = data3D[1:NumComp_Vec[i],1:length(ix),1:length(iy),1:length(iz)];

    Data_3D_Arrays = [Data_3D_Arrays; data3D_full]

    return coord_x_full, coord_y_full, coord_z_full, Data_3D_Arrays, Name_Vec, NumComp_Vec, ix, iy, iz
end

 # Parses a line of a *.vtr file & retrieve Type/Name/NumberOfComponents/Offset
 function Parse_VTR_Line(line)
    line_strip  = lstrip(line)

    # Retrieve Type
    if findfirst("type", line_strip) != nothing
        id_start    = findfirst("type", line_strip)[1]+6

        line_strip  = line_strip[id_start:end]
        id_end      = findfirst("\"", line_strip)[1]-1
        Type        = line_strip[1:id_end]
        line_strip  = line_strip[id_end:end]
    else
        Type=nothing;
    end

    # Retrieve Name
    if findfirst("Name", line_strip) != nothing
        id_start    = findfirst("Name", line_strip)[1]+6
        line_strip  = line_strip[id_start:end]
        id_end      = findfirst("\"", line_strip)[1]-1
        Name        = line_strip[1:id_end]
        line_strip  = line_strip[id_end:end]
    else
        Name=nothing
    end

    # Retrieve number of components
    if findfirst("NumberOfComponents", line_strip) != nothing
        id_start    = findfirst("NumberOfComponents", line_strip)[1]+20
        line_strip  = line_strip[id_start:end]
        id_end      = findfirst("\"", line_strip)[1]-1
        NumberOfComponents     = parse(Int64,line_strip[1:id_end])
        line_strip  = line_strip[id_end:end]
    else
        NumberOfComponents=nothing
    end

    # Offset
    if findfirst("offset", line_strip) != nothing
        id_start    = findfirst("offset", line_strip)[1]+8
        line_strip  = line_strip[id_start:end]
        id_end      = findfirst("\"", line_strip)[1]-1
        Offset     = parse(Int64,line_strip[1:id_end])
    else
        Offset=nothing;
    end
    return Type, Name, NumberOfComponents, Offset
end


function getArray(data, PieceExtent, NumComp)
    data        =   reshape(data, (NumComp, PieceExtent[2]-PieceExtent[1]+1,  PieceExtent[4]-PieceExtent[3]+1,  PieceExtent[6]-PieceExtent[5]+1))
    return data
end

function ReadBinaryData(file::IOStream, start_bin::Int64, Offset::Int64, BytesToRead; DataType=Float32)

    seekstart(file);                                # go to start
    skip(file, start_bin+Offset)                    # move to beginning of raw binary data
    buffer      =   read(file,BytesToRead)          # Read necesaary bytes
    data        =   reinterpret(DataType,buffer)    # Transfer to buffer

    data        =   Float64.(data[1:end]);        # Transfer to Float64
    return data
end

"""
    Data::ParaviewData = ReadData_PVTR(fname, dir)

Reads a parallel, rectilinear, `*.vts` file with the name `fname` and located in `dir` and create a 3D `Data` struct from it.

# Example
```julia
julia> Data = ReadData_PVTR("Haaksbergen.pvtr", "./Timestep_00000005_3.35780500e-01/")
ParaviewData
  size  : (33, 33, 33)
  x     ϵ [ -3.0 : 3.0]
  y     ϵ [ -2.0 : 2.0]
  z     ϵ [ -2.0 : 0.0]
  fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
```
"""
function  ReadData_PVTR(fname, dir)
    file = open(joinpath(dir,fname), "r")

    header = true
    num = 1;
    FullSize= (1,1,1);
    num_data_sets = 1;
    Data_3D=[]; coord_x=[]; coord_y=[]; coord_z=[]; NumComp=[]; Names=[]
    while header==true

        line        = readline(file)
        line_strip  = lstrip(line)
        if startswith(line_strip, "<PRectilinearGrid")
            id_start    = findfirst("WholeExtent=", line_strip)[1]+13
            line_strip  = line_strip[id_start:end]
            id_end      = findfirst("\"", line_strip)[1]-1
            line_piece  = line_strip[1:id_end]

            WholeExtent = parse.(Int64,split(line_piece))
            FullSize    = (WholeExtent[2],WholeExtent[4],WholeExtent[6])
        end


        if startswith(line_strip, "<Piece")
            id_start    = findfirst("Source=", line_strip)[1]+8
            line_strip  = line_strip[id_start:end]
            id_end      = findfirst("\"", line_strip)[1]-1
            fname_piece = line_strip[1:id_end]

            if num_data_sets==1
                coord_x, coord_y, coord_z, Data_3D, Names, NumComp, ix,iy,iz = ReadData_VTR(joinpath(dir,fname_piece), FullSize);
            else
                coord_x1, coord_y1, coord_z1, Data_3D1, Names, NumComp, ix,iy,iz  = ReadData_VTR(joinpath(dir,fname_piece), FullSize);
                coord_x[ix]   = coord_x1[ix];
                coord_y[iy]   = coord_y1[iy];
                coord_z[iz]   = coord_z1[iz];

                Data_3D = Data_3D+Data_3D1;

            end
            num_data_sets += 1
        end


        if startswith(line_strip, "</PRectilinearGrid")
            header=false;
        end
    end

    # Create a named-Tuple out of the fields
    NamesSymbol = [];
    for i=1:length(Names)
        id = findfirst(" ", Names[i])
        if id == nothing
            Names_Strip = Names[i]
        else
            Names_Strip = Names[i][1:findfirst(" ", Names[i])[1]-1];
        end
        NamesSymbol = [NamesSymbol; Names_Strip]
    end

  #  NamesSymbol =   [Names[i][1:findfirst(" ", Names[i])[1]-1] for i=1:length(Names)]
    Names1      =   Symbol.(NamesSymbol)

    Data_Array = [];
    num     =   1;
    for i=1:length(NumComp)
        data        =   Data_3D[num:num+NumComp[i]-1,:,:,:];
        data_arrays =   [data[i,:,:,:] for i=1:size(data,1)]
        data_tuple  =   tuple(data_arrays...)

        if size(data,1)>1
            Data_NamedTuple = NamedTuple{(Names1[i],)}((data_tuple,))
        else
            Data_NamedTuple = NamedTuple{(Names1[i],)}((data_tuple[1],))
        end
        Data_Array = [Data_Array; Data_NamedTuple]

        num = num+NumComp[i];
    end

    # Merge vector with tuples into a NamedTuple
    fields = Data_Array[1];
    for i=2:length(Data_Array)
        fields = merge(fields, Data_Array[i])
    end

    # Create a ParaviewData struct from it.
    X,Y,Z       =   XYZGrid(coord_x, coord_y, coord_z)
    DataC       =   ParaviewData(X,Y,Z, fields);

    return DataC
end

"""
    Save_LaMEMTopography(Topo::CartData, filename::String)

This writes a topography file `Topo` for use in LaMEM, which should have size `(nx,ny,1)` and contain the field `:Topography`
"""
function Save_LaMEMTopography(Topo::CartData, filename::String)

    if (size(Topo.z.val,3) != 1)
        error("Not a valid `CartData' Topography file (size in 3rd dimension should be 1)")
    end
    if !haskey(Topo.fields,:Topography)
        error("The topography `CartData` structure requires a field :Topography")
    end

    # get grid properties
    nx = Float64(size(Topo.fields.Topography,1));
    ny = Float64(size(Topo.fields.Topography,2));
    x0 = ustrip(Topo.x.val[1,1,1]);
    y0 = ustrip(Topo.y.val[1,1,1]);

    # LaMEM wants a uniform grid, so interpolate if necessary
    if length(unique(trunc.(diff(Topo.x.val[:,1,1]), digits=8))) > 1 || length(unique(trunc.(diff(Topo.y.val[1,:,1]), digits=8))) > 1
        x1       = ustrip(Topo.x.val[end,1,1]);
        y1       = ustrip(Topo.y.val[1,end,1]);
        dx       = (x1-x0) / (nx-1);
        dy       = (y1-y0) / (ny-1);

        itp      = LinearInterpolation((Topo.x.val[:,1,1], Topo.y.val[1,:,1]), ustrip.(Topo.fields.Topography[:,:,1]));
        Topo_itp = [itp(x,y) for x in x0:dx:x1, y in y0:dy:y1];

        # Code the topograhic data into a vector
        Topo_vec = [ nx;ny;x0;y0;dx;dy; Topo_itp[:]]
    else
        dx = ustrip(Topo.x.val[2,2,1]) - x0
        dy = ustrip(Topo.y.val[2,2,1]) - y0
        # Code the topograhic data into a vector
        Topo_vec = [ nx;ny;x0;y0;dx;dy; ustrip.(Topo.fields.Topography[:])]
    end

    # Write as PetscBinary file
    PetscBinaryWrite_Vec(filename, Topo_vec)

    println("Written LaMEM topography file: $(filename)")

    return nothing
end

"""
    CreatePartitioningFile(LaMEM_input::String, NumProc::Int64; LaMEM_dir::String=pwd(), LaMEM_options::String="", MPI_dir="", verbose=true)

This executes LaMEM for the input file `LaMEM_input` & creates a parallel partitioning file for `NumProc` processors.
The directory where the LaMEM binary is can be specified; if not it is assumed to be in the current directory.
Likewise for the `mpiexec` directory (if not specified it is assumed to be available on the command line).

"""
function CreatePartitioningFile(LaMEM_input::String,NumProc::Int64; LaMEM_dir::String=pwd(), LaMEM_options="", MPI_dir="", verbose=true)

    # Create string to execute LaMEM
    mpi_str     =  MPI_dir*"mpiexec -n $(NumProc) "
    LaMEM_str   =  LaMEM_dir*"/"*"LaMEM -ParamFile "*LaMEM_input*" -mode save_grid "
    str         =  mpi_str*LaMEM_str

    if verbose==true
        println("Executing command: $str")
    end
    # Run
    exit=run(`sh -c $str`, wait=false);

    # Retrieve newest file
    if success(exit)
        files=readdir(glob"ProcessorPartitioning_*.bin")
        time_modified = zeros(length(files))
        for (i,file) in enumerate(files)
            time_modified[i] = stat(file).mtime
        end
        id          = findall(time_modified.==maximum(time_modified))   # last modified
        PartFile    = files[id]
        if verbose==true
            println("Successfully generated PartitioningFile: $(PartFile[1])")
        end
    else
        error("Something went wrong with executing command ")
    end

    return PartFile[1]

end

"""
    X,Y,Z = coordinate_grids(Data::LaMEM_grid)

Returns 3D coordinate arrays
"""
function coordinate_grids(Data::LaMEM_grid)

    return Data.X, Data.Y, Data.Z
end
