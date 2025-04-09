using Base: Int64, Float64, NamedTuple
using Printf
using Glob
using Interpolations

import Base: show, size

# LaMEM I/O
#
# These are routines that help to create a LaMEM marker files from a ParaviewData structure, which can be used to perform geodynamic simulations
# We also include routines with which we can read LaMEM *.pvtr files into julia

export LaMEM_grid, read_LaMEM_inputfile
export save_LaMEM_markers_parallel, save_LaMEM_topography
export get_processor_partitioning, read_data_VTR, read_data_PVTR, create_partitioning_file
export crop_bounds, get_proc_bound, get_proc_grid, get_particles_distribution, get_LaMEM_grid_info, get_processor_partitioning_info, check_markers_directory, setup_model_domain, LaMEM_partitioning_info

"""
Structure that holds information about the LaMEM partitioning
"""
struct LaMEM_partitioning_info <: AbstractGeneralGrid

    # Number of processors in each direction
    nProcX::Int64
    nProcY::Int64
    nProcZ::Int64
    # Number of nodes in each direction
    nNodeX::Int64
    nNodeY::Int64
    nNodeZ::Int64
    # Coordinates of the nodes end of each processor
    xc
    yc
    zc

end

"""
Structure that holds information about the LaMEM particles distribution for partitioning
"""
struct particles_distribution <: AbstractGeneralGrid

    x_start
    x_end
    y_start
    y_end
    z_start
    z_end

end

"""
Structure that holds information about the LaMEM grid (usually read from an input file).
"""
struct LaMEM_grid <: AbstractGeneralGrid
    # number of markers per element
    nmark_x::Int64
    nmark_y::Int64
    nmark_z::Int64
    # total number of markers
    nump_x::Int64
    nump_y::Int64
    nump_z::Int64
    # total number of elements in grid
    nel_x::Int64
    nel_y::Int64
    nel_z::Int64
    # extent of the grid
    W::Float64
    L::Float64
    H::Float64
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
size(d::LaMEM_grid) = (d.nump_x, d.nump_y, d.nump_z)

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
    Below = below_surface(Data_LaMEM::LaMEM_grid, DataSurface_Cart::CartData)

Determines if points within the 3D `LaMEM_grid` structure are below the Cartesian surface DataSurface_Cart
"""
function below_surface(Grid::LaMEM_grid, DataSurface_Cart::CartData)
    return above_surface(CartData(Grid, (Z = Grid.Z,)), DataSurface_Cart; above = false)
end

"""
    Above = above_surface(Data_LaMEM::LaMEM_grid, DataSurface_Cart::CartData)

Determines if points within the 3D `LaMEM_grid` structure are above the Cartesian surface DataSurface_Cart
"""
function above_surface(Grid::LaMEM_grid, DataSurface_Cart::CartData)
    return above_surface(CartData(Grid, (Z = Grid.Z,)), DataSurface_Cart; above = true)
end


"""
    value = ParseValue_LaMEM_InputFile(file,keyword,type; args::String=nothing)

Extracts a certain `keyword` from a LaMEM input `file` and convert it to a certain type.
Optionally, you can also pass command-line arguments which will override the value read from the input file.

# Example 1:
```julia-repl
julia> nmark_z = ParseValue_LaMEM_InputFile("SaltModels.dat","nmark_z",Int64)
```

# Example 2:
```julia-repl
julia> nmark_z = ParseValue_LaMEM_InputFile("SaltModels.dat","nmark_z",Int64, args="-nel_x 128 -coord_x -4,4")
```

"""
function ParseValue_LaMEM_InputFile(file, keyword, type; args::Union{String, Nothing} = nothing)
    value = nothing
    for line in eachline(file)
        line_strip = lstrip(line)       # strip leading tabs/spaces

        # Strip comments
        ind = findfirst("#", line)
        if isnothing(ind)
            # no comments
        else
            line_strip = line_strip[1:(ind[1] - 2)]
        end
        line_strip = rstrip(line_strip)       # strip last tabs/spaces

        if startswith(line_strip, keyword)
            ind = findfirst("=", line_strip)
            if type == String
                value = split(line_strip)[3:end]
            else
                value = parse.(type, split(line_strip)[3:end])

                if length(value) == 1
                    value = value[1]
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
function ParseValue_CommandLineArgs(args, keyword, type, value)
    args_vec = split(args, "-" * keyword)

    if length(args_vec) == 2
        # we found the keyword
        args_vec_keyword = split(args_vec[2])
        str = args_vec_keyword[1]               # first block after keyword is what we want
        str_strip = replace(str, "," => " ")    # in case we have an array of values
        value = parse.(type, split(str_strip))  # puts an array of values in a vector

        if length(value) == 1
            value = value[1]
        end
    end

    return value
end


"""
    Grid::LaMEM_grid = read_LaMEM_inputfile(file, args::Union{String,Nothing}=nothing)

Parses a LaMEM input file and stores grid information in the `Grid` structure.
Optionally, you can pass LaMEM command-line arguments as well.

# Example 1
```julia-repl
julia> Grid = read_LaMEM_inputfile("SaltModels.dat")
LaMEM Grid:
nel         : (32, 32, 32)
marker/cell : (3, 3, 3)
markers     : (96, 96, 96)
x           ϵ [-3.0 : 3.0]
y           ϵ [-2.0 : 2.0]
z           ϵ [-2.0 : 0.0]
```

# Example 2 (with command-line arguments)
```julia-repl
julia> Grid = read_LaMEM_inputfile("SaltModels.dat", args="-nel_x 64 -coord_x -4,4")
LaMEM Grid:
  nel         : (64, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (192, 96, 96)
  x           ϵ [-4.0 : 4.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
```

"""
function read_LaMEM_inputfile(file; args::Union{String, Nothing} = nothing)

    # read information from file
    nmark_x = ParseValue_LaMEM_InputFile(file, "nmark_x", Int64, args = args)
    nmark_y = ParseValue_LaMEM_InputFile(file, "nmark_y", Int64, args = args)
    nmark_z = ParseValue_LaMEM_InputFile(file, "nmark_z", Int64, args = args)

    nel_x = ParseValue_LaMEM_InputFile(file, "nel_x", Int64, args = args)
    nel_y = ParseValue_LaMEM_InputFile(file, "nel_y", Int64, args = args)
    nel_z = ParseValue_LaMEM_InputFile(file, "nel_z", Int64, args = args)

    coord_x = ParseValue_LaMEM_InputFile(file, "coord_x", Float64, args = args)
    coord_y = ParseValue_LaMEM_InputFile(file, "coord_y", Float64, args = args)
    coord_z = ParseValue_LaMEM_InputFile(file, "coord_z", Float64, args = args)

    nseg_x = ParseValue_LaMEM_InputFile(file, "nseg_x", Int64, args = args)
    nseg_y = ParseValue_LaMEM_InputFile(file, "nseg_y", Int64, args = args)
    nseg_z = ParseValue_LaMEM_InputFile(file, "nseg_z", Int64, args = args)

    bias_x = ParseValue_LaMEM_InputFile(file, "bias_x", Float64, args = args)
    bias_y = ParseValue_LaMEM_InputFile(file, "bias_y", Float64, args = args)
    bias_z = ParseValue_LaMEM_InputFile(file, "bias_z", Float64, args = args)

    # compute information from file
    W = coord_x[end] - coord_x[1]
    L = coord_y[end] - coord_y[1]
    H = coord_z[end] - coord_z[1]

    nel_x_tot = sum(nel_x)
    nel_y_tot = sum(nel_y)
    nel_z_tot = sum(nel_z)

    nump_x = nel_x_tot * nmark_x
    nump_y = nel_y_tot * nmark_y
    nump_z = nel_z_tot * nmark_z

    # Create 1D coordinate vectors (either regular or refined)
    xn, x = Create1D_grid_vector(coord_x, nel_x, nmark_x, nseg_x, bias_x)
    yn, y = Create1D_grid_vector(coord_y, nel_y, nmark_y, nseg_y, bias_y)
    zn, z = Create1D_grid_vector(coord_z, nel_z, nmark_z, nseg_z, bias_z)

    # node grid
    Xn, Yn, Zn = xyz_grid(xn, yn, zn)

    # marker grid
    X, Y, Z = xyz_grid(x, y, z)

    # finish Grid
    Grid = LaMEM_grid(
        nmark_x, nmark_y, nmark_z,
        nump_x, nump_y, nump_z,
        nel_x_tot, nel_y_tot, nel_z_tot,
        W, L, H,
        coord_x, coord_y, coord_z,
        x, y, z,
        X, Y, Z,
        xn, yn, zn,
        Xn, Yn, Zn
    )

    return Grid
end

"""
Returns 1D coordinate vectors of grid points and of marker locations for a regular spacing
"""
function Create1D_grid_vector(coord::Vector{Float64}, nel::Int64, nmark::Int64, nseg::Union{Nothing, Int64}, bias::Union{Nothing, Float64})
    W = coord[end] - coord[1]
    Δ = W / nel
    xn = range(coord[1], coord[end], length = nel + 1)    # coordinates of the normals to the cells

    nump = nmark * nel
    Δ_m = W / nump
    x = range(coord[1] + Δ_m / 2, coord[end] - Δ_m / 2, length = nump)
    return xn, x
end

"""
Returns 1D coordinate vectors of grid points and of marker locations for a regular spacing
"""
function Create1D_grid_vector(coord::Vector{T}, nel::Vector{I}, nmark::I, nseg::I, bias::Union{Nothing, T, Vector{T}}) where {T <: Float64, I <: Int64}
    if isnothing(bias)
        bias = ones(length(nel))
    end

    xn = make1DCoords(nseg, nel, coord, bias)
    x = make1DMarkerCoords(xn, nmark)

    return xn, x
end

function make1DMarkerCoords(xn::Array{Float64, 1}, nmark::Int64)
    # preallocate
    nel = length(xn) - 1
    nump = nel * nmark
    x = zeros(Float64, nump)

    # compute coordinates
    for i in 1:nel
        # start of cell
        x0 = xn[i]
        # markers spacing inside cell
        dx = (xn[i + 1] - x0) / nmark

        # compute position
        for j in 1:nmark
            x[nmark * i - (nmark - j)] = x0 + dx / 2 + (j - 1) * dx
        end
    end

    return x
end

function make1DCoords(nseg::Int64, nel, coord::Array{Float64, 1}, bias)
    # preallocate
    nel_tot = sum(nel)
    x = zeros(Float64, nel_tot + 1)

    for i in 1:nseg
        # indices of this segment in the coordinate vector
        if i == 1
            indE = nel[1] + 1
        else
            indE = sum(nel[1:i]) + 1
        end
        indS = indE - nel[i]

        # compute coordinates
        x[indS:indE] = makeCoordSegment(coord[i], coord[i + 1], nel[i], bias[i])
    end

    return x
end

function makeCoordSegment(xStart::Float64, xEnd::Float64, numCells::Int64, bias::Float64)
    # average cell size
    avgSize = (xEnd - xStart) / numCells

    # uniform case
    if bias == 1.0
        x = Array(xStart:avgSize:xEnd)
        # non-uniform case
    else
        x = zeros(Float64, numCells + 1)
        # cell size limits
        begSize = 2.0 * avgSize / (1.0 + bias)
        endSize = bias * begSize

        # cell size increment (negative for bias < 1)
        dx = (endSize - begSize) / (numCells - 1)

        # generate coordinates
        x[1] = xStart
        for i in 2:(numCells + 1)
            x[i] = x[i - 1] + begSize + (i - 2) * dx
        end

        # overwrite last coordinate
        x[end] = xEnd
    end

    return x
end

# Print an overview of the LaMEM Grid struct:
function Base.show(io::IO, d::LaMEM_grid)
    println(io, "LaMEM Grid: ")
    println(io, "  nel         : ($(d.nel_x), $(d.nel_y), $(d.nel_z))")
    println(io, "  marker/cell : ($(d.nmark_x), $(d.nmark_y), $(d.nmark_z))")
    println(io, "  markers     : ($(d.nump_x), $(d.nump_y), $(d.nump_z))")
    println(io, "  x           ϵ [$(d.coord_x[1]) : $(d.coord_x[end])]")
    println(io, "  y           ϵ [$(d.coord_y[1]) : $(d.coord_y[end])]")
    return println(io, "  z           ϵ [$(d.coord_z[1]) : $(d.coord_z[end])]")
end

"""
    save_LaMEM_markers_parallel(Grid::CartData; PartitioningFile=empty, directory="./markers", verbose=true, is64bit=false)

Saves a LaMEM marker file from the `CartData` structure `Grid`. It must have a field called `Phases`, holding phase information (as integers) and optionally a field `Temp` with temperature info.
It is possible to provide a LaMEM partitioning file `PartitioningFile`. If not, output is assumed to be for one processor. By default it is assumed that the partitioning file was generated on a 32bit PETSc installation. If `Int64` was used instead, set the flag.

The size of `Grid` should be consistent with what is provided in the LaMEM input file. In practice, the size of the mesh can be retrieved from a LaMEM input file using `read_LaMEM_inputfile`.

# Example

```
julia> Grid    = read_LaMEM_inputfile("LaMEM_input_file.dat")
julia> Phases  = zeros(Int32,size(Grid.X));
julia> Temp    = ones(Float64,size(Grid.X));
julia> Model3D = CartData(Grid, (Phases=Phases,Temp=Temp))
julia> save_LaMEM_markers_parallel(Model3D)
Writing LaMEM marker file -> ./markers/mdb.00000000.dat
```
If you want to create a LaMEM input file for multiple processors:
```
julia> save_LaMEM_markers_parallel(Model3D, PartitioningFile="ProcessorPartitioning_4cpu_1.2.2.bin")
Writing LaMEM marker file -> ./markers/mdb.00000000.dat
Writing LaMEM marker file -> ./markers/mdb.00000001.dat
Writing LaMEM marker file -> ./markers/mdb.00000002.dat
Writing LaMEM marker file -> ./markers/mdb.00000003.dat
```

"""
function save_LaMEM_markers_parallel(Grid::CartData; PartitioningFile = empty, directory = "./markers", verbose = true, is64bit = false)

    x = ustrip.(Grid.x.val[:, 1, 1])
    y = ustrip.(Grid.y.val[1, :, 1])
    z = ustrip.(Grid.z.val[1, 1, :])

    if haskey(Grid.fields, :Phases)
        Phases = Grid.fields[:Phases]
    else
        error("You must provide the field :Phases in the structure")
    end

    if haskey(Grid.fields, :Temp)
        Temp = Grid.fields[:Temp]
    else
        if verbose
            println("Field :Temp is not provided; setting it to zero")
        end
        Temp = zeros(size(Phases))
    end

    if PartitioningFile == empty
        # in case we run this on 1 processor only
        Nprocx = 1
        Nprocy = 1
        Nprocz = 1
        xc, yc, zc = x, y, z
    else
        Nprocx, Nprocy, Nprocz,
            xc, yc, zc,
            nNodeX, nNodeY, nNodeZ = get_processor_partitioning(PartitioningFile, is64bit = is64bit)
        if verbose
            @show  Nprocx, Nprocy, Nprocz, xc, yc, zc, nNodeX, nNodeY, nNodeZ
        end
    end

    Nproc = Nprocx * Nprocy * Nprocz
    num, num_i, num_j, num_k = get_numscheme(Nprocx, Nprocy, Nprocz)

    xi, ix_start, ix_end = get_ind(x, xc, Nprocx)
    yi, iy_start, iy_end = get_ind(y, yc, Nprocy)
    zi, iz_start, iz_end = get_ind(z, zc, Nprocz)

    x_start = ix_start[num_i[:]]
    y_start = iy_start[num_j[:]]
    z_start = iz_start[num_k[:]]
    x_end = ix_end[num_i[:]]
    y_end = iy_end[num_j[:]]
    z_end = iz_end[num_k[:]]

    # Loop over all processors partition
    for n in 1:Nproc
        # Extract coordinates for current processor

        part_x = ustrip.(Grid.x.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_y = ustrip.(Grid.y.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_z = ustrip.(Grid.z.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_phs = Phases[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]]
        part_T = Temp[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]]
        num_particles = size(part_x, 1) * size(part_x, 2) * size(part_x, 3)
        println("x_start[n]= $x_start[n]","x_end[n]= $x_end[n]","y_start[n]= $y_start[n]","y_end[n]= $y_end[n]","z_start[n]= $z_start[n]","z_end[n]= $z_end[n]")
        println("size(part_x) = $(size(part_x))","size(part_y) = $(size(part_y))","size(part_z) = $(size(part_z))","size(part_phs) = $(size(part_phs))","size(part_T) = $(size(part_T))")
 
        # Information vector per processor
        num_prop = 5       # number of properties we save [x/y/z/phase/T]
        lvec_info = num_particles

        lvec_prtcls = zeros(Float64, num_prop * num_particles)

        lvec_prtcls[1:num_prop:end] = part_x[:]
        lvec_prtcls[2:num_prop:end] = part_y[:]
        lvec_prtcls[3:num_prop:end] = part_z[:]
        lvec_prtcls[4:num_prop:end] = part_phs[:]
        lvec_prtcls[5:num_prop:end] = part_T[:]

        # Write output files
        if ~isdir(directory)
            mkdir(directory)
        end         # Create dir if not existent
        fname = @sprintf "%s/mdb.%1.8d.dat"  directory (n - 1)    # Name
        if verbose
            println("Writing LaMEM marker file -> $fname")                   # print info
        end
        lvec_output = [lvec_info; lvec_prtcls]           # one vec with info about length

        PetscBinaryWrite_Vec(fname, lvec_output)            # Write PETSc vector as binary file

    end
    return
end


# Internal routine to retrieve indices of local portion of the grid
function get_ind(x, xc, Nprocx)
    if Nprocx == 1
        xi = length(x)
        ix_start = [1]
        ix_end = [length(x)]
    else

        xi = zeros(Int64, Nprocx)
        for k in 1:Nprocx
            if k == 1
                xi[k] = length(x[(x .>= xc[k]) .& (x .<= xc[k + 1])])
            else
                xi[k] = length(x[(x .> xc[k]) .& (x .<= xc[k + 1])])
            end
        end
        ix_start = cumsum([0; xi[1:(end - 1)]]) .+ 1
        ix_end = cumsum(xi[1:end])
    end


    return xi, ix_start, ix_end
end

# Same as get_ind but without the need for the x vector
function get_ind2(dx,xc,Nprocx)

    if Nprocx == 1

        xi       = [xc[end] - xc[1]]/dx;
        ix_start = [1];
        ix_end   = [xc[end] - xc[1]]/dx;

    else
        xi = zeros(Int64, Nprocx)
        for k = 1:Nprocx
            xi[k] = round((xc[k+1] - xc[k])/dx);
        end

        ix_start = cumsum([0; xi[1:(end - 1)]]) .+ 1
        ix_end = cumsum(xi[1:end])

    end

    return xi,ix_start,ix_end

end

# Internal routine
function get_numscheme(Nprocx, Nprocy, Nprocz)
    n = zeros(Int64, Nprocx * Nprocy * Nprocz)
    nix = zeros(Int64, Nprocx * Nprocy * Nprocz)
    njy = zeros(Int64, Nprocx * Nprocy * Nprocz)
    nkz = zeros(Int64, Nprocx * Nprocy * Nprocz)

    num = 0
    for k in 1:Nprocz
        for j in 1:Nprocy
            for i in 1:Nprocx
                num = num + 1
                n[num] = num
                nix[num] = i
                njy[num] = j
                nkz[num] = k
            end
        end
    end

    return n, nix, njy, nkz
end


# Internal routine, to write a PETSc vector (as Float64)
"""
    PetscBinaryWrite_Vec(filename, A)

Writes a vector `A` to disk, such that it can be read with `PetscBinaryRead` (which assumes a Big Endian type)

"""
function PetscBinaryWrite_Vec(filename, A)

    # Note: use "hton" to transfer to Big Endian type, which is what PETScBinaryRead expects
    return open(filename, "w+") do f
        n = length(A)
        nummark = A[1]            # number of markers

        write(f, hton(Float64(1211214)))     # header (not actually used)
        write(f, hton(Float64(nummark)))     # info about # of markers written

        for i in 2:n
            write(f, hton(Float64(A[i])))   # Write data itself
        end

    end


end

"""
    nProcX,nProcY,nProcZ, xc,yc,zc, nNodeX,nNodeY,nNodeZ = get_processor_partitioning(filename; is64bit=false)

Reads a LaMEM processor partitioning file, used to create marker files, and returns the parallel layout.
By default this is done for a 32bit PETSc installation, which will fail if you actually use a 64bit version.

"""
function get_processor_partitioning(filename; is64bit = false)

    if is64bit
        typ = Int64
    else
        typ = Int32
    end
    io = open(filename, "r")


    nProcX = ntoh(read(io, typ))
    nProcY = ntoh(read(io, typ))
    nProcZ = ntoh(read(io, typ))

    nNodeX = ntoh(read(io, typ))
    nNodeY = ntoh(read(io, typ))
    nNodeZ = ntoh(read(io, typ))

    iX = [ntoh(read(io, typ)) for i in 1:(nProcX + 1)]
    iY = [ntoh(read(io, typ)) for i in 1:(nProcY + 1)]
    iZ = [ntoh(read(io, typ)) for i in 1:(nProcZ + 1)]

    CharLength = ntoh(read(io, Float64))
    xcoor = [ntoh(read(io, Float64)) for i in 1:nNodeX] .* CharLength
    ycoor = [ntoh(read(io, Float64)) for i in 1:nNodeY] .* CharLength
    zcoor = [ntoh(read(io, Float64)) for i in 1:nNodeZ] .* CharLength

    xc = xcoor[iX .+ 1]
    yc = ycoor[iY .+ 1]
    zc = zcoor[iZ .+ 1]

    close(io)

    return nProcX, nProcY, nProcZ,
        xc, yc, zc,
        nNodeX, nNodeY, nNodeZ

end


"""
    coord, Data_3D_Arrays, Name_Vec = read_data_VTR(fname)

Reads a VTR (structured grid) VTK file `fname` and extracts the coordinates, data arrays and names of the data.
In general, this only contains a piece of the data, and one should open a `*.pvtr` file to retrieve the full data
"""
function read_data_VTR(fname, FullSize)
    file = open(fname, "r")

    header = true
    num = 1
    CoordOffset = zeros(Int64, 3)
    Offset_Vec = []
    Name_Vec = [];     Type_Vec = []
    NumComp_Vec = [];     PieceExtent = []; WholeExtent = []
    while header == true

        line = readline(file)
        line_strip = lstrip(line)
        if startswith(line_strip, "<RectilinearGrid WholeExtent")
            id_start = findfirst("\"", line_strip)[1] + 1
            id_end = findlast("\"", line_strip)[1] - 1
            WholeExtent = parse.(Int64, split(line_strip[id_start:id_end]))
        end
        if startswith(line_strip, "<Piece Extent=")
            id_start = findfirst("\"", line_strip)[1] + 1
            id_end = findlast("\"", line_strip)[1] - 1
            PieceExtent = parse.(Int64, split(line_strip[id_start:id_end]))


        end
        if startswith(line_strip, "<Coordinates>")
            # Read info where the coordinates are stored
            Type, Name, NumberOfComponents, CoordOffset[1] = Parse_VTR_Line(readline(file)); num += 1
            Type, Name, NumberOfComponents, CoordOffset[2] = Parse_VTR_Line(readline(file)); num += 1
            Type, Name, NumberOfComponents, CoordOffset[3] = Parse_VTR_Line(readline(file)); num += 1
        end

        if startswith(line_strip, "<PointData>")
            line_strip = lstrip(readline(file))
            while ~startswith(line_strip, "</PointData>")
                Type, Name, NumberOfComponents, Offset = Parse_VTR_Line(line_strip);  num += 1

                Offset_Vec = [Offset_Vec;  Offset]
                Name_Vec = [Name_Vec;    Name]
                Type_Vec = [Type_Vec;   Type]
                NumComp_Vec = [NumComp_Vec; NumberOfComponents]
                line_strip = lstrip(readline(file))
            end
        end

        if startswith(line_strip, "<CellData>")
            line_strip = lstrip(readline(file))
            while ~startswith(line_strip, "</CellData>")
                Type, Name, NumberOfComponents, Offset = Parse_VTR_Line(line_strip);  num += 1

                Offset_Vec = [Offset_Vec;  Offset]
                Name_Vec = [Name_Vec;    Name]
                Type_Vec = [Type_Vec;   Type]
                NumComp_Vec = [NumComp_Vec; NumberOfComponents]
                line_strip = lstrip(readline(file))
                # if we have cell Data, for some reason we need to increment this by one.
                PieceExtent[1:2:end] .+= 1
            end

        end

        if startswith(line_strip, "<AppendedData ")
            header = false
        end

        num += 1
    end

    # Skip to beginning of raw data (linebreak)
    skip(file, 5)
    start_bin = position(file)      # start of binary data

    # Determine the end of the raw data
    seekend(file)
    skip(file, -29)

    end_bin = position(file)

    # Start with reading the coordinate arrays:
    coord_x = ReadBinaryData(file, start_bin, CoordOffset[1], (PieceExtent[2] - PieceExtent[1] + 1) * sizeof(Float32))
    coord_y = ReadBinaryData(file, start_bin, CoordOffset[2], (PieceExtent[4] - PieceExtent[3] + 1) * sizeof(Float32))
    coord_z = ReadBinaryData(file, start_bin, CoordOffset[3], (PieceExtent[6] - PieceExtent[5] + 1) * sizeof(Float32))


    # Read data arrays:
    Data_3D_Arrays = []
    ix = PieceExtent[1]:PieceExtent[2]
    iy = PieceExtent[3]:PieceExtent[4]
    iz = PieceExtent[5]:PieceExtent[6]
    numPoints = length(ix) * length(iy) * length(iz)

    coord_x_full = zeros(Float64, FullSize[1])
    coord_y_full = zeros(Float64, FullSize[2])
    coord_z_full = zeros(Float64, FullSize[3])

    coord_x_full[ix] = coord_x[1:length(ix)]
    coord_y_full[iy] = coord_y[1:length(iy)]
    coord_z_full[iz] = coord_z[1:length(iz)]

    for i in 1:(length(Name_Vec) - 1)

        data3D = ReadBinaryData(file, start_bin, Offset_Vec[i], numPoints * NumComp_Vec[i] * sizeof(Float32))
        data3D = getArray(data3D, PieceExtent, NumComp_Vec[i])

        data3D_full = zeros(Float64, NumComp_Vec[i], FullSize[1], FullSize[2], FullSize[3])  # Generate full data

        # ugly hack to make it work with parallel files
        ix_left = ix; ix_right = 1:length(ix_left)
        iy_left = iy; iy_right = 1:length(iy_left)
        iz_left = iz; iz_right = 1:length(iz_left)
        if ix_left[1] > 1
            ix_left = ix_left[2:end]; ix_right = ix_right[2:end]
        end
        if iy_left[1] > 1
            iy_left = iy_left[2:end]; iy_right = iy_right[2:end]
        end
        if iz_left[1] > 1
            iz_left = iz_left[2:end]; iz_right = iz_right[2:end]
        end

        data3D_full[1:NumComp_Vec[i], ix_left, iy_left, iz_left] = data3D[1:NumComp_Vec[i], ix_right, iy_right, iz_right]
        #data3D_full[1:NumComp_Vec[i], ix, iy, iz]        = data3D;

        Data_3D_Arrays = [Data_3D_Arrays; data3D_full]
    end
    i = length(Name_Vec)

    if Type_Vec[i] == "UInt8"
        data3D = ReadBinaryData(file, start_bin, Offset_Vec[i], numPoints * NumComp_Vec[i] * sizeof(UInt8), DataType = UInt8)
    else
        data3D = ReadBinaryData(file, start_bin, Offset_Vec[i], numPoints * NumComp_Vec[i] * sizeof(Float32))
    end

    data3D = getArray(data3D, PieceExtent, NumComp_Vec[i])
    data3D_full = zeros(Float64, NumComp_Vec[i], FullSize[1], FullSize[2], FullSize[3])  # Generate full d
    data3D_full[1:NumComp_Vec[i], ix, iy, iz] = data3D[1:NumComp_Vec[i], 1:length(ix), 1:length(iy), 1:length(iz)]

    Data_3D_Arrays = [Data_3D_Arrays; data3D_full]

    return coord_x_full, coord_y_full, coord_z_full, Data_3D_Arrays, Name_Vec, NumComp_Vec, ix, iy, iz
end

# Parses a line of a *.vtr file & retrieve Type/Name/NumberOfComponents/Offset
function Parse_VTR_Line(line)
    line_strip = lstrip(line)

    # Retrieve Type
    if findfirst("type", line_strip) != nothing
        id_start = findfirst("type", line_strip)[1] + 6

        line_strip = line_strip[id_start:end]
        id_end = findfirst("\"", line_strip)[1] - 1
        Type = line_strip[1:id_end]
        line_strip = line_strip[id_end:end]
    else
        Type = nothing
    end

    # Retrieve Name
    if findfirst("Name", line_strip) != nothing
        id_start = findfirst("Name", line_strip)[1] + 6
        line_strip = line_strip[id_start:end]
        id_end = findfirst("\"", line_strip)[1] - 1
        Name = line_strip[1:id_end]
        line_strip = line_strip[id_end:end]
    else
        Name = nothing
    end

    # Retrieve number of components
    if findfirst("NumberOfComponents", line_strip) != nothing
        id_start = findfirst("NumberOfComponents", line_strip)[1] + 20
        line_strip = line_strip[id_start:end]
        id_end = findfirst("\"", line_strip)[1] - 1
        NumberOfComponents = parse(Int64, line_strip[1:id_end])
        line_strip = line_strip[id_end:end]
    else
        NumberOfComponents = nothing
    end

    # Offset
    if findfirst("offset", line_strip) != nothing
        id_start = findfirst("offset", line_strip)[1] + 8
        line_strip = line_strip[id_start:end]
        id_end = findfirst("\"", line_strip)[1] - 1
        Offset = parse(Int64, line_strip[1:id_end])
    else
        Offset = nothing
    end
    return Type, Name, NumberOfComponents, Offset
end


function getArray(data, PieceExtent, NumComp)
    data = reshape(data, (NumComp, PieceExtent[2] - PieceExtent[1] + 1, PieceExtent[4] - PieceExtent[3] + 1, PieceExtent[6] - PieceExtent[5] + 1))
    return data
end

function ReadBinaryData(file::IOStream, start_bin::Int64, Offset::Int64, BytesToRead; DataType = Float32)

    seekstart(file)                                 # go to start
    skip(file, start_bin + Offset)                    # move to beginning of raw binary data
    buffer = read(file, BytesToRead)          # Read necesaary bytes
    data = reinterpret(DataType, buffer)    # Transfer to buffer

    data = Float64.(data[1:end])         # Transfer to Float64
    return data
end

"""
    Data::ParaviewData = read_data_PVTR(fname, dir)

Reads a parallel, rectilinear, `*.vts` file with the name `fname` and located in `dir` and create a 3D `Data` struct from it.

# Example
```julia-repl
julia> Data = read_data_PVTR("Haaksbergen.pvtr", "./Timestep_00000005_3.35780500e-01/")
ParaviewData
  size  : (33, 33, 33)
  x     ϵ [ -3.0 : 3.0]
  y     ϵ [ -2.0 : 2.0]
  z     ϵ [ -2.0 : 0.0]
  fields: (:phase, :density, :visc_total, :visc_creep, :velocity, :pressure, :temperature, :dev_stress, :strain_rate, :j2_dev_stress, :j2_strain_rate, :plast_strain, :plast_dissip, :tot_displ, :yield, :moment_res, :cont_res)
```
"""
function read_data_PVTR(fname, dir)
    file = open(joinpath(dir, fname), "r")

    header = true
    num = 1
    FullSize = (1, 1, 1)
    num_data_sets = 1
    Data_3D = []; coord_x = []; coord_y = []; coord_z = []; NumComp = []; Names = []
    while header == true

        line = readline(file)
        line_strip = lstrip(line)
        if startswith(line_strip, "<PRectilinearGrid")
            id_start = findfirst("WholeExtent=", line_strip)[1] + 13
            line_strip = line_strip[id_start:end]
            id_end = findfirst("\"", line_strip)[1] - 1
            line_piece = line_strip[1:id_end]

            WholeExtent = parse.(Int64, split(line_piece))
            FullSize = (WholeExtent[2], WholeExtent[4], WholeExtent[6])
        end


        if startswith(line_strip, "<Piece")
            id_start = findfirst("Source=", line_strip)[1] + 8
            line_strip = line_strip[id_start:end]
            id_end = findfirst("\"", line_strip)[1] - 1
            fname_piece = line_strip[1:id_end]

            if num_data_sets == 1
                coord_x, coord_y, coord_z, Data_3D, Names, NumComp, ix, iy, iz = read_data_VTR(joinpath(dir, fname_piece), FullSize)
            else
                coord_x1, coord_y1, coord_z1, Data_3D1, Names, NumComp, ix, iy, iz = read_data_VTR(joinpath(dir, fname_piece), FullSize)
                coord_x[ix] = coord_x1[ix]
                coord_y[iy] = coord_y1[iy]
                coord_z[iz] = coord_z1[iz]

                Data_3D = Data_3D + Data_3D1

            end
            num_data_sets += 1
        end


        if startswith(line_strip, "</PRectilinearGrid")
            header = false
        end
    end

    # Create a named-Tuple out of the fields
    NamesSymbol = []
    for i in 1:length(Names)
        id = findfirst(" ", Names[i])
        if id == nothing
            Names_Strip = Names[i]
        else
            Names_Strip = Names[i][1:(findfirst(" ", Names[i])[1] - 1)]
        end
        NamesSymbol = [NamesSymbol; Names_Strip]
    end

    #  NamesSymbol =   [Names[i][1:findfirst(" ", Names[i])[1]-1] for i=1:length(Names)]
    Names1 = Symbol.(NamesSymbol)

    Data_Array = []
    num = 1
    for i in 1:length(NumComp)
        data = Data_3D[num:(num + NumComp[i] - 1), :, :, :]
        data_arrays = [data[i, :, :, :] for i in 1:size(data, 1)]
        data_tuple = tuple(data_arrays...)

        if size(data, 1) > 1
            Data_NamedTuple = NamedTuple{(Names1[i],)}((data_tuple,))
        else
            Data_NamedTuple = NamedTuple{(Names1[i],)}((data_tuple[1],))
        end
        Data_Array = [Data_Array; Data_NamedTuple]

        num = num + NumComp[i]
    end

    # Merge vector with tuples into a NamedTuple
    fields = Data_Array[1]
    for i in 2:length(Data_Array)
        fields = merge(fields, Data_Array[i])
    end

    # Create a ParaviewData struct from it.
    X, Y, Z = xyz_grid(coord_x, coord_y, coord_z)
    DataC = ParaviewData(X, Y, Z, fields)

    return DataC
end

"""
    save_LaMEM_topography(Topo::CartData, filename::String)

This writes a topography file `Topo` for use in LaMEM, which should have size `(nx,ny,1)` and contain the field `:Topography`
"""
function save_LaMEM_topography(Topo::CartData, filename::String)

    if (size(Topo.z.val, 3) != 1)
        error("Not a valid `CartData' Topography file (size in 3rd dimension should be 1)")
    end
    if !haskey(Topo.fields, :Topography)
        error("The topography `CartData` structure requires a field :Topography")
    end

    # get grid properties
    nx = Float64(size(Topo.fields.Topography, 1))
    ny = Float64(size(Topo.fields.Topography, 2))
    x0 = ustrip(Topo.x.val[1, 1, 1])
    y0 = ustrip(Topo.y.val[1, 1, 1])

    # LaMEM wants a uniform grid, so interpolate if necessary
    if length(unique(trunc.(diff(Topo.x.val[:, 1, 1]), digits = 8))) > 1 || length(unique(trunc.(diff(Topo.y.val[1, :, 1]), digits = 8))) > 1
        x1 = ustrip(Topo.x.val[end, 1, 1])
        y1 = ustrip(Topo.y.val[1, end, 1])
        dx = (x1 - x0) / (nx - 1)
        dy = (y1 - y0) / (ny - 1)

        itp = LinearInterpolation((Topo.x.val[:, 1, 1], Topo.y.val[1, :, 1]), ustrip.(Topo.fields.Topography[:, :, 1]))
        Topo_itp = [itp(x, y) for x in x0:dx:x1, y in y0:dy:y1]

        # Code the topograhic data into a vector
        Topo_vec = [ nx;ny;x0;y0;dx;dy; Topo_itp[:]]
    else
        dx = ustrip(Topo.x.val[2, 2, 1]) - x0
        dy = ustrip(Topo.y.val[2, 2, 1]) - y0
        # Code the topograhic data into a vector
        Topo_vec = [ nx;ny;x0;y0;dx;dy; ustrip.(Topo.fields.Topography[:])]
    end

    # Write as PetscBinary file
    PetscBinaryWrite_Vec(filename, Topo_vec)

    println("Written LaMEM topography file: $(filename)")

    return nothing
end

"""
    create_partitioning_file(LaMEM_input::String, NumProc::Int64; LaMEM_dir::String=pwd(), LaMEM_options::String="", MPI_dir="", verbose=true)

This executes LaMEM for the input file `LaMEM_input` & creates a parallel partitioning file for `NumProc` processors.
The directory where the LaMEM binary is can be specified; if not it is assumed to be in the current directory.
Likewise for the `mpiexec` directory (if not specified it is assumed to be available on the command line).

"""
function create_partitioning_file(LaMEM_input::String, NumProc::Int64; LaMEM_dir::String = pwd(), LaMEM_options = "", MPI_dir = "", verbose = true)

    # Create string to execute LaMEM
    mpi_str = MPI_dir * "mpiexec -n $(NumProc) "
    LaMEM_str = LaMEM_dir * "/" * "LaMEM -ParamFile " * LaMEM_input * " -mode save_grid "
    str = mpi_str * LaMEM_str

    if verbose == true
        println("Executing command: $str")
    end
    # Run
    exit = run(`sh -c $str`, wait = false)

    # Retrieve newest file
    if success(exit)
        files = readdir(glob"ProcessorPartitioning_*.bin")
        time_modified = zeros(length(files))
        for (i, file) in enumerate(files)
            time_modified[i] = stat(file).mtime
        end
        id = findall(time_modified .== maximum(time_modified))   # last modified
        PartFile = files[id]
        if verbose == true
            println("Successfully generated PartitioningFile: $(PartFile[1])")
        end
    else
        error("Something went wrong with executing command ")
    end

    return PartFile[1]

end

"""
    X,Y,Z = coordinate_grids(Data::LaMEM_grid; cell=false)

Returns 3D coordinate arrays
"""
function coordinate_grids(Data::LaMEM_grid; cell = false)
    X, Y, Z = Data.X, Data.Y, Data.Z
    if cell
        X, Y, Z = average_q1(X), average_q1(Y), average_q1(Z)
    end

    return X, Y, Z
end

function Base.show(io::IO, d::LaMEM_partitioning_info)

    println(io, "LaMEM Partitioning info: ")
    println(io, "  nProcX : $(d.nProcX)")
    println(io, "  nProcY : $(d.nProcY)")
    println(io, "  nProcZ : $(d.nProcZ)")
    println(io, "  nNodeX : $(d.nNodeX)")
    println(io, "  nNodeY : $(d.nNodeY)")
    println(io, "  nNodeZ : $(d.nNodeZ)")
    println(io, "  xc     : $(d.xc)")
    println(io, "  yc     : $(d.yc)")

    return println(io, "  zc     : $(d.zc)")

end

function check_markers_directory(directory)

    if !isdir(directory)
        mkdir(directory)
    end

end

function get_LaMEM_grid_info(file; args::Union{String, Nothing} = nothing)

    # read information from file
    nmark_x = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nmark_x", Int64, args = args)
    nmark_y = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nmark_y", Int64, args = args)
    nmark_z = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nmark_z", Int64, args = args)

    nel_x = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nel_x", Int64, args = args)
    nel_y = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nel_y", Int64, args = args)
    nel_z = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "nel_z", Int64, args = args)

    parsed_x = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "coord_x", Float64, args = args)
    parsed_y = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "coord_y", Float64, args = args)
    parsed_z = GeophysicalModelGenerator.ParseValue_LaMEM_InputFile(file, "coord_z", Float64, args = args)

    # compute information from file
    W = parsed_x[end] - parsed_x[1]
    L = parsed_y[end] - parsed_y[1]
    H = parsed_z[end] - parsed_z[1]

    nel_x_tot = sum(nel_x)
    nel_y_tot = sum(nel_y)
    nel_z_tot = sum(nel_z)

    nump_x = nel_x_tot * nmark_x
    nump_y = nel_y_tot * nmark_y
    nump_z = nel_z_tot * nmark_z

    # finish Grid
    Grid = LaMEM_grid(
        nmark_x, nmark_y, nmark_z,
        nump_x, nump_y, nump_z,
        nel_x,   nel_y,   nel_z,
        W, L, H,
        parsed_x, parsed_y, parsed_z,
        [],[],[],
        [],[],[],
        [],[],[],
        [],[],[]
    )

    return Grid

end


"""
    p_dist = get_particles_distribution(Grid, P)
 Get the distribution of particles in the grid
    Grid: LaMEM_grid
    P:    LaMEM_partitioning_info
 Returns a LaMEM_partitioning_info object with the distribution of particles in the grid
"""
function get_particles_distribution(Grid,P)

    # get number of processors and processor coordnate bounds
    nProcX = P.nProcX;
    nProcY = P.nProcY;
    nProcZ = P.nProcZ;
    xc     = P.xc;
    yc     = P.yc;
    zc     = P.zc;

    (num, num_i, num_j, num_k) = get_numscheme(nProcX, nProcY, nProcZ);

    dx     = Grid.W/Grid.nump_x;
    dy     = Grid.L/Grid.nump_y;
    dz     = Grid.H/Grid.nump_z;

    # % Get particles of respective procs
    # % xi - amount of particles in x direction in each core
    # % ix_start - indexes where they start for each core
    (xi,ix_start,ix_end) = get_ind2(dx,xc,nProcX);
    (yi,iy_start,iy_end) = get_ind2(dy,yc,nProcY);
    (zi,iz_start,iz_end) = get_ind2(dz,zc,nProcZ);

    x_start = ix_start[num_i[:]]
    y_start = iy_start[num_j[:]]
    z_start = iz_start[num_k[:]]
    x_end = ix_end[num_i[:]]
    y_end = iy_end[num_j[:]]
    z_end = iz_end[num_k[:]]

    p_dist = particles_distribution(x_start,x_end,y_start,y_end,z_start,z_end);

    return p_dist

end

"""
    Grid = get_proc_grid(Grid_info, p_dist, proc_bounds, proc_num, RandomNoise)
 Get the local grid for the current processor
    Grid_info:   LaMEM_grid
    p_dist:      LaMEM_partitioning_info
    proc_bounds: bounds of the current processor
    proc_num:    processor number
    RandomNoise: add random noise to the grid (0/1)
 Returns a LaMEM_grid object with the local grid for the current processor
"""
function get_proc_grid(Grid_info,p_dist,proc_bounds,proc_num,RandomNoise)

    x_proc_bound  = proc_bounds[1];
    y_proc_bound  = proc_bounds[2];
    z_proc_bound  = proc_bounds[3];

    loc_nump_x = p_dist.x_end[proc_num] - p_dist.x_start[proc_num] + 1
    loc_nump_y = p_dist.y_end[proc_num] - p_dist.y_start[proc_num] + 1
    loc_nump_z = p_dist.z_end[proc_num] - p_dist.z_start[proc_num] + 1
    
    loc_nel_x = loc_nump_x/Grid_info.nmark_x
    loc_nel_y = loc_nump_y/Grid_info.nmark_y
    loc_nel_z = loc_nump_z/Grid_info.nmark_z

    x  = range(x_proc_bound[1], x_proc_bound[2], length=loc_nump_x)
    y  = range(y_proc_bound[1], y_proc_bound[2], length=loc_nump_y)
    z  = range(z_proc_bound[1], z_proc_bound[2], length=loc_nump_z)

    # marker grid
    X, Y, Z = GeophysicalModelGenerator.xyz_grid(x, y, z)

    W = x_proc_bound[2] - x_proc_bound[1]
    L = y_proc_bound[2] - y_proc_bound[1]
    H = z_proc_bound[2] - z_proc_bound[1]

    if RandomNoise == 1
        dx = x[2]   - x[1]
        dy = y[2]   - y[1]
        dz = z[2]   - z[1]
        dXNoise = zeros(size(X)) + dx;
        dYNoise = zeros(size(Y)) + dy;
        dZNoise = zeros(size(Z)) + dz;
    
        dXNoise = dXNoise.*(rand(size(dXNoise))-0.5);
        dYNoise = dYNoise.*(rand(size(dYNoise))-0.5);
        dZNoise = dZNoise.*(rand(size(dZNoise))-0.5);
    
        Xpart   = X + dXNoise;
        Ypart   = Y + dYNoise;
        Zpart   = Z + dZNoise;
    
        X       = Xpart;
        Y       = Ypart;
        Z       = Zpart;
        x       = X(1,:,1);
        y       = Y(:,1,1);
        z       = Z(1,1,:);
    
    end

    Grid = LaMEM_grid(
        Grid_info.nmark_x, Grid_info.nmark_y, Grid_info.nmark_z,
        loc_nump_x, loc_nump_y, loc_nump_z,
        loc_nel_x, loc_nel_y, loc_nel_z,
        W, L, H,
        x, y, z,
        x, y, z,
        X, Y, Z,
        [], [], [],
        [], [], []
    )
    return Grid

end

"""
proc_bounds = get_proc_bound(Grid, p_dist, proc_num)
 Get the bounds of the current processor in x, y, z direction
    Grid:       LaMEM_grid
    p_dist:     LaMEM_partitioning_info
    proc_num:   processor number
 Returns a 3 element vector with maximum and minimum values[[x_min, x_max],[y_min, y_max],[z_min, z_max]] for current processor proc_num

# Example for a model with 8 MPI ranks for current processor number 2
And gives us coordinates for the current processor
```julia
Grid_example = LaMEM_grid(
        3, 3, 3,
        12, 12, 12,
        4,   4,   4,
        10, 5, 5,
        [10.0, 18.0], [20,28], [30,38],
        [],[],[],
        [],[],[],
        [],[],[],
        [],[],[]
    )
P_example      = setup_model_domain(Grid_example.coord_x, Grid_example.coord_y, Grid_example.coord_z, Grid_example.nel_x, Grid_example.nel_x, Grid_example.nel_x, 8)
p_dist_example = get_particles_distribution(Grid_example,P_example)
proc_bounds    = get_proc_bound(Grid_example,p_dist_example,2)
```
"""
function get_proc_bound(Grid,p_dist,proc_num)

    dx           = Grid.W/Grid.nump_x;
    dy           = Grid.L/Grid.nump_y;
    dz           = Grid.H/Grid.nump_z;

    parsed_x     = Grid.coord_x
    parsed_y     = Grid.coord_y
    parsed_z     = Grid.coord_z

    model_x      = [ parsed_x[1] + dx/2, parsed_x[end] - dx/2 ]
    model_y      = [ parsed_y[1] + dy/2, parsed_y[end] - dy/2 ]
    model_z      = [ parsed_z[1] + dz/2, parsed_z[end] - dz/2 ]

    x_left       = model_x[1];
    y_front      = model_y[1];
    z_bot        = model_z[1];

    x_start      = p_dist.x_start;
    x_end        = p_dist.x_end;
    y_start      = p_dist.y_start;
    y_end        = p_dist.y_end;
    z_start      = p_dist.z_start;
    z_end        = p_dist.z_end;

    x_proc_bound = [ x_left  + dx*( x_start[proc_num] - 1 ), x_left  + dx*( x_end[proc_num] - 1 ) ];
    y_proc_bound = [ y_front + dy*( y_start[proc_num] - 1 ), y_front + dy*( y_end[proc_num] - 1 ) ];
    z_proc_bound = [ z_bot   + dz*( z_start[proc_num] - 1 ), z_bot   + dz*( z_end[proc_num] - 1 ) ];

    return [ x_proc_bound, y_proc_bound, z_proc_bound ]

end

"""
    crop_bounds(uncropped_bounds, proc_bounds, x, y, z)
    Crop boundaries from the whole model to only the extent of the current processor
    uncropped_bounds: 3 element vector with maximum and minimum values[[x_min, x_max],[y_min, y_max],[z_min, z_max]] for the whole model
    proc_bounds:      3 element vector with maximum and minimum values[[x_min, x_max],[y_min, y_max],[z_min, z_max]] for the current processor
    x,y,z:           3 element vector with maximum and minimum values[[x_min, x_max],[y_min, y_max],[z_min, z_max]] for the current processor
    Returns a 3 element vector with maximum and minimum values[[x_min, x_max],[y_min, y_max],[z_min, z_max]] for the current processor
"""
function crop_bounds(uncropped_bounds, proc_bounds, x, y, z)

    # Crop boundaries from the whole model to only the extent of the current processor
    vecs       = [x, y, z]  
    new_bounds = [zeros(size(vecs[i])) for i in eachindex(vecs)]
    for i in eachindex(vecs)
        vec = vecs[i]
        test_bound = uncropped_bounds[i]
        mask_bound = proc_bounds[i]

        new_bound = Float64[]

        if test_bound[1] < test_bound[2]
            if test_bound[1] <= mask_bound[1]
                if test_bound[2] >= mask_bound[2]
                    new_bound = [mask_bound[1], mask_bound[2]]
                elseif test_bound[2] >= mask_bound[1]
                    new_bound = [mask_bound[1], test_bound[2]]
                end
            end

            if test_bound[1] >= mask_bound[1]
                if test_bound[2] <= mask_bound[2]
                    new_bound = [closest_val(test_bound[1], vec), closest_val(test_bound[2], vec)]
                elseif test_bound[1] <= mask_bound[2] && test_bound[2] >= mask_bound[2]
                    new_bound = [test_bound[1], mask_bound[2]]
                end
            end
        else
            error("Wrong coordinates assignment")
        end

        if isempty(new_bound)
            return []
        end
        new_bounds[i] = new_bound
    end

    return new_bounds

end

function closest_val(val, vec)
    return  vec[argmin(abs.(vec .- val))]
end

"""
    decompose_mpi_ranks(total_ranks::Int, nx::Int, ny::Int, nz::Int) -> Tuple{Int,Int,Int}

Decompose total number of MPI ranks into a 3D processor grid (px, py, pz),
optimizing for cell aspect ratio closest to 1.0.
"""
function decompose_mpi_ranks(total_ranks::Int, nx::Int, ny::Int, nz::Int)
    # Get all factors of total_ranks
    factors = get_factors(total_ranks)
    
    # Initialize best configuration
    best_px = 1
    best_py = 1
    best_pz = 1
    best_metric = Inf
    
    # Try all possible combinations of factors
    for px in factors
        remaining = total_ranks ÷ px
        rem_factors = get_factors(remaining)
        
        for py in rem_factors
            pz = remaining ÷ py
            
            # Skip invalid combinations
            if px * py * pz != total_ranks
                continue
            end
            
            # Calculate local grid sizes
            local_nx = nx / px
            local_ny = ny / py
            local_nz = nz / pz
            
            # Calculate aspect ratios (always ≥ 1.0)
            ar_xy = max(local_nx/local_ny, local_ny/local_nx)
            ar_xz = max(local_nx/local_nz, local_nz/local_nx)
            ar_yz = max(local_ny/local_nz, local_nz/local_ny)
            
            # Metric: average deviation from aspect ratio of 1.0
            metric = (ar_xy + ar_xz + ar_yz) / 3.0
            
            # Update best configuration if this one is better
            # If metrics are equal, prefer larger px
            if metric < best_metric || 
               (isapprox(metric, best_metric, rtol=1e-10) && px > best_px)
                best_metric = metric
                best_px = px
                best_py = py
                best_pz = pz
            end
        end
    end
    
    # Calculate and print aspect ratios for chosen decomposition
    local_nx = nx / best_px
    local_ny = ny / best_py
    local_nz = nz / best_pz

    println("Maximum aspect ratio: $(max(local_nx/local_ny, local_ny/local_nx,
                                    local_nx/local_nz, local_nz/local_nx,
                                    local_ny/local_nz, local_nz/local_ny))")
    
    return (best_px, best_py, best_pz)
end

"""
Get all factors of a number n
"""
function get_factors(n::Int)
    factors = Int[]
    for i in 1:isqrt(n)
        if n % i == 0
            push!(factors, i)
            if i != n÷i
                push!(factors, n÷i)
            end
        end
    end
    sort!(factors)
    return factors
end

"""
    setup_model_domain(coord_x::Vector{Float64}, 
                      coord_y::Vector{Float64}, 
                      coord_z::Vector{Float64},
                      nx::Int, ny::Int, nz::Int, 
                      n_ranks::Int) -> ModelDomain

Setup model domain decomposition using domain boundaries and resolution.

Parameters:
- coord_x, coord_y, coord_z: 2-element vectors specifying [min, max] for each direction
- nx, ny, nz: Number of cells in each direction
- n_ranks: Total number of MPI ranks

Returns:
- ModelDomain struct containing all domain decomposition information
"""
function setup_model_domain(coord_x::AbstractVector{<:Real}, 
                            coord_y::AbstractVector{<:Real}, 
                            coord_z::AbstractVector{<:Real},
                            nx::Int, ny::Int, nz::Int,
                            n_ranks::Int)
    
    # Verify input vectors have correct size
    if any(length.([coord_x, coord_y, coord_z]) .!= 2)
        error("coord_x, coord_y, and coord_z must be 2-element vectors [min, max]")
    end
    
    # Generate full coordinate vectors
    nnodx = nx + 1
    nnody = ny + 1
    nnodz = nz + 1
    
    xcoor = collect(range(coord_x[1], coord_x[2], length=nnodx))
    ycoor = collect(range(coord_y[1], coord_y[2], length=nnody))
    zcoor = collect(range(coord_z[1], coord_z[2], length=nnodz))
    
    # Decompose MPI ranks into 3D processor grid
    Nprocx, Nprocy, Nprocz = decompose_mpi_ranks(n_ranks, nx, ny, nz)
    
    # Calculate subdomain divisions
    function calculate_domain_divisions(N::Int, nproc::Int)
        base_size = div(N, nproc)
        remainder = N % nproc
        
        indices = zeros(Int, nproc + 1)
        indices[1] = 1
        
        for i in 1:nproc
            local_size = base_size + (i <= remainder ? 1 : 0)
            indices[i + 1] = indices[i] + local_size
        end
        
        return indices
    end
    
    # Calculate divisions for each direction
    ix = calculate_domain_divisions(nx, Nprocx)
    iy = calculate_domain_divisions(ny, Nprocy)
    iz = calculate_domain_divisions(nz, Nprocz)
    
    P = LaMEM_partitioning_info(
        Nprocx, Nprocy, Nprocz,
        nnodx, nnody, nnodz, 
        xcoor[ix], ycoor[iy],zcoor[iz]
        )

        xcoor=[]
        ycoor=[]
        zcoor=[]

    return P

end
