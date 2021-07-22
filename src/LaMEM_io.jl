using Base: Int64, Float64, NamedTuple
using Printf

# LaMEM I/O
# 
# These are routines that help to create a LaMEM marker files from a CartData structure, which can be used to perform geodynamic simulations
# We also include routines with which we can read LaMEM *.pvtr files into julia 

export LaMEM_grid, ReadLaMEM_InputFile
export Save_LaMEMMarkersParallel, GetProcessorPartitioning, ReadData_VTR, ReadData_PVTR

"""
Structure that holds information about the LaMEM grid (usually read from an input file).
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

# Example
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

Parses a LaMEM input file and stores grid information in the `Grid` structure.

# Example
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

    X,Y,Z   =   XYZGrid(x,y,z); # create 3D grid using regular spacing
    
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
    Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers", verbose=true)

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
function Save_LaMEMMarkersParallel(Grid::CartData; PartitioningFile=empty, directory="./markers", verbose=true)

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

        # uggly hack to make it work with parallel files    
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
    Data::CartData = ReadData_PVTR(fname, dir)

Reads a parallel, rectilinear, `*.vts` file with the name `fname` and located in `dir` and create a 3D `Data` struct from it.

# Example
```julia
julia> Data = ReadData_PVTR("Haaksbergen.pvtr", "./Timestep_00000005_3.35780500e-01/")
CartData 
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

    # Create a CartData struct from it.
    X,Y,Z       =   XYZGrid(coord_x, coord_y, coord_z)
    DataC       =   CartData(X,Y,Z, fields);

    return DataC
end
