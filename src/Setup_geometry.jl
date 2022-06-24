using Base: Int64, Float64, NamedTuple
using Printf
using Parameters        # helps setting default parameters in structures
using SpecialFunctions: erfc  

# Setup_geometry
# 
# These are routines that help to create input geometries, such as slabs with a given angle
#

export  AddBox!, AddSphere!, AddEllipsoid!, AddCylinder!,
        ConstantTemp, LinearTemp, HalfspaceCoolingTemp, SpreadingRateTemp,
        ConstantPhase, LithosphericPhases, 
        Compute_ThermalStructure, Compute_Phase


"""
    AddBox!(Phase, Temp, Grid::AbstractGeneralGrid; xlim=Tuple{2}, [ylim=Tuple{2}], zlim=Tuple{2},
            Origin=nothing, StrikeAngle=0, DipAngle=0,
            phase = ConstantPhase(1),
            T=nothing )

Adds a box with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- Phase - Phase array (consistent with Grid)
- Temp  - Temperature array (consistent with Grid)
- Grid -  grid structure (usually obtained with ReadLaMEM_InputFile, but can also be other grid types)
- xlim -  left/right coordinates of box
- ylim -  front/back coordinates of box [optional; if not specified we use the whole box]
- zlim -  bottom/top coordinates of box
- Origin - the origin, used to rotate the box around. Default is the left-front-top corner
- StrikeAngle - strike angle of slab
- DipAngle - dip angle of slab
- phase - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()` 
- T - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()` 


Examples
========

Example 1) Box with constant phase and temperature & a dip angle of 10 degrees:
```julia 
julia> Grid = ReadLaMEM_InputFile("test_files/SaltModels.dat")
LaMEM Grid: 
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=ConstantTemp(1000))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> Write_Paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview 
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"   
```

Example 2) Box with halfspace cooling profile
```julia 
julia> Grid = ReadLaMEM_InputFile("test_files/SaltModels.dat")
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> AddBox!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=ConstantTemp(1000))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> Write_Paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview 
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"   
```
"""
function AddBox!(Phase, Temp, Grid::AbstractGeneralGrid;                 # required input
                xlim=Tuple{2}, ylim=nothing, zlim=Tuple{2},     # limits of the box
                Origin=nothing, StrikeAngle=0, DipAngle=0,      # origin & dip/strike
                phase = ConstantPhase(1),                       # Sets the phase number(s) in the box
                T=nothing )                                     # Sets the thermal structure (various fucntions are available)
    
    # Retrieve 3D data arrays for the grid
    X,Y,Z = coordinate_grids(Grid)

    # Limits of block                
    if ylim==nothing 
        ylim = (minimum(Y), maximum(Y)) 
    end
    
    if Origin==nothing 
        Origin = (xlim[1], ylim[1], zlim[2])  # upper-left corner
    end

    # Perform rotation of 3D coordinates:
    Xrot = X .- Origin[1];
    Yrot = Y .- Origin[2];
    Zrot = Z .- Origin[3];

    Rot3D!(Xrot,Yrot,Zrot, StrikeAngle, DipAngle)


    # Set phase number & thermal structure in the full domain
    ztop = zlim[2] - Origin[3]
    zbot = zlim[1] - Origin[3]
    ind = findall(  (Xrot .>= (xlim[1] - Origin[1])) .& (Xrot .<= (xlim[2] - Origin[1])) .&  
                    (Yrot .>= (ylim[1] - Origin[2])) .& (Yrot .<= (ylim[2] - Origin[2])) .&  
                    (Zrot .>= zbot) .& (Zrot .<= ztop)  )
    

    # Compute thermal structure accordingly. See routines below for different options
    if T != nothing
        Temp[ind] = Compute_ThermalStructure(Temp[ind], Xrot[ind], Yrot[ind], Zrot[ind], T)
    end

    # Set the phase. Different routines are available for that - see below.
    Phase[ind] = Compute_Phase(Phase[ind], Temp[ind], Xrot[ind], Yrot[ind], Zrot[ind], phase)
    
    return nothing
end

"""
    AddSphere!(Phase, Temp, Grid::AbstractGeneralGrid; cen=Tuple{3}, radius=Tuple{1},
            phase = ConstantPhase(1).
            T=nothing )

Adds a sphere with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- Phase - Phase array (consistent with Grid)
- Temp  - Temperature array (consistent with Grid)
- Grid - LaMEM grid structure (usually obtained with ReadLaMEM_InputFile)
- cen - center coordinates of sphere
- radius - radius of sphere
- phase - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()` 
- T - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()` 


Example
========

Sphere with constant phase and temperature:
```julia 
julia> Grid = ReadLaMEM_InputFile("test_files/SaltModels.dat")
LaMEM Grid: 
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> AddSphere!(Phases,Temp,Grid, cen=(0,0,-1), radius=0.5, phase=ConstantPhase(2), T=ConstantTemp(800))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> Write_Paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview 
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"   
```
"""
function AddSphere!(Phase, Temp, Grid::AbstractGeneralGrid;      # required input
    cen=Tuple{3}, radius=Tuple{1},                         # center and radius of the sphere
    phase = ConstantPhase(1),                           # Sets the phase number(s) in the sphere
    T=nothing )                                         # Sets the thermal structure (various fucntions are available)
          
    # Retrieve 3D data arrays for the grid
    X,Y,Z = coordinate_grids(Grid)

    # Set phase number & thermal structure in the full domain
    ind = findall(((X .- cen[1]).^2 + (Y .- cen[2]).^2 + (Z .- cen[3]).^2).^0.5 .< radius)

    # Compute thermal structure accordingly. See routines below for different options
    if T != nothing
        Temp[ind] = Compute_ThermalStructure(Temp[ind], X[ind], Y[ind], Z[ind], T)
    end

    # Set the phase. Different routines are available for that - see below.
    Phase[ind] = Compute_Phase(Phase[ind], Temp[ind], X[ind], Y[ind], Z[ind], phase)

    return nothing
end

"""
    AddEllipsoid!(Phase, Temp, Grid::AbstractGeneralGrid; cen=Tuple{3}, axes=Tuple{3},
            Origin=nothing, StrikeAngle=0, DipAngle=0,
            phase = ConstantPhase(1).
            T=nothing )

Adds an Ellipsoid with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- Phase - Phase array (consistent with Grid)
- Temp  - Temperature array (consistent with Grid)
- Grid - LaMEM grid structure (usually obtained with ReadLaMEM_InputFile)
- cen - center coordinates of sphere
- axes - semi-axes of ellipsoid in X,Y,Z
- Origin - the origin, used to rotate the box around. Default is the left-front-top corner
- StrikeAngle - strike angle of slab
- DipAngle - dip angle of slab
- phase - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()` 
- T - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()` 


Example
========

Ellipsoid with constant phase and temperature, rotated 90 degrees and tilted by 45 degrees:
```julia 
julia> Grid = ReadLaMEM_InputFile("test_files/SaltModels.dat")
LaMEM Grid: 
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> AddEllipsoid!(Phases,Temp,Grid, cen=(-1,-1,-1), axes=(0.2,0.1,0.5), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(3), T=ConstantTemp(600))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> Write_Paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview 
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"   
```
"""
function AddEllipsoid!(Phase, Temp, Grid::AbstractGeneralGrid;      # required input
    cen=Tuple{3}, axes=Tuple{3},                           # center and semi-axes of the ellpsoid
    Origin=nothing, StrikeAngle=0, DipAngle=0,             # origin & dip/strike
    phase = ConstantPhase(1),                              # Sets the phase number(s) in the box
    T=nothing )                                            # Sets the thermal structure (various fucntions are available)

    if Origin==nothing 
        Origin = cen  # center
    end

    # Retrieve 3D data arrays for the grid
    X,Y,Z = coordinate_grids(Grid)

    # Perform rotation of 3D coordinates:
    Xrot = X .- Origin[1];
    Yrot = Y .- Origin[2];
    Zrot = Z .- Origin[3];

    Rot3D!(Xrot,Yrot,Zrot, StrikeAngle, DipAngle)

    # Set phase number & thermal structure in the full domain
    x2     = axes[1]^2
    y2     = axes[2]^2
    z2     = axes[3]^2
    cenRot = cen .- Origin
    ind = findall((((Xrot .- cenRot[1]).^2)./x2 + ((Yrot .- cenRot[2]).^2)./y2 +
                   ((Zrot .- cenRot[3]).^2)./z2) .^0.5 .<= 1)

    # Compute thermal structure accordingly. See routines below for different options
    if T != nothing
        Temp[ind] = Compute_ThermalStructure(Temp[ind], Xrot[ind], Yrot[ind], Zrot[ind], T)
    end

    # Set the phase. Different routines are available for that - see below.
    Phase[ind] = Compute_Phase(Phase[ind], Temp[ind], Xrot[ind], Yrot[ind], Zrot[ind], phase)

    return nothing
end

"""
    AddCylinder!(Phase, Temp, Grid::AbstractGeneralGrid; base=Tuple{3}, cap=Tuple{3}, radius=Tuple{1},
            phase = ConstantPhase(1).
            T=nothing )

Adds a cylinder with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- Phase - Phase array (consistent with Grid)
- Temp  - Temperature array (consistent with Grid)
- Grid - Grid structure (usually obtained with ReadLaMEM_InputFile)
- base - center coordinate of bottom of cylinder
- cap - center coordinate of top of cylinder
- radius - radius of the cylinder
- phase - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()` 
- T - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()` 


Example
========

Cylinder with constant phase and temperature:
```julia 
julia> Grid = ReadLaMEM_InputFile("test_files/SaltModels.dat")
LaMEM Grid: 
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> AddCylinder!(Phases,Temp,Grid, base=(-1,-1,-1.5), cap=(1,1,-0.5), radius=0.25, phase=ConstantPhase(4), T=ConstantTemp(400))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> Write_Paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview 
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"   
```
"""
function AddCylinder!(Phase, Temp, Grid::AbstractGeneralGrid;   # required input
    base=Tuple{3}, cap=Tuple{3}, radius=Tuple{1},               # center and radius of the sphere
    phase = ConstantPhase(1),                           # Sets the phase number(s) in the sphere
    T=nothing )                                         # Sets the thermal structure (various fucntions are available)
    
    # axis vector of cylinder
    axVec = cap .- base
    ax2   = (axVec[1]^2 + axVec[2]^2 + axVec[3]^2)

    # Retrieve 3D data arrays for the grid
    X,Y,Z = coordinate_grids(Grid)

    # distance between grid points and cylinder base
    dx_b  = X .- base[1]
    dy_b  = Y .- base[2]
    dz_b  = Z .- base[3]

    # find normalized parametric coordinate of a point-axis projection
    t     = (axVec[1] .* dx_b .+ axVec[2] .* dy_b .+ axVec[3] .* dz_b) ./ ax2

    # find distance vector between point and axis
    dx    = dx_b .- t.*axVec[1]
    dy    = dy_b .- t.*axVec[2]
    dz    = dz_b .- t.*axVec[3]

    # Set phase number & thermal structure in the full domain
    ind = findall((t .>= 0.0) .& (t .<= 1.0) .& ((dx.^2 + dy.^2 + dz.^2).^0.5 .<= radius))

    # Compute thermal structure accordingly. See routines below for different options
    if T != nothing
        Temp[ind] = Compute_ThermalStructure(Temp[ind], X[ind], Y[ind], Z[ind], T)
    end

    # Set the phase. Different routines are available for that - see below.
    Phase[ind] = Compute_Phase(Phase[ind], Temp[ind], X[ind], Y[ind], Z[ind], phase)

    return nothing
end

# Internal function that rotates the coordinates
function Rot3D!(X,Y,Z, StrikeAngle, DipAngle)

    # rotation matrixes
    roty = [cosd(-DipAngle) 0 sind(-DipAngle) ; 0 1 0 ; -sind(-DipAngle) 0  cosd(-DipAngle)];
    rotz = [cosd(StrikeAngle) -sind(StrikeAngle) 0 ; sind(StrikeAngle) cosd(StrikeAngle) 0 ; 0 0 1]

    for i in eachindex(X)
        CoordVec = [X[i], Y[i], Z[i]]
        CoordRot =  rotz*CoordVec; 
        CoordRot =  roty*CoordRot; 
        X[i] = CoordRot[1];
        Y[i] = CoordRot[2];
        Z[i] = CoordRot[3];
    end

    return nothing
end


abstract type AbstractThermalStructure end  


"""
    ConstantTemp(T=1000)
    
Sets a constant temperature inside the box

Parameters
===
- T : the value
"""
@with_kw_noshow mutable struct ConstantTemp <: AbstractThermalStructure
    T = 1000
end

function Compute_ThermalStructure(Temp, X, Y, Z, s::ConstantTemp)
    Temp .= s.T
    return Temp
end


"""
    LinearTemp(Ttop=0, Tbot=1000)
    
Set a linear temperature structure from top to bottom

Parameters
===
- Ttop : the value @ the top
- Tbot : the value @ the bottom

"""
@with_kw_noshow mutable struct LinearTemp <: AbstractThermalStructure
    Ttop = 0
    Tbot = 1350
end

function Compute_ThermalStructure(Temp, X, Y, Z, s::LinearTemp)
    @unpack Ttop, Tbot  = s

    dz   = Z[end]-Z[1];
    dT   = Tbot - Ttop

    Temp = abs.(Z./dz).*dT .+ Ttop
    return Temp
end

"""
    HalfspaceCoolingTemp(Tsurface=0, Tmantle=1350, Age=60, Adiabat=0)
    
Sets a halfspace temperature structure in plate

Parameters
========
- Tsurface : surface temperature [C]
- Tmantle : mantle temperature [C]
- Age : Thermal Age of plate [Myrs]
- Adiabat : Mantle Adiabat [K/km]

"""
@with_kw_noshow mutable struct HalfspaceCoolingTemp <: AbstractThermalStructure
    Tsurface = 0       # top T
    Tmantle = 1350     # bottom T
    Age  = 60          # thermal age of plate [in Myrs]
    Adiabat = 0        # Adiabatic gradient in K/km
end

function Compute_ThermalStructure(Temp, X, Y, Z, s::HalfspaceCoolingTemp)
    @unpack Tsurface, Tmantle, Age, Adiabat  = s

    kappa       =   1e-6;
    SecYear     =   3600*24*365
    dz          =   Z[end]-Z[1];
    ThermalAge  =   Age*1e6*SecYear;

    MantleAdiabaticT    =   Tmantle .+ Adiabat*abs.(Z);   # Adiabatic temperature of mantle
    
    for i in eachindex(Temp)
        Temp[i] =   (Tsurface .- Tmantle)*erfc((abs.(Z[i])*1e3)./(2*sqrt(kappa*ThermalAge))) + MantleAdiabaticT[i];
    end
    return Temp
end


"""
    SpreadingRateTemp(Tsurface=0, Tmantle=1350, Adiabat=0, MORside="left",SpreadingVel=3, AgeRidge=0, maxAge=80)
    
Sets a halfspace temperature structure within the box, combined with a spreading rate (which implies that the plate age varies)

Parameters
========
- Tsurface : surface temperature [C]
- Tmantle : mantle temperature [C]
- Adiabat : Mantle Adiabat [K/km]
- MORside : side of the box where the MOR is located ["left","right","front","back"]
- SpreadingVel : spreading velocity [cm/yr]
- AgeRidge : thermal age of the ridge [Myrs]
- maxAge : maximum thermal Age of plate [Myrs]

"""
@with_kw_noshow mutable struct SpreadingRateTemp <: AbstractThermalStructure
    Tsurface = 0       # top T
    Tmantle = 1350     # bottom T
    Adiabat = 0        # Adiabatic gradient in K/km
    MORside = "left"   # side of box where the MOR is located
    SpreadingVel = 3   # spreading velocity [cm/yr]
    AgeRidge = 0       # Age of the ridge [Myrs]
    maxAge  = 60       # maximum thermal age of plate [Myrs]
end

function Compute_ThermalStructure(Temp, X, Y, Z, s::SpreadingRateTemp)
    @unpack Tsurface, Tmantle, Adiabat, MORside, SpreadingVel, AgeRidge, maxAge  = s

    kappa       =   1e-6;
    SecYear     =   3600*24*365
    dz          =   Z[end]-Z[1];
  

    MantleAdiabaticT    =   Tmantle .+ Adiabat*abs.(Z);   # Adiabatic temperature of mantle
    
    if MORside=="left"
        Distance = X .- X[1,1,1]; 
    elseif MORside=="right"
        Distance = X[end,1,1] .- X; 
    elseif MORside=="front"
        Distance = Y .- Y[1,1,1]; 
    elseif MORside=="back"
        Distance = Y[1,end,1] .- Y; 
    else
        error("unknown side")
    end
    
    for i in eachindex(Temp)
        ThermalAge    =   abs(Distance[i]*1e3*1e2)/SpreadingVel + AgeRidge*1e6;   # Thermal age in years
        if ThermalAge>maxAge*1e6
            ThermalAge = maxAge*1e6
        end
        
        ThermalAge    =   ThermalAge*SecYear;

        Temp[i] = (Tsurface .- Tmantle)*erfc((abs.(Z[i])*1e3)./(2*sqrt(kappa*ThermalAge))) + MantleAdiabaticT[i];
    end
    return Temp
end



abstract type AbstractPhaseNumber end  


"""
    ConstantPhase(phase=1)
    
Sets a constant phase inside the box

Parameters
===
- phase : the value
"""
@with_kw_noshow mutable struct ConstantPhase <: AbstractPhaseNumber
    phase = 1
end

function Compute_Phase(Phase, Temp, X, Y, Z, s::ConstantPhase)
    Phase .= s.phase
    return Phase
end



"""
    LithosphericPhases(Layers=[10 20 30], Phases=[1 2 3 4], Tlab=nothing )
    
This allows defining a layered lithosphere. Layering is defined from the top downwards.

Parameters
===
- Layers : the thickness of the layers from top downwards
- Phases : the phases of the layers from top down. Note that this array 
- Tlab   : Temperature of the lithosphere asthenosphere boundary. If specified, the phases at locations with T>Tlab is set to Phases[end]

"""
@with_kw_noshow mutable struct LithosphericPhases <: AbstractPhaseNumber
    Layers  = [10., 20., 50.]
    Phases  = [1,   2  , 3,  4]
    Tlab    = nothing
end

function Compute_Phase(Phase, Temp, X, Y, Z, s::LithosphericPhases)
    @unpack Layers, Phases, Tlab  = s

    Phase .= Phases[end]
    Ztop  = 0
    for i=1:length(Layers)
        Zbot = Ztop-Layers[i]
        ind = findall( ( Z .>= Zbot) .&  (Z .<= Ztop) );
        Phase[ind] .= Phases[i]

        Ztop = Zbot
    end

    # set phase to mantle if requested
    if Tlab != nothing
        ind = findall(Temp .> Tlab)
        Phase[ind] .= Phases[end]
    end

    return Phase
end