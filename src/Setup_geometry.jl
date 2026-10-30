using Printf
using Parameters        # helps setting default parameters in structures
using SpecialFunctions: erfc
using GeoParams
using StaticArrays

import Base: show
# Setup_geometry
#
# These are routines that help to create input geometries, such as slabs with a given angle
#

export add_box!, add_sphere!, add_ellipsoid!, add_cylinder!, add_layer!, add_polygon!, add_plate!, add_slab!, add_stripes!, add_volcano!, add_fault!,
    make_volc_topo,
    ConstantTemp, LinearTemp, HalfspaceCoolingTemp, SpreadingRateTemp, LithosphericTemp, LinearWeightedTemperature,
    McKenzie_subducting_slab,
    ConstantPhase, LithosphericPhases,
    Trench, compute_slab_surface,
    compute_thermal_structure, compute_phase

# Extrema of the coordinates that are not NaN; NaN marks absent points, such as empty particle slots
finite_extrema(A) = (mapreduce(a -> isnan(a) ? typemax(a) : a, min, A), mapreduce(a -> isnan(a) ? typemin(a) : a, max, A))

"""
    add_stripes!(Phase, Grid::AbstractGeneralGrid;
        stripAxes       = (1,1,0),
        stripeWidth     =  0.2,
        stripeSpacing   =  1,
        Origin          =  nothing,
        StrikeAngle     =  0,
        DipAngle        =  10,
        phase           =  ConstantPhase(3),
        stripePhase     =  ConstantPhase(4),
        cell            = false)

Adds stripes to a pre-defined phase (e.g. added using add_box!)


Parameters
====

- `Phase` - Phase array (consistent with Grid)
- `Grid` -  grid structure (usually obtained with read_LaMEM_inputfile, but can also be other grid types)
- `stripAxes` - sets the axis for which we want the stripes. Default is (1,1,0) i.e. X, Y and not Z
- `stripeWidth` - width of the stripe
- `stripeSpacing` - space between two stripes
- `Origin` - the origin, used to rotate the box around. Default is the left-front-top corner
- `StrikeAngle` - strike angle
- `DipAngle` - dip angle
- `phase` - specifies the phase we want to apply stripes to
- `stripePhase` - specifies the stripe phase
- `cell` - if true, `Phase` and `Temp` are defined on centers

Example
========

Example: Box with striped phase and constant temperature & a dip angle of 10 degrees:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_box!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=ConstantTemp(1000))
julia> add_stripes!(Phases, Grid, stripAxes=(1,1,1), stripeWidth=0.2, stripeSpacing=1, Origin=nothing, StrikeAngle=0, DipAngle=10, phase=ConstantPhase(3), stripePhase=ConstantPhase(4))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```
"""
function add_stripes!(
        Phase, Grid::AbstractGeneralGrid;                # required input
        stripAxes = (1, 1, 0),                          # activate stripes along dimensions x, y and z when set to 1
        stripeWidth = 0.2,                             # full width of a stripe
        stripeSpacing = 1,                               # spacing between two stripes centers
        Origin = nothing,                         # origin
        StrikeAngle = 0,                               # strike
        DipAngle = 0,                               # dip angle
        phase = ConstantPhase(3),                # phase to be striped
        stripePhase = ConstantPhase(4),                # stripe phase
        cell = false
    )                          # if true, Phase and Temp are defined on cell centers

    # warnings
    if stripeWidth >= stripeSpacing / 2.0
        print("WARNING: stripeWidth should be strictly < stripeSpacing/2.0, otherwise phase is overwritten by the stripePhase\n")
    elseif sum(stripAxes .== 0) == 3
        print("WARNING: at least one axis should be set to 1 e.g. stripAxes = (1,0,0), otherwise no stripes will be added\n")
    end

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # sets origin
    if isnothing(Origin)
        Origin = last.(finite_extrema.((X, Y, Z)))  # upper-left corner
    end

    # vec, so that a 2D Phase array works with a grid of size (nx, 1, nz)
    length(Phase) == length(X) || throw(DimensionMismatch("Phase has $(length(Phase)) entries but the grid has $(length(X)) points"))
    p = (; r = rotation(Origin, StrikeAngle, DipAngle), stripeSpacing, stripeWidth, stripAxes, phase = phase.phase, stripePhase = stripePhase.phase)
    P = vec(Phase)   # one wrapper on both sides, so that the broadcast needs no copy
    P .= stripe_phase.(Ref(p), P, vec(X), vec(Y), vec(Z))

    return nothing
end

in_stripe(p, x) = abs(rem(x, p.stripeSpacing, RoundNearest)) <= p.stripeWidth / 2.0

function stripe_phase(p, ph, x, y, z)
    ph == p.phase || return ph
    xr, yr, zr = rotate(p.r, x, y, z)
    striped = (p.stripAxes[1] == 1 && in_stripe(p, xr)) || (p.stripAxes[2] == 1 && in_stripe(p, yr)) || (p.stripAxes[3] == 1 && in_stripe(p, zr))
    return striped ? p.stripePhase : ph
end


"""
    add_box!(Phase, Temp, Grid::AbstractGeneralGrid; xlim::Tuple = (20,100), [ylim::Tuple = (1,10)], zlim::Tuple = (10,80),
            Origin=nothing, StrikeAngle=0, DipAngle=0,
            phase = ConstantPhase(1),
            T=nothing,
            segments=nothing,
            cell=false )

Adds a box with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid` -  grid structure (can be any of the grid types in `GMG`)
- `xlim` -  left/right coordinates of box
- `ylim` -  front/back coordinates of box [optional; if not specified we use the whole box]
- `zlim` -  bottom/top coordinates of box
- `Origin` - the origin, used to rotate the box around. Default is the left-front-top corner
- `StrikeAngle` - strike angle of slab
- `DipAngle` - dip angle of slab
- `phase` - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T` - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`,`LithosphericTemp()`
- `segments` - optional parameter to define multiple ridge segments within the box
- `cell` - if true, `Phase` and `Temp` are defined on centers

Examples
========

Example 1) Box with constant phase and temperature & a dip angle of 10 degrees:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_box!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=ConstantTemp(1000))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```

Example 2) Box with halfspace cooling profile
```julia-repl
julia> Grid = CartData(xyz_grid(-1000:10:1000,0,-660:10:0))
julia> Phases = zeros(Int32,   size(Grid));
julia> Temp   = zeros(Float64, size(Grid));
julia> add_box!(Phases,Temp,Grid, xlim=(0,500), zlim=(-50,0), phase=ConstantPhase(3), DipAngle=10, T=HalfspaceCoolingTemp(Age=30))
julia> Grid = addfield(Grid, (;Phases,Temp));       # Add to Cartesian model
julia> write_paraview(Grid,"LaMEM_ModelSetup")  # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```

Example 3) Box with ridge thermal structure
```julia-repl
julia> Grid = CartData(xyz_grid(-1000:10:1000, -1000:10:1000, -660:5:0))
julia> Phases = fill(2, size(Grid));
julia> Temp   = fill(1350.0, size(Grid));
julia> segments = [((-500.0, -1000.0), (-500.0, 0.0)),
                    ((-250.0, 0.0), (-250.0, 200.0)),
                    ((-750.0, 200.0), (-750.0, 1000.0))];
julia> lith = LithosphericPhases(Layers=[15 55], Phases=[1 2], Tlab=1250);
julia> add_box!(Phases, Temp, Grid; xlim=(-1000.0, 0.0), ylim=(-500.0, 500.0),
                zlim=(-80.0, 0.0), phase=lith,
                T=SpreadingRateTemp(SpreadingVel=3), segments=segments)
julia> Grid = addfield(Grid, (; Phases, Temp));       # Add to Cartesian model
julia> write_paraview(Grid, "Ridge_Thermal_Structure")  # Save model to Paraview
1-element Vector{String}:
 "Ridge_Thermal_Structure.vts"
"""
function add_box!(
        Phase, Temp, Grid::AbstractGeneralGrid;       # required input
        xlim::Tuple = (20, 100), ylim = nothing, zlim::Tuple = (10, 80),     # limits of the box
        Origin = nothing, StrikeAngle = 0, DipAngle = 0,      # origin & dip/strike
        phase = ConstantPhase(1),                       # Sets the phase number(s) in the box
        T = nothing,                              # Sets the thermal structure (various functions are available)
        segments = nothing,                       # Allows defining multiple ridge segments
        cell = false
    )                            # if true, Phase and Temp are defined on cell centers


    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # ensure that the input arrays have the correct size
    #@assert size(X) == size(Phase) == size(Temp)

    # Limits of block
    if isnothing(ylim)
        ylim = finite_extrema(Y)
    end

    # McKenzie temperature is computed relative to [xmin,ymin,zmax] of the box
    if Origin !== nothing && isa(T, McKenzie_subducting_slab)
        if Origin[1] != xlim[1] || Origin[2] != ylim[1] || Origin[3] != zlim[2]
            error("Origin is not set up correctly for McKenzie_subducting_slab. Leave Origin blank or set Origin = [xmin,ymin,zmax] of the box")
        end
    end

    if isnothing(Origin)
        Origin = (xlim[1], ylim[1], zlim[2])  # upper-left corner
    end

    # Set phase number & thermal structure in the full domain
    ztop = maximum(zlim) - Origin[3]
    zbot = minimum(zlim) - Origin[3]
    xmin, xmax = extrema(xlim) .- Origin[1]
    ymin, ymax = extrema(ylim) .- Origin[2]
    p = (; r = rotation(Origin, StrikeAngle, DipAngle), xmin, xmax, ymin, ymax, zbot, ztop)
    frame = (; xmin, xmax, ymin, ymax, zbot, ztop)

    if isnothing(segments)
        set_shape!(Phase, Temp, X, Y, Z, box_coordinates, p, frame, phase, T)
    else
        set_shape!(Phase, Temp, X, Y, Z, ridge_coordinates, p, frame, phase, T; T_args = (segments,))
    end

    return nothing
end

function box_coordinates(p, x, y, z)
    xyz = rotate(p.r, x, y, z)
    xr, yr, zr = xyz
    return p.xmin <= xr <= p.xmax && p.ymin <= yr <= p.ymax && p.zbot <= zr <= p.ztop ? xyz : nothing
end

# ridge segments are given in absolute x, y coordinates
function ridge_coordinates(p, x, y, z)
    xyz = box_coordinates(p, x, y, z)
    return isnothing(xyz) ? nothing : (x, y, xyz[3])
end

"""
    add_box!(Phase, Temp, Grid::AbstractGeneralGrid,
    bounds::AbstractVector{<:AbstractVector{<:Real}};
    Origin = nothing, StrikeAngle = 0, DipAngle = 0,
    phase = ConstantPhase(1),
    T = nothing,
    segments = nothing,
    cell = false     )
Add box function but getting bounds as a vector bounds in a way of [[xmin,xmax],[ymin,ymax],[zmin,zmax]]. If bounds is empty return nothing 
"""
function add_box!(
        Phase, Temp, Grid::AbstractGeneralGrid,           # required input
        bounds::Union{Vector{Any}, AbstractVector{<:AbstractVector{<:Real}}};    # limits of the box
        Origin = nothing, StrikeAngle = 0, DipAngle = 0,     # origin & dip/strike
        phase = ConstantPhase(1),                            # Sets the phase number(s) in the box
        T = nothing,                                         # Sets the thermal structure (various functions are available)
        segments = nothing,                                  # Allows defining multiple ridge segments
        cell = false                                         # if true, Phase and Temp are defined on cell centers
    )

    isempty(bounds) && return nothing

    return add_box!(
        Phase, Temp, Grid;
        xlim = Tuple(bounds[1]), ylim = Tuple(bounds[2]), zlim = Tuple(bounds[3]),
        Origin, StrikeAngle, DipAngle, phase, T, segments, cell
    )
end

"""
    add_layer!(Phase, Temp, Grid::AbstractGeneralGrid; xlim::Tuple = (1,100), [ylim::Tuple = (0,20)], zlim::Tuple = (0,-100),
            phase = ConstantPhase(1),
            T=nothing, cell=false )


Adds a layer with phase & temperature structure to a 3D model setup. The most common use would be to add a lithospheric layer to a model setup.
This simplifies creating model geometries in geodynamic models


Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid` -  grid structure (usually obtained with read_LaMEM_inputfile, but can also be other grid types)
- `xlim` -  left/right coordinates of box
- `ylim` -  front/back coordinates of box
- `zlim` -  bottom/top coordinates of box
- `phase` - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T` - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`


Examples
========

Example 1) Layer with constant phase and temperature
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_layer!(Phases,Temp,Grid, zlim=(-50,0), phase=ConstantPhase(3), T=ConstantTemp(1000))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```

Example 2) Box with halfspace cooling profile
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_layer!(Phases,Temp,Grid, zlim=(-50,0), phase=ConstantPhase(3), T=HalfspaceCoolingTemp())
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```
"""
function add_layer!(
        Phase, Temp, Grid::AbstractGeneralGrid;     # required input
        xlim = nothing, ylim = nothing, zlim = nothing,       # limits of the layer
        phase = ConstantPhase(1),                       # Sets the phase number(s) in the box
        T = nothing,                                      # Sets the thermal structure (various functions are available)
        cell = false
    )                                 # if true, Phase and Temp are defined on cell centers

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Limits of block
    if all(isnothing, (xlim, ylim, zlim))
        error("You need to specify at least one of the limits (xlim, ylim, zlim)")
    end

    if isnothing(xlim)
        xlim = finite_extrema(X)
    end
    if isnothing(ylim)
        ylim = finite_extrema(Y)
    end
    if isnothing(zlim)
        zlim = finite_extrema(Z)
    end

    # Set phase number & thermal structure in the full domain
    (x1, x2), (y1, y2), (z1, z2) = xlim, ylim, zlim
    p = (; x1, x2, y1, y2, z1, z2)
    frame = (; xmin = x1, xmax = x2, ymin = y1, ymax = y2, zbot = minimum(zlim), ztop = maximum(zlim))
    set_shape!(Phase, Temp, X, Y, Z, layer_coordinates, p, frame, phase, T)

    return nothing
end

layer_coordinates(p, x, y, z) = p.x1 <= x <= p.x2 && p.y1 <= y <= p.y2 && p.z1 <= z <= p.z2 ? (x, y, z) : nothing


"""
    add_sphere!(Phase, Temp, Grid::AbstractGeneralGrid; cen::Tuple = (0,0,-1), radius::Number,
            phase = ConstantPhase(1).
            T=nothing, cell=false )


Adds a sphere with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid` - LaMEM grid structure (usually obtained with read_LaMEM_inputfile)
- `cen` - center coordinates of sphere
- `radius` - radius of sphere
- `phase` - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T` - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`
- `cell` - if true, `Phase` and `Temp` are defined on cell centers


Example
========

Sphere with constant phase and temperature:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_sphere!(Phases,Temp,Grid, cen=(0,0,-1), radius=0.5, phase=ConstantPhase(2), T=ConstantTemp(800))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```
"""
function add_sphere!(
        Phase, Temp, Grid::AbstractGeneralGrid;    # required input
        cen::Tuple = (0, 0, -1), radius::Number,              # center and radius of the sphere
        phase = ConstantPhase(1),                                   # Sets the phase number(s) in the sphere
        T = nothing, cell = false
    )                         # Sets the thermal structure (various functions are available)

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Set phase number & thermal structure in the full domain
    p = (; cen, radius)
    frame = (; xmin = cen[1] - radius, xmax = cen[1] + radius, ymin = cen[2] - radius, ymax = cen[2] + radius, zbot = cen[3] - radius, ztop = cen[3] + radius)
    set_shape!(Phase, Temp, X, Y, Z, sphere_coordinates, p, frame, phase, T)

    return nothing
end

sphere_coordinates(p, x, y, z) = sqrt((x - p.cen[1])^2 + (y - p.cen[2])^2 + (z - p.cen[3])^2) < p.radius ? (x, y, z) : nothing

"""
    add_ellipsoid!(Phase, Temp, Grid::AbstractGeneralGrid; cen::Tuple = (-1,-1,-1), axes::Tuple = (0.2,0.1,0.5),
            Origin=nothing, StrikeAngle=0, DipAngle=0,
            phase = ConstantPhase(1).
            T=nothing, cell=false )

Adds an Ellipsoid with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid` - LaMEM grid structure (usually obtained with read_LaMEM_inputfile)
- `cen` - center coordinates of sphere
- `axes` - semi-axes of ellipsoid in X,Y,Z
- `Origin` - the origin, used to rotate the box around. Default is the left-front-top corner
- `StrikeAngle` - strike angle of slab
- `DipAngle` - dip angle of slab
- `phase` - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T` - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`
- `cell` - if true, `Phase` and `Temp` are defined on cell centers

Example
========

Ellipsoid with constant phase and temperature, rotated 90 degrees and tilted by 45 degrees:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_ellipsoid!(Phases,Temp,Grid, cen=(-1,-1,-1), axes=(0.2,0.1,0.5), StrikeAngle=90, DipAngle=45, phase=ConstantPhase(3), T=ConstantTemp(600))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```
"""
function add_ellipsoid!(
        Phase, Temp, Grid::AbstractGeneralGrid;     # required input
        cen::Tuple = (-1, -1, -1), axes::Tuple = (0.2, 0.1, 0.5),   # center and semi-axes of the ellpsoid
        Origin = nothing, StrikeAngle = 0, DipAngle = 0,                      # origin & dip/strike
        phase = ConstantPhase(1),                                       # Sets the phase number(s) in the box
        T = nothing, cell = false
    )                             # Sets the thermal structure (various functions are available)

    if isnothing(Origin)
        Origin = cen  # center
    end

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Set phase number & thermal structure in the full domain
    x2 = axes[1]^2
    y2 = axes[2]^2
    z2 = axes[3]^2
    cenRot = cen .- Origin
    p = (; r = rotation(Origin, StrikeAngle, DipAngle), cenRot, x2, y2, z2)
    frame = (;
        xmin = cenRot[1] - axes[1], xmax = cenRot[1] + axes[1], ymin = cenRot[2] - axes[2], ymax = cenRot[2] + axes[2],
        zbot = cenRot[3] - axes[3], ztop = cenRot[3] + axes[3],  # see #135
    )
    set_shape!(Phase, Temp, X, Y, Z, ellipsoid_coordinates, p, frame, phase, T)

    return nothing
end

@inline function ellipsoid_coordinates(p, x, y, z)
    xyz = rotate(p.r, x, y, z)
    xr, yr, zr = xyz
    return sqrt((xr - p.cenRot[1])^2 / p.x2 + (yr - p.cenRot[2])^2 / p.y2 + (zr - p.cenRot[3])^2 / p.z2) <= 1 ? xyz : nothing
end

"""
    add_cylinder!(Phase, Temp, Grid::AbstractGeneralGrid; base::Tuple = (-1,-1,-1.5), cap::Tuple = (-1,-1,-0.5), radius::Number,
            phase = ConstantPhase(1),
            T=nothing, cell=false )


Adds a cylinder with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models


Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid` - Grid structure (usually obtained with `read_LaMEM_inputfile`)
- `base` - center coordinate of bottom of cylinder
- `cap` - center coordinate of top of cylinder
- `radius` - radius of the cylinder
- `phase` - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T` - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`
- `cell` - if true, `Phase` and `Temp` are defined on cell centers


Example
========

Cylinder with constant phase and temperature:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_cylinder!(Phases,Temp,Grid, base=(-1,-1,-1.5), cap=(1,1,-0.5), radius=0.25, phase=ConstantPhase(4), T=ConstantTemp(400))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```
"""
function add_cylinder!(
        Phase, Temp, Grid::AbstractGeneralGrid;  # required input
        base::Tuple = (-1, -1, -1.5), cap::Tuple = (-1, -1, -0.5), radius::Number,    # center and radius of the sphere
        phase = ConstantPhase(1),                                   # Sets the phase number(s) in the sphere
        T = nothing, cell = false
    )                             # Sets the thermal structure (various functions are available)

    # axis vector of cylinder
    axVec = cap .- base
    ax2 = (axVec[1]^2 + axVec[2]^2 + axVec[3]^2)

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Set phase number & thermal structure in the full domain
    p = (; base, axVec, ax2, radius)

    # Vertical extent of the cylinder, including the tilt of its circular end caps; see #135
    zhalf = radius * sqrt(1 - axVec[3]^2 / ax2)
    frame = (;
        xmin = min(base[1], cap[1]) - radius, xmax = max(base[1], cap[1]) + radius,
        ymin = min(base[2], cap[2]) - radius, ymax = max(base[2], cap[2]) + radius,
        zbot = min(base[3], cap[3]) - zhalf, ztop = max(base[3], cap[3]) + zhalf,
    )
    set_shape!(Phase, Temp, X, Y, Z, cylinder_coordinates, p, frame, phase, T)

    return nothing
end

function in_cylinder(p, x, y, z)
    # distance between the point and the cylinder base
    dx_b, dy_b, dz_b = x - p.base[1], y - p.base[2], z - p.base[3]
    # normalized parametric coordinate of the point-axis projection
    t = (p.axVec[1] * dx_b + p.axVec[2] * dy_b + p.axVec[3] * dz_b) / p.ax2
    # distance vector between the point and the axis
    dx, dy, dz = dx_b - t * p.axVec[1], dy_b - t * p.axVec[2], dz_b - t * p.axVec[3]
    return 0 <= t <= 1 && sqrt(dx^2 + dy^2 + dz^2) <= p.radius
end

cylinder_coordinates(p, x, y, z) = in_cylinder(p, x, y, z) ? (x, y, z) : nothing


# Sets the thermal structure `T` (unless it is `nothing`) and then `phase` at the points (X, Y, Z) that lie
# inside a shape; `Phase` and `Temp` have one entry per point, in the same order.
# `shape(p, x, y, z)` returns the coordinates of a point in the frame of the shape, or `nothing` for a point
# outside it; `p` holds the parameters of the shape. `frame` is the extent of the shape in that frame (xmin, xmax, ymin, ymax, zbot, ztop),
# which thermal structures use as reference. `T_args` are extra arguments of `T` (ridge segments).
function set_shape!(Phase, Temp, X, Y, Z, shape, p, frame, phase, T; T_args = ())
    length(Phase) == length(X) || throw(DimensionMismatch("Phase has $(length(Phase)) entries but the grid has $(length(X)) points"))
    length(Temp) == length(X) || throw(DimensionMismatch("Temp has $(length(Temp)) entries but the grid has $(length(X)) points"))
    if needs_points(T)
        set_shape_gathered!(Phase, Temp, X, Y, Z, shape, p, frame, phase, T, T_args)
    else
        cT = isnothing(T) ? nothing : prepare_temperature(T, T_args..., frame, nothing)
        set_points!(Phase, Temp, X, Y, Z, shape, p, phase, T, cT)
    end
    return nothing
end

function set_points!(Phase, Temp, X, Y, Z, shape, p, phase, T, cT)
    for (iP, iT, i) in zip(eachindex(Phase), eachindex(Temp), eachindex(X, Y, Z))
        xyz = shape(p, X[i], Y[i], Z[i])
        isnothing(xyz) && continue
        xl, yl, zl = xyz
        if !isnothing(T)
            Temp[iT] = point_temperature(T, cT, xl, yl, zl)
        end
        Phase[iP] = point_phase(phase, Phase[iP], Temp[iT], xl, yl, zl)
    end
    return nothing
end

# For thermal structures that need the points of the region: the phase is set first, as the thermal
# structure depends on it, then the temperature, then the phase again with the new temperature.
function set_shape_gathered!(Phase, Temp, X, Y, Z, shape, p, frame, phase, T, T_args)
    xyz = shape.(Ref(p), X, Y, Z)
    ind = findall(!isnothing, vec(xyz))
    isempty(ind) && return nothing
    xl, yl, zl = (getindex.(xyz[ind], k) for k in 1:3)
    P, Tp = vec(Phase)[ind], vec(Temp)[ind]
    P = compute_phase(P, Tp, xl, yl, zl, phase)
    Tp = point_temperatures!(Tp, xl, yl, zl, T, prepare_temperature(T, T_args..., frame, (xl, yl, zl, P)))
    vec(Temp)[ind] = Tp
    vec(Phase)[ind] = compute_phase(P, Tp, xl, yl, zl, phase)
    return nothing
end


"""
        add_polygon!(Phase, Temp, Grid::AbstractGeneralGrid; xlim=(), ylim=(0.0,0.8), zlim=(), phase = ConstantPhase(1), T=nothing, cell=false )


Adds a polygon with phase & temperature structure to a 3D model setup.  This simplifies creating model geometries in geodynamic models

Parameters
====
- `Phase` - Phase array (consistent with Grid)
- `Temp`  - Temperature array (consistent with Grid)
- `Grid`  - Grid structure (usually obtained with read_LaMEM_inputfile)
- `xlim`  - `x`-coordinate of the polygon points, same ordering as zlim, number of points unlimited
- `ylim`  - `y`-coordinate, limitation in length possible (two values (start and stop))
- `zlim`  - `z`-coordinate of the polygon points, same ordering as xlim, number of points unlimited
- `phase` - specifies the phase of the box. See `ConstantPhase()`
- `T`     - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`
- `cell`  - if true, `Phase` and `Temp` are defined on cell centers

Example
========

Polygon with constant phase and temperature:

```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X));
julia> Temp   = zeros(Float64, size(Grid.X));
julia> add_polygon!(Phase, Temp, Cart; xlim=(0,0, 1.6, 2.0),ylim=(0,0.8), zlim=(0,-1,-2,0), phase = ConstantPhase(8), T=ConstantTemp(30))
julia> Model3D = ParaviewData(Grid, (Phases=Phases,Temp=Temp)); # Create Cartesian model
julia> write_paraview(Model3D,"LaMEM_ModelSetup")           # Save model to paraview
1-element Vector{String}:
 "LaMEM_ModelSetup.vts"
```

"""
function add_polygon!(
        Phase, Temp, Grid::AbstractGeneralGrid;   # required input
        xlim = (), ylim = (0.0, 0.8), zlim = (),           # limits of the box
        phase = ConstantPhase(1),                                   # Sets the phase number(s) in the box
        T = nothing, cell = false
    )                             # Sets the thermal structure (various functions are available)


    xlim_ = Float64.(collect(xlim))
    ylim_ = Float64.(collect(ylim))
    zlim_ = Float64.(collect(zlim))


    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # points within the polygon in the x-z plane and between the y limits
    p = (; xlim_, ylim_, zlim_)
    (xmin, xmax), (zbot, ztop) = extrema(xlim_), extrema(zlim_)
    frame = (; xmin, xmax, ymin = ylim_[1], ymax = ylim_[2], zbot, ztop)
    set_shape!(Phase, Temp, X, Y, Z, polygon_coordinates, p, frame, phase, T)

    return nothing
end

in_polygon(p, x, y, z) = inpoly_either(p.xlim_, p.zlim_, x, z) && p.ylim_[1] <= y <= p.ylim_[2]

polygon_coordinates(p, x, y, z) = in_polygon(p, x, y, z) ? (x, y, z) : nothing

"""
    add_polygon!(Phase, Temp, Grid::AbstractGeneralGrid,
    bounds::AbstractVector{<:AbstractVector{<:Real}};
    phase = ConstantPhase(1),
    T = nothing,
    cell = false     )
Add polygon function but getting bounds as a vector bounds in a way of [[xmin,xmax],[ymin,ymax],[zmin,zmax]]. If bounds is empty return nothing 
"""
function add_polygon!(
        Phase, Temp, Grid::AbstractGeneralGrid,           # required input
        bounds::Union{Vector{Any}, AbstractVector{<:AbstractVector{<:Real}}};    # limits of the box
        phase = ConstantPhase(1),                            # Sets the phase number(s) in the box
        T = nothing,                                         # Sets the thermal structure (various functions are available)
        cell = false                                         # if true, Phase and Temp are defined on cell centers
    )

    return if !isempty(bounds)
        xlim = (Tuple(bounds[1]))
        ylim = (Tuple(bounds[2]))
        zlim = (Tuple(bounds[3]))

        add_polygon!(
            Phase, Temp, Grid;       # required input
            xlim = xlim, ylim = ylim, zlim = zlim,     # limits of the box
            phase = phase,                              # Sets the phase number(s) in the box
            T = T,                                  # Sets the thermal structure (various functions are available)
            cell = cell
        )
    end

end

"""
        add_plate!(Phase, Temp, Grid::AbstractGeneralGrid; xlim=(), ylim=(), zlim::Tuple = (0.0,0.8), phase = ConstantPhase(1), T=nothing, segments=nothing, cell=false )
Adds a tectonic plate with phase and temperature structure to a 3D model setup.
This function enables the definition of tectonic plates in the xy plane and projects them along the z-axis, providing a flexible approach to model complex plate geometries.
Parameters
==========
- `Phase`  - Phase array (consistent with Grid)
- `Temp`   - Temperature array (consistent with Grid)
- `Grid`   - Grid structure (usually obtained with read_LaMEM_inputfile)
- `xlim`   - `x`-coordinate of the polygon points, same ordering as ylim, number of points unlimited
- `ylim`   - `y`-coordinate of the polygon points, same ordering as xlim, number of points unlimited
- `zlim`   - `z`-coordinate range for projecting the polygon (start and stop, two values)
- `phase`  - Specifies the phase of the plate. See `ConstantPhase()`
- `T`      - Specifies the temperature of the plate. See `ConstantTemp()`, `LinearTemp()`, `HalfspaceCoolingTemp()`, `SpreadingRateTemp()`
- `segments` - Optional. Allows for thermal segmentation within the polygon. Useful for ridge systems or complex thermal structures.
- `cell`   - If true, `Phase` and `Temp` are defined on cell centers
Example
========
Tectonic plate in the xy plane with phase and temperature structure:
```julia-repl
julia> Grid = CartData(xyz_grid(x, y, z))
Grid:
  nel         : (512, 512, 128)
  marker/cell : (1, 1, 1)
  markers     : (512, 512, 128)
  x           ϵ [-1000.0 : 0.0]
  y           ϵ [-1000.0 : 1000.0]
  z           ϵ [-660.0 : 0.0]
julia> Phases = zeros(Int32,   size(Grid.X))
julia> Temp   = zeros(Float64, size(Grid.X))
julia> segments = [
           ((-500.0, -1000.0), (-500.0, 0.0)),  # Segment 1
           ((-250.0, 0.0), (-250.0, 200.0)),    # Segment 2
           ((-750.0, 200.0), (-750.0, 1000.0))  # Segment 3
       ]
julia> lith = LithosphericPhases(Layers=[15 55], Phases=[1 2], Tlab=1250)
julia> add_plate!(Phases, Temp, Grid;
           xlim=(-1000.0, -750.0, -250.0, 0.0, -250.0, -750.0),
           ylim=(0.0, 500.0, 500.0, 0.0, -500.0, -500.0),
           zlim=(-150.0, 0.0),
           phase=lith,
           T=SpreadingRateTemp(SpreadingVel=3),
           segments=segments)
julia> Grid = addfield(Grid, (; Phases, Temp))  # Add fields
julia> write_paraview(Grid, "Plate")  # Save model to Paraview
1-element Vector{String}:
 "Plate.vts"
"""

function add_plate!(
        Phase, Temp, Grid::AbstractGeneralGrid;
        xlim = (), ylim = (), zlim::Tuple = (0.0, 0.8),
        phase = ConstantPhase(1),
        T = nothing, segments = nothing, cell = false
    )

    xlim_ = collect(xlim)
    ylim_ = collect(ylim)
    zlim_ = collect(zlim)

    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # points within the polygon in the x-y plane and between the z limits
    p = (; xlim_, ylim_, zlim_)
    (xmin, xmax), (ymin, ymax), (zbot, ztop) = extrema(xlim_), extrema(ylim_), extrema(zlim_)
    frame = (; xmin, xmax, ymin, ymax, zbot, ztop)
    set_shape!(Phase, Temp, X, Y, Z, plate_coordinates, p, frame, phase, T; T_args = isnothing(segments) ? () : (segments,))

    return nothing
end

in_plate(p, x, y, z) = p.zlim_[1] <= z <= p.zlim_[2] && inpoly_either(p.xlim_, p.ylim_, x, y)

plate_coordinates(p, x, y, z) = in_plate(p, x, y, z) ? (x, y, z) : nothing

"""
    add_plate!(Phase, Temp, Grid::AbstractGeneralGrid,
    bounds::AbstractVector{<:AbstractVector{<:Real}};
    phase = ConstantPhase(1),
    T = nothing,
    segments = nothing,
    cell = false     )
Add plate function but getting bounds as a vector bounds in a way of [[xmin,xmax],[ymin,ymax],[zmin,zmax]]. If bounds is empty return nothing 
"""
function add_plate!(
        Phase, Temp, Grid::AbstractGeneralGrid,           # required input
        bounds::Union{Vector{Any}, AbstractVector{<:AbstractVector{<:Real}}};    # limits of the box
        phase = ConstantPhase(1),                            # Sets the phase number(s) in the box
        T = nothing,                                         # Sets the thermal structure (various functions are available)
        segments = nothing,                                  # Allows defining multiple ridge segments
        cell = false                                         # if true, Phase and Temp are defined on cell centers
    )

    isempty(bounds) && return nothing

    return add_plate!(
        Phase, Temp, Grid;
        xlim = Tuple(bounds[1]), ylim = Tuple(bounds[2]), zlim = Tuple(bounds[3]),
        phase, T, segments, cell
    )
end

"""
    xrot, yrot, zrot = Rot3D(X::Number,Y::Number,Z::Number, cosStrikeAngle, sindStrikeAngle, cosDipAngle, sinDipAngle)

    Perform rotation for a point in 3D space
"""
# Coordinates relative to the origin `O`, rotated by `StrikeAngle` and `DipAngle` (in degrees)
function rotation(O, StrikeAngle, DipAngle)
    sinStrike, cosStrike = sincosd(StrikeAngle)
    sinDip, cosDip = sincosd(-DipAngle)   # Rot3D takes the dip angle with the opposite sign
    return (; O, cosStrike, sinStrike, cosDip, sinDip)
end

# Applies the rotation `r` returned by `rotation`
rotate(r, x, y, z) = Rot3D(x - r.O[1], y - r.O[2], z - r.O[3], r.cosStrike, r.sinStrike, r.cosDip, r.sinDip)

function Rot3D(X::_T, Y::_T, Z::_T, cosStrikeAngle::_T, sindStrikeAngle::_T, cosDipAngle::_T, sinDipAngle::_T) where {_T <: Number}

    # rotation matrixes
    #roty = [cosd(-DipAngle) 0 sind(-DipAngle) ; 0 1 0 ; -sind(-DipAngle) 0  cosd(-DipAngle)];
    roty = @SMatrix [cosDipAngle 0 sinDipAngle ; 0 1 0 ; -sinDipAngle 0  cosDipAngle]        # note that dip-angle is changed from before!
    rotz = @SMatrix [cosStrikeAngle -sindStrikeAngle 0 ; sindStrikeAngle cosStrikeAngle 0 ; 0 0 1]

    CoordVec = @SVector [X, Y, Z]
    CoordRot = rotz * CoordVec
    CoordRot = roty * CoordRot

    return CoordRot[1], CoordRot[2], CoordRot[3]
end

"""
add_volcano!(
    Phases, Temp, Grid::CartData;
    volcanic_phase,
    center,
    height,
    radius,
    crater,
    base,
    background,
    T,
)

Adds a volcano topography (cones and truncated cones)

Parameters
====
- Phases - Phase array (consistent with Grid)
- Temp - Temperature array (consistent with Grid)
- Grid - CartData

Optional Parameters
====
- volcanic_phase - phase number of the volcano,
- center - x- and -coordinates of center of volcano
- height - height of volcano
- radius - radius of volcano
- T - temperature structure, applied to every point with a phase > 0 (all rock, not only the volcano; phase 0 is taken as air), with depth measured from the new topography, so that isotherms follow the volcano surface
- crater - this will create a truncated cone and the option defines the radius of the flat top
- base - this sets the flat topography around the volcano
- background - this allows loading in a topography and only adding the volcano on top (also allows stacking of several cones to get a volcano with different slopes)
"""
function add_volcano!(
        Phases,
        Temp,
        Grid::AbstractGeneralGrid;
        volcanic_phase = 1,
        center = (0, 0, 0),
        height = 0.0,
        radius = 0.0,
        crater = 0.0,
        base = 0.0,
        background = nothing,
        T = HalfspaceCoolingTemp(Age = 0)
    )
    if isnothing(background) && !needs_points(T)
        X, Y, Z = coordinate_grids(Grid)
        length(Phases) == length(X) || throw(DimensionMismatch("Phases has $(length(Phases)) entries but the grid has $(length(X)) points"))
        cone = (; center, height, radius, crater, base)
        return add_volcano_points!(vec(Phases), vec(Temp), vec(X), vec(Y), vec(Z), cone, volcanic_phase, T)
    end
    Grid isa CartData || error("add_volcano! with a background topography or a $(nameof(typeof(T))) thermal structure requires a CartData grid")

    H = make_volc_topo(
        Grid;
        center = center,
        height = height,
        radius = radius,
        crater = crater,
        base = base,
        background = background
    )

    # Broadcasts rather than loops: Grid.z.val is not inferable, so a loop here would dispatch per point
    Z = Grid.z.val
    P, Tp = reshape(Phases, size(Z)), reshape(Temp, size(Z))   # also for 2D arrays on a grid of size (nx, 1, nz)
    depth = max.(H .- Z, 0)
    P[(Z .< H) .& (Z .≥ base)] .= volcanic_phase
    ind = findall(P .> 0)

    if !isempty(ind) && !isnothing(T)
        x, y, d = Grid.x.val[ind], Grid.y.val[ind], depth[ind]
        cT = prepare_temperature(T, volcano_frame(center, radius, maximum(d)), (x, y, d, P[ind]))
        Tp[ind] = point_temperatures!(Tp[ind], x, y, d, T, cT)
    end

    return nothing
end

"""
make_volc_topo(Grid::LaMEM_grid; center::Array{Float64, 1}, height::Float64, radius::Float64, crater::Float64,
            base=0.0m, background=nothing)

Creates a generic volcano topography (cones and truncated cones)


Parameters
====
- Grid - LaMEM grid (created by read_LaMEM_inputfile)
- center - x- and -coordinates of center of volcano
- height - height of volcano
- radius - radius of volcano

Optional Parameters
====
- crater - this will create a truncated cone and the option defines the radius of the flat top
- base - this sets the flat topography around the volcano
- background - this allows loading in a topography and only adding the volcano on top (also allows stacking of several cones to get a volcano with different slopes)


Example
========

Cylinder with constant phase and temperature:
```julia-repl
julia> Grid = read_LaMEM_inputfile("test_files/SaltModels.dat")
LaMEM Grid:
  nel         : (32, 32, 32)
  marker/cell : (3, 3, 3)
  markers     : (96, 96, 96)
  x           ϵ [-3.0 : 3.0]
  y           ϵ [-2.0 : 2.0]
  z           ϵ [-2.0 : 0.0]
julia> Topo = make_volc_topo(Grid, center=[0.0,0.0], height=0.4, radius=1.5, crater=0.5, base=0.1)
CartData
    size    : (33, 33, 1)
    x       ϵ [ -3.0 : 3.0]
    y       ϵ [ -2.0 : 2.0]
    z       ϵ [ 0.1 : 0.4]
    fields  : (:Topography,)
  attributes: ["note"]
julia> Topo = make_volc_topo(Grid, center=[0.0,0.0], height=0.8, radius=0.5, crater=0.0, base=0.4, background=Topo.fields.Topography)
CartData
    size    : (33, 33, 1)
    x       ϵ [ -3.0 : 3.0]
    y       ϵ [ -2.0 : 2.0]
    z       ϵ [ 0.1 : 0.8]
    fields  : (:Topography,)
  attributes: ["note"]
julia> write_paraview(Topo,"VolcanoTopo")           # Save topography to paraview
Saved file: VolcanoTopo.vts
```
"""
function make_volc_topo(
        Grid::LaMEM_grid;
        center::Array{Float64, 1},
        height::Float64,
        radius::Float64,
        crater = 0.0,
        base = 0.0,
        background = nothing
    )

    # create nondimensionalization object
    CharUnits = SI_units(length = 1000m)

    # get node grid
    X = Grid.Xn[:, :, 1]
    Y = Grid.Yn[:, :, 1]
    nx = size(X, 1)
    ny = size(X, 2)

    # compute radial distance to volcano center
    DX = X .- center[1]
    DY = Y .- center[2]
    RD = (DX .^ 2 .+ DY .^ 2) .^ 0.5

    # get radial distance from crater rim
    RD .-= crater

    # find position relative to crater rim
    dr = radius - crater
    pos = (-RD ./ dr .+ 1)

    ## assign topography
    H = zeros(Float64, (nx, ny))
    # check if there is a background supplied
    if background === nothing
        H .= base
    else
        background = nondimensionalize(background, CharUnits)
        if size(background) == size(X)
            H .= background
        elseif size(background) == size(reshape(X, nx, ny, 1))
            H .= background[:, :, 1]
        else
            error("Size of background must be ", string(nx), "x", string(ny))
        end
    end
    ind = findall(x -> 0.0 <= x < 1.0, pos)
    H[ind] .= pos[ind] .* (height - base) .+ base
    ind = findall(x -> x >= 1.0, pos)
    H[ind] .= height

    # dimensionalize
    Topo = dimensionalize(H, km, CharUnits)

    # build and return CartData
    return CartData(reshape(X, nx, ny, 1), reshape(Y, nx, ny, 1), reshape(Topo, nx, ny, 1), (Topography = reshape(Topo, nx, ny, 1),))
end

function make_volc_topo(
        Grid::CartData;
        center = (0, 0, 0),
        height = 0.0,
        radius = 0.0,
        crater = 0.0,
        base = 0.0,
        background = nothing
    )
    # get node grid
    X = @views Grid.x.val[:, :, 1]
    Y = @views Grid.y.val[:, :, 1]
    nx = size(X, 1)
    ny = size(X, 2)

    ## assign topography
    H = zeros(nx, ny)
    # check if there is a background supplied
    if background === nothing
        H .= base
    else
        # background = nondimensionalize(background, CharUnits)
        if size(background) == size(X)
            H .= background
        elseif size(background) == size(reshape(X, nx, ny, 1))
            H .= @views background[:, :, 1]
        else
            error("Size of background must be ", nx, "x", ny)
        end
    end

    H .= volcano_height.(X, Y, H; center, height, radius, crater, base)

    return H
end

# add_volcano! without a background topography, point by point: the volcanic phase is set where a point lies
# below the cone, then `T` at every point with a phase > 0, with depth measured downwards from the cone.
# Points with NaN coordinates (empty particle slots) are left unchanged.
function add_volcano_points!(Phases, Temp, X, Y, Z, cone, volcanic_phase, T)
    Phases .= ((ph, x, y, z) -> z < volcano_h(cone, x, y) && z ≥ cone.base ? volcanic_phase : ph).(Phases, X, Y, Z)
    isnothing(T) && return nothing

    # a lazy broadcast, so that the reduction allocates no array
    zbot = maximum(Broadcast.instantiate(Broadcast.broadcasted((ph, x, y, z) -> ph > 0 && !isnan(z) ? volcano_depth(cone, x, y, z) : -Inf, Phases, X, Y, Z)))
    zbot == -Inf && return nothing
    cT = prepare_temperature(T, volcano_frame(cone.center, cone.radius, zbot), nothing)
    Temp .= ((t, ph, x, y, z) -> ph > 0 && !isnan(z) ? point_temperature(T, cT, x, y, volcano_depth(cone, x, y, z)) : t).(Temp, Phases, X, Y, Z)
    return nothing
end

volcano_h(cone, x, y) = volcano_height(x, y, cone.base; cone...)

volcano_depth(cone, x, y, z) = max(volcano_h(cone, x, y) - z, 0)

# Frame of the thermal structure of a volcano: the footprint of the cone, with depth measured downwards from the
# topography, so that the surface is at depth 0
volcano_frame(center, radius, zbot) =
    (; xmin = center[1] - radius, xmax = center[1] + radius, ymin = center[2] - radius, ymax = center[2] + radius, ztop = 0.0, zbot)

# Topography at (x, y) of a volcanic cone with a crater, on top of the topography `h0` away from the cone
function volcano_height(x, y, h0; center, height, radius, crater, base)
    RD = √((x - center[1])^2 + (y - center[2])^2) - crater   # radial distance from the crater rim
    pos = -RD / (radius - crater) + 1                          # position relative to the crater rim
    return 0 ≤ pos < 1 ? pos * (height - base) + base : pos ≥ 1 ? height : h0
end

abstract type AbstractThermalStructure end


"""
    ConstantTemp(T=1000)

Sets a constant temperature inside the box

Parameters
===
- T : the value
"""
@with_kw_noshow struct ConstantTemp{A} <: AbstractThermalStructure
    T::A = 1000
end

# Every thermal structure is evaluated point by point. `prepare_temperature(s, frame, points)` computes,
# once per call, the reference values a structure needs from `frame`, the extent of the region
# (xmin, xmax, ymin, ymax, zbot, ztop); `points = (X, Y, Z, Phase)` holds the points of the region for
# structures that depend on them (`LithosphericTemp`), and is `nothing` otherwise.
# `point_temperature` returns the temperature at one point. `compute_thermal_structure` applies both to
# all points, taking the extent of the points as frame where `ztop` and `zbot` are not given.
compute_thermal_structure(Temp, X, Y, Z, Phase, s::AbstractThermalStructure; ztop = nothing, zbot = nothing) =
    point_temperatures!(Temp, X, Y, Z, s, prepare_temperature(s, points_frame(X, Y, Z; ztop, zbot), (X, Y, Z, Phase)))

function points_frame(X, Y, Z; ztop = nothing, zbot = nothing)
    (xmin, xmax), (ymin, ymax) = extrema(X), extrema(Y)
    return (; xmin, xmax, ymin, ymax, ztop = isnothing(ztop) ? maximum(Z) : ztop, zbot = isnothing(zbot) ? minimum(Z) : zbot)
end

# Thermal structures whose preparation needs the points of the region, not only its extent
needs_points(s) = false

function point_temperatures!(Temp, X, Y, Z, s, c)
    for i in eachindex(Temp, X, Y, Z)
        Temp[i] = point_temperature(s, c, X[i], Y[i], Z[i])
    end
    return Temp
end

prepare_temperature(s::ConstantTemp, frame, points) = nothing
point_temperature(s::ConstantTemp, c, x, y, z) = s.T


"""
    LinearTemp(Ttop=0, Tbot=1000)

Set a linear temperature structure from top to bottom.

`Ttop` is applied at the top of the region the structure is added to (e.g. the top of the
box or polygon) and `Tbot` at its bottom, wherever these are located; the temperature
varies linearly in between.

Parameters
===
- Ttop : the value @ the top
- Tbot : the value @ the bottom

"""
@with_kw_noshow struct LinearTemp{A, B} <: AbstractThermalStructure
    Ttop::A = 0
    Tbot::B = 1350
end

# Anchored at the top and bottom of the region
prepare_temperature(s::LinearTemp, frame, points) = (; frame.ztop, dz = frame.ztop - frame.zbot, dT = s.Tbot - s.Ttop)

# a region without vertical extent has no gradient to apply
point_temperature(s::LinearTemp, c, x, y, z) = c.dz == 0 ? s.Ttop : abs((z - c.ztop) / c.dz) * c.dT + s.Ttop

"""
    HalfspaceCoolingTemp(Tsurface=0, Tmantle=1350, Age=60, Adiabat=0)

Sets a halfspace temperature structure in plate.

`Tsurface` is applied at the top of the region the structure is added to (e.g. the top of the
box or polygon), and the depth used in the cooling model is measured from there.

Parameters
========
- Tsurface : surface temperature [C]
- Tmantle : mantle temperature [C]
- Age : Thermal Age of plate [Myrs]
- Adiabat : Mantle Adiabat [K/km]

"""
@with_kw_noshow struct HalfspaceCoolingTemp{A, B, C, D} <: AbstractThermalStructure
    Tsurface::A = 0       # top T
    Tmantle::B = 1350     # bottom T
    Age::C = 60           # thermal age of plate [in Myrs]
    Adiabat::D = 0        # Adiabatic gradient in K/km
end

# Depth is measured from the top of the region (the surface of the cooling halfspace)
function prepare_temperature(s::HalfspaceCoolingTemp, frame, points)
    kappa = 1.0e-6
    SecYear = 3600 * 24 * 365
    ThermalAge = s.Age * 1.0e6 * SecYear
    return (; frame.ztop, scale = 2 * sqrt(kappa * ThermalAge))
end

point_temperature(s::HalfspaceCoolingTemp, c, x, y, z) =
    (s.Tsurface - s.Tmantle) * erfc((abs(z - c.ztop) * 1.0e3) / c.scale) + (s.Tmantle + s.Adiabat * abs(z))


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

Note: the thermal age at the mid oceanic ridge is set to 1 year to avoid division by zero

"""
@with_kw_noshow struct SpreadingRateTemp{A, B, C, D, E, F} <: AbstractThermalStructure
    Tsurface::A = 0         # top T
    Tmantle::B = 1350       # bottom T
    Adiabat::C = 0          # Adiabatic gradient in K/km
    MORside::String = "left"  # side of box where the MOR is located
    SpreadingVel::D = 3     # spreading velocity [cm/yr]
    AgeRidge::E = 0         # Age of the ridge [Myrs]
    maxAge::F = 60          # maximum thermal age of plate [Myrs]
end

# The plate age grows with the distance from the ridge, which lies at the edge of the region on the side
# `MORside`: the minimum x ("left"), maximum x ("right"), minimum y ("front") or maximum y ("back").
function prepare_temperature(s::SpreadingRateTemp, frame, points)
    side = findfirst(==(s.MORside), ("left", "right", "front", "back"))
    isnothing(side) && error("unknown side")
    ref = (frame.xmin, frame.xmax, frame.ymin, frame.ymax)[side]
    return (; frame.ztop, side, ref, segments = nothing)
end

# Ridge segments: the plate age grows with the distance from the segment whose region contains the point.
function prepare_temperature(s::SpreadingRateTemp, segments, frame, points)
    delimiters = [(segments[i][2], segments[i + 1][1]) for i in 1:(length(segments) - 1)]
    return (; frame.ztop, segments, delimiters)
end

function point_temperature(s::SpreadingRateTemp, c, x, y, z)
    kappa = 1.0e-6
    SecYear = 3600 * 24 * 365
    if isnothing(c.segments)
        Distance = c.side == 1 ? x - c.ref : c.side == 2 ? c.ref - x : c.side == 3 ? y - c.ref : c.ref - y
        ThermalAge = abs(Distance * 1.0e3 * 1.0e2) / s.SpreadingVel + s.AgeRidge * 1.0e6    # Thermal age in years
    else
        region = determine_region(x, y, c.delimiters, c.segments)
        x1, y1 = c.segments[region][1]
        x2, y2 = c.segments[region][2]
        Distance = perpendicular_distance_to_segment(x, y, x1, y1, x2, y2)
        ThermalAge = abs(Distance * 1.0e5) / s.SpreadingVel + s.AgeRidge * 1.0e6  # Thermal age in years
    end
    if ThermalAge > s.maxAge * 1.0e6
        ThermalAge = s.maxAge * 1.0e6
    end
    ThermalAge = ThermalAge * SecYear
    if ThermalAge == 0
        ThermalAge = 1.0e-6   # doesn't like zero
    end
    return (s.Tsurface - s.Tmantle) * erfc((abs(z - c.ztop) * 1.0e3) / (2 * sqrt(kappa * ThermalAge))) + (s.Tmantle + s.Adiabat * abs(z))
end


"""
    compute_thermal_structure(Temp, X, Y, Z, Phase, s::SpreadingRateTemp, segments)

Calculates the temperature distribution across the plate considering multiple ridge segments.

This function computes the thermal structure based on the perpendicular distance from each point to its corresponding ridge segment, and applies a thermal model using a spreading velocity and thermal age.

Parameters
==========
- Temp     : Temperature field to be updated (array)
- X, Y, Z  : Coordinates of the points (arrays)
- Phase    : Phase of the material (unused)
- s        : SpreadingRateTemp object containing the thermal and spreading parameters
- segments : List of ridge segments, where each segment is defined by two tuples representing the start and end coordinates (x1, y1) and (x2, y2) for each segment.

Note
====
The temperature at each point is calculated using the thermal age, which is determined by the distance from the point to the nearest ridge segment and the spreading velocity.

The function works with one or more segments, and the ridge can be placed at any position within the box, not only at its boundary.

The thermal age is capped at `maxAge` years, and the temperature is adjusted based on the distance to the ridge and the corresponding thermal gradient.

"""
compute_thermal_structure(Temp, X, Y, Z, Phase, s::SpreadingRateTemp, segments::AbstractVector; ztop = nothing, zbot = nothing) =
    point_temperatures!(Temp, X, Y, Z, s, prepare_temperature(s, segments, points_frame(X, Y, Z; ztop, zbot), (X, Y, Z, Phase)))


# Supporting functions for multi-segment ridge functionality

# Function to calculate the perpendicular distance from a point to a segment
function perpendicular_distance_to_segment(x, y, x1, y1, x2, y2)
    num = abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    den = sqrt((y2 - y1)^2 + (x2 - x1)^2)
    return num / den
end

# Function to determine the side of a point with respect to a line (adjusted for segment direction)
function side_of_line(x, y, x1, y1, x2, y2, direction)
    side = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    return direction == :left ? side > 0 : side < 0
end

# Function to determine in which region a point lies (based on delimiters)
function determine_region(px, py, delimiters, segments)
    for i in 1:length(delimiters)
        x1, y1 = delimiters[i][1]
        x2, y2 = delimiters[i][2]

        # Determine the direction of the segments
        direction = x2 < x1 ? :left : :right

        # Check the side of the line considering the direction
        if side_of_line(px, py, x1, y1, x2, y2, direction)
            return i  # Region corresponding to segment i
        end
    end
    return length(segments)  # Last region
end

"""
    LithosphericTemp(Tsurface=0.0, Tpot=1350.0, dTadi=0.5,
                        ubound="const", lbound="const, utbf = 50.0e-3, ltbf = 10.0e-3,
                        age = 120.0, dtfac = 0.9, nz = 201,
                        rheology = example_CLrheology()
                    )

Calculates a 1D temperature profile [C] for variable thermal parameters including radiogenic heat source and
    linearly interpolates the temperature profile onto the box. The thermal parameters are defined in
    rheology and the structure of the lithosphere is define by LithosphericPhases().


Parameters
========
- Tsurface  : surface temperature [C]
- Tpot      : potential mantle temperature [C]
- dTadi     : adiabatic gradient [K/km]
- ubound    : Upper thermal boundary condition ["const","flux"]
- lbound    : Lower thermal boundary condition ["const","flux"]
- utbf      : Upper thermal heat flux [W/m]; if ubound == "flux"
- ltbf      : Lower thermal heat flux [W/m]; if lbound == "flux"
- age       : age of the lithosphere [Ma]
- dtfac     : Diffusion stability criterion to calculate T_age
- nz        : Grid spacing for the 1D profile within the box
- rheology  : Structure containing the thermal parameters for each phase [default example_CLrheology]

"""
@with_kw_noshow struct LithosphericTemp <: AbstractThermalStructure
    Tsurface = 0.0      # top T [C]
    Tpot = 1350.0       # potential T [C]
    dTadi = 0.5         # adiabatic gradient in K/km
    ubound = "const"    # Upper thermal boundary condition
    lbound = "const"    # lower thermal boundary condition
    utbf = 50.0e-3      # q [W/m^2]; if ubound = "flux"
    ltbf = 10.0e-3      # q [W/m^2]; if lbound = "flux"
    age = 120.0         # Lithospheric age [Ma]
    dtfac = 0.9         # Diffusion stability criterion
    nz = 201
    rheology = example_CLrheology()
end

struct Thermal_parameters{A}
    ρ::A
    Cp::A
    k::A
    ρCp::A
    H::A
    function Thermal_parameters(ni)
        ρ = zeros(ni)
        Cp = zeros(ni)
        k = zeros(ni)
        ρCp = zeros(ni)
        H = zeros(ni)
        return new{typeof(ρ)}(ρ, Cp, k, ρCp, H)
    end
end

needs_points(s::LithosphericTemp) = true

function prepare_temperature(s::LithosphericTemp, frame, points)
    X, Y, Z, Phase = points
    @unpack Tsurface, Tpot, dTadi, ubound, lbound, utbf, ltbf, age,
        dtfac, nz, rheology = s

    # Create 1D depth profile within the box
    z = LinRange(maximum(Z) * 1.0e3, minimum(Z) * 1.0e3, nz)    # [m]
    dz = z[2] - z[1]                                         # Gride resolution
    # Initialize 1D arrays for explicit solver
    T = zeros(nz)
    phase = Int64.(zeros(nz))

    # Assign phase id from Phase to 1D phase array
    phaseid = (minimum(Phase):1:maximum(Phase))
    zsurf = maximum(Z[findall(Phase .== phaseid[1])]) * 1.0e3
    zbase = zeros(length(phaseid)) # base of each layer

    # for each phase id
    for i in 1:length(phaseid)
        # Calculate layer thickness from Phase array
        zbase[i] = minimum(Z[findall(Phase .== phaseid[i])]) * 1.0e3
    end
    for i in 1:length(phaseid)
        # Assign phase ids
        ztop = i === 1 ? zsurf : zbase[i - 1]
        ind = findall((z .>= zbase[i]) .& (z .<= ztop))
        phase[ind] .= phaseid[i]
    end

    # Setup initial T-profile
    Tpot = Tpot + 273.15                   # Potential temp [K]
    Tsurface = Tsurface + 273.15               # Surface temperature [ K ]
    T = @. Tpot + abs.(z ./ 1.0e3) * dTadi  # Initial T-profile [ K ]
    T[1] = Tsurface

    args = (;)
    thermal_parameters = Thermal_parameters(nz)

    ## Update thermal parameters ======================================== #
    compute_density!(thermal_parameters.ρ, rheology, phase, args)
    compute_heatcapacity!(thermal_parameters.Cp, rheology, phase, args)
    compute_conductivity!(thermal_parameters.k, rheology, phase, args)
    thermal_parameters.ρCp .= @. thermal_parameters.Cp * thermal_parameters.ρ
    compute_radioactive_heat!(thermal_parameters.H, rheology, phase, args)

    # Thermal diffusivity [ m^2/s ]
    κ = maximum(thermal_parameters.k) /
        minimum(thermal_parameters.ρ) / minimum(thermal_parameters.Cp)
    ## =================================================================== #
    ## Time stability criterion ========================================= #
    tfac = 60.0 * 60.0 * 24.0 * 365.25   # Seconds per year
    age = age * 1.0e6 * tfac          # Age in seconds
    dtexp = dz^2.0 / 2.0 / κ            # Stability criterion for explicit
    dt = dtfac * dtexp             # [s]
    nit = Int64(ceil(age / dt))     # Number of iterations
    time = zeros(nit)              # Time array

    for i in 1:nit
        if i > 1
            time[i] = time[i - 1] + dt
        end
        SolveDiff1Dexplicit_vary!(
            T,
            thermal_parameters,
            ubound, lbound,
            utbf, ltbf,
            dz,
            dt
        )
    end

    # Flat extrapolation: the extreme points can lie one ulp outside the nodes after the m/km round trip
    interp_linear_T = linear_interpolation(-z ./ 1.0e3, T .- 273.15; extrapolation_bc = Flat())

    return (; itp = interp_linear_T)
end

point_temperature(s::LithosphericTemp, c, x, y, z) = c.itp(-z)

function SolveDiff1Dexplicit_vary!(
        T,
        thermal_parameters,
        ubound, lbound,
        utbf, ltbf,
        di,
        dt
    )
    nz = length(T)
    T0 = T

    if ubound == "const"
        T[1] = T0[1]
    elseif ubound == "flux"
        kB = (thermal_parameters.k[2] + thermal_parameters.k[1]) / 2.0
        kA = (thermal_parameters.k[1] + thermal_parameters.k[1]) / 2.0
        a = (dt * (kA + kB)) / (di^2.0 * thermal_parameters.ρCp[1])
        b = 1 - (dt * (kA + kB)) / (di^2.0 * thermal_parameters.ρCp[1])
        c = (dt * 2.0 * utbf) / (di * thermal_parameters.ρCp[1])
        T[1] = a * T0[2] + b * T0[1] + c +
            thermal_parameters.H[1] * dt / thermal_parameters.ρCp[1]
    end
    if lbound == "const"
        T[nz] = T0[nz]
    elseif lbound == "flux"
        kB = (thermal_parameters.k[nz] + thermal_parameters.k[nz]) / 2.0
        kA = (thermal_parameters.k[nz] + thermal_parameters.k[nz - 1]) / 2.0
        a = (dt * (kA + kB)) / (di^2.0 * thermal_parameters.ρCp[nz])
        b = 1 - (dt * (kA + kB)) / (di^2.0 * thermal_parameters.ρCp[nz])
        c = -(dt * 2.0 * ltbf) / (di * thermal_parameters.ρCp[nz])
        T[nz] = a * T0[nz - 1] + b * T0[nz] + c
    end

    kAi = @. (thermal_parameters.k[1:(end - 2)] + thermal_parameters.k[2:(end - 1)]) / 2.0
    kBi = @. (thermal_parameters.k[2:(end - 1)] + thermal_parameters.k[3:end]) / 2.0
    ai = @. (kBi * dt) / (di^2.0 * thermal_parameters.ρCp[2:(end - 1)])
    bi = @. 1.0 - (dt * (kAi + kBi)) / (di^2.0 * thermal_parameters.ρCp[2:(end - 1)])
    ci = @. (kAi * dt) / (di^2.0 * thermal_parameters.ρCp[2:(end - 1)])
    T[2:(end - 1)] = @. ai * T0[3:end] + bi * T0[2:(end - 1)] + ci * T0[1:(end - 2)] +
        thermal_parameters.H[2:(end - 1)] * dt / thermal_parameters.ρCp[2:(end - 1)]
    return T
end

function example_CLrheology(;
        ρM = 3.0e3,           # Density [ kg/m^3 ]
        CpM = 1.0e3,          # Specific heat capacity [ J/kg/K ]
        kM = 2.3,             # Thermal conductivity [ W/m/K ]
        HM = 0.0,             # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
        ρUC = 2.7e3,          # Density [ kg/m^3 ]
        CpUC = 1.0e3,         # Specific heat capacity [ J/kg/K ]
        kUC = 3.0,            # Thermal conductivity [ W/m/K ]
        HUC = 617.0e-12,      # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
        ρLC = 2.9e3,          # Density [ kg/m^3 ]
        CpLC = 1.0e3,         # Specific heat capacity [ J/kg/K ]
        kLC = 2.0,            # Thermal conductivity [ W/m/K ]
        HLC = 43.0e-12,       # Radiogenic heat source per mass [H] = W/kg; [H] = [Q/rho]
    )

    rheology = (
        # Name              = "UpperCrust",
        SetMaterialParams(;
            Phase = 1,
            Density = ConstantDensity(; ρ = ρUC),
            HeatCapacity = ConstantHeatCapacity(; Cp = CpUC),
            Conductivity = ConstantConductivity(; k = kUC),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = HUC * ρUC),     # [H] = W/m^3
        ),
        # Name              = "LowerCrust",
        SetMaterialParams(;
            Phase = 2,
            Density = ConstantDensity(; ρ = ρLC),
            HeatCapacity = ConstantHeatCapacity(; Cp = CpLC),
            Conductivity = ConstantConductivity(; k = kLC),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = HLC * ρLC),     # [H] = W/m^3
        ),
        # Name              = "LithosphericMantle",
        SetMaterialParams(;
            Phase = 3,
            Density = ConstantDensity(; ρ = ρM),
            HeatCapacity = ConstantHeatCapacity(; Cp = CpM),
            Conductivity = ConstantConductivity(; k = kM),
            RadioactiveHeat = ConstantRadioactiveHeat(; H_r = HM * ρM),       # [H] = W/m^3
        ),
    )
    return rheology
end

abstract type AbstractPhaseNumber end

# Thermal and phase structures are shown as a call of their constructor, without type parameters
function show(io::IO, s::Union{AbstractThermalStructure, AbstractPhaseNumber})
    print(io, nameof(typeof(s)), "(")
    for (i, f) in enumerate(fieldnames(typeof(s)))
        i > 1 && print(io, ", ")
        show(io, getfield(s, f))
    end
    return print(io, ")")
end


"""
    ConstantPhase(phase=1)

Sets a constant phase inside the box

Parameters
===
- phase : the value
"""
@with_kw_noshow struct ConstantPhase{P} <: AbstractPhaseNumber
    phase::P = 1
end

# Every phase structure is evaluated point by point by `point_phase`, from the previous phase, the
# temperature and the coordinates of the point; `compute_phase` applies it to all points.
function compute_phase(Phase, Temp, X, Y, Z, s::AbstractPhaseNumber; kwargs...)
    for i in eachindex(Phase, Temp, X, Y, Z)
        Phase[i] = point_phase(s, Phase[i], Temp[i], X[i], Y[i], Z[i]; kwargs...)
    end
    return Phase
end

point_phase(s::ConstantPhase, ph, T, x, y, z) = s.phase


"""
    LithosphericPhases(Layers=[10 20 15], Phases=[1 2 3 4], Tlab=nothing )

This allows defining a layered lithosphere. Layering is defined from the top downwards.

Parameters
===
- Layers : The thickness of each layer, ordered from top to bottom. The thickness of the last layer does not have to be specified.
- Phases : The phases of the layers, ordered from top to bottom.
- Tlab   : Temperature of the lithosphere asthenosphere boundary. If specified, the phases at locations with T>Tlab are set to Phases[end].

"""
@with_kw_noshow struct LithosphericPhases{L, P, T} <: AbstractPhaseNumber
    Layers::L = [10.0, 20.0, 15.0]
    Phases::P = [1, 2, 3, 4]
    Tlab::T = nothing
end


"""
    Phase = compute_phase(Phase, Temp, X, Y, Z, s::LithosphericPhases; Ztop = 0)

or

    Phase = compute_phase(Phase, Temp, Grid::AbstractGeneralGrid, s::LithosphericPhases)

This copies the layered lithosphere onto the Phase matrix.

Parameters
===
- Phase - Phase array
- Temp  - Temperature array
- X     - x-coordinate array (consistent with Phase and Temp)
- Y     - y-coordinate array (consistent with Phase and Temp)
- Z     - Vertical coordinate array (consistent with Phase and Temp)
- s     - LithosphericPhases
- Ztop  - Vertical coordinate of top of model box
- Grid  - Grid structure (usually obtained with read_LaMEM_inputfile)
"""
compute_phase(Phase, Temp, X, Y, Z, s::LithosphericPhases)

# Layers are counted downward from Ztop; where two layers share a boundary, the lower one applies.
function point_phase(s::LithosphericPhases, ph, T, x, y, z; Ztop = 0)
    p = s.Phases[end]
    for i in eachindex(s.Layers)
        Zbot = Ztop - s.Layers[i]
        if Zbot <= z <= Ztop
            p = s.Phases[i]
        end
        Ztop = Zbot
    end

    # set phase to mantle if requested
    if !isnothing(s.Tlab) && T > s.Tlab
        p = s.Phases[end]
    end
    return p
end

# allow AbstractGeneralGrid instead of Z and Ztop
compute_phase(Phase, Temp, Grid::LaMEM_grid, s::LithosphericPhases) = compute_phase(Phase, Temp, Grid.X, Grid.Y, Grid.Z, s::LithosphericPhases, Ztop = maximum(Grid.coord_z))


"""
    McKenzie_subducting_slab

Thermal structure by McKenzie for a subducted slab that is fully embedded in the mantle.

Parameters
===
- `Tsurface`:     Top T [C]
- `Tmantle`:      Bottom T [C]
- `Adiabat`:      Adiabatic gradient in K/km
- `v_cm_yr`:      Subduction velocity [cm/yr]
- `κ`:            Thermal diffusivity [m2/s]
- `it`:           Number iterations employed in the harmonic summation

"""
@with_kw_noshow struct McKenzie_subducting_slab <: AbstractThermalStructure
    Tsurface::Float64 = 20.0       # top T
    Tmantle::Float64 = 1350.0     # bottom T
    Adiabat::Float64 = 0.4        # Adiabatic gradient in K/km
    v_cm_yr::Float64 = 2.0        # velocity of subduction [cm/yr]
    κ::Float64 = 1.0e-6       # Thermal diffusivity [m2/s]
    it::Int64 = 36         # number of harmonic summation (look Mckenzie formula)
end

"""
    compute_thermal_structure(Temp, X, Y, Z, Phase, s::McKenzie_subducting_slab)

Compute the temperature field of a `McKenzie_subducting_slab`. Uses the analytical solution
of McKenzie (1969) ["Speculations on the consequences and causes of plate motions"]. Depth is measured from
the top of the region, and the slab thickness is the height of the region.

Parameters

=============================
- `Temp`:  Temperature array
- `X`:    X Array
- `Y`:    Y Array
- `Z`:    Z Array
- `Phase`: Phase array
- `s`:    `McKenzie_subducting_slab`
"""
compute_thermal_structure(Temp, X, Y, Z, Phase, s::McKenzie_subducting_slab)

# The McKenzie model is defined with Z = 0 at the top of the region, and its thickness is that of the region
function prepare_temperature(s::McKenzie_subducting_slab, frame, points)
    Thickness = (frame.ztop - frame.zbot)

    # Convert subduction velocity from cm/yr -> m/s;
    convert_velocity = 1 / (100.0 * 365.25 * 60.0 * 60.0 * 24.0)
    v_s = s.v_cm_yr * convert_velocity

    # calculate the thermal Reynolds number
    Re = (v_s * Thickness * 1000) / 2 / s.κ      # factor 1000 to transfer Thickness from km to m

    # coefficients 2a, r and iπ of the terms of the harmonic sum
    terms = ntuple(i -> (2 * ((-1.0)^i / (i * pi)), Re - (Re^2 + i^2.0 * pi^2.0)^0.5, i * pi), s.it)
    return (; Zref = frame.ztop, sc = 1 / Thickness, terms)
end

function point_temperature(s::McKenzie_subducting_slab, c, x, y, z)
    σ = 1.0
    zshift = z - c.Zref
    for (a2, r, iπ) in c.terms
        σ += a2 * exp(r * x * c.sc) * sin(iπ * (1 - abs(zshift * c.sc)))
    end
    return s.Tsurface + (s.Tmantle - s.Tsurface) * σ + s.Adiabat * abs(z)
end

"""
    LinearWeightedTemperature

Structure that defined a linear average temperature between two temperature fields as a function of distance

Parameters
===
- w_min:        Minimum weight
- w_max:        Maximum weight
- crit_dist:    Critical distance
- dir:          Direction of the averaging (`:X`, `:Y` or `:Z`)
- F1:           First temperature field
- F2:           Second temperature field

"""
@with_kw_noshow struct LinearWeightedTemperature{T1 <: AbstractThermalStructure, T2 <: AbstractThermalStructure} <: AbstractThermalStructure
    w_min::Float64 = 0.0
    w_max::Float64 = 1.0
    crit_dist::Float64 = 100.0
    dir::Symbol = :X
    F1::T1 = ConstantTemp()
    F2::T2 = ConstantTemp()
end

"""
    compute_thermal_structure(Temp, X, Y, Z, Phase, s::LinearWeightedTemperature)

Weight average along distance

Do a weight average between two field along a specified direction

Given a distance (could be any array, from X,Y) -> the weight of F1 increase from the origin, while F2 decreases.

This function has been conceived for averaging the solution of McKenzie and half space cooling models, but it
can be used to smooth the temperature field from continent ocean:
- Select the boundary to apply;
- transform the coordinate such that dist represent the perpendicular direction along which you want to apply this smoothening
  and in a such way that 0.0 is the point in which the weight of F1 is equal to 0.0;
- Select the points that belongs to this area
- compute the thermal fields {F1} {F2}
- then modify F.
"""
compute_thermal_structure(Temp, X, Y, Z, Phase, s::LinearWeightedTemperature)

prepare_temperature(s::LinearWeightedTemperature, frame, points) =
    (; c1 = prepare_temperature(s.F1, frame, points), c2 = prepare_temperature(s.F2, frame, points))
needs_points(s::LinearWeightedTemperature) = needs_points(s.F1) || needs_points(s.F2)

function point_temperature(s::LinearWeightedTemperature, c, x, y, z)
    dist = s.dir === :X ? x : s.dir === :Y ? y : z
    weight = s.w_min + (s.w_max - s.w_min) / (s.crit_dist) * (dist)
    if weight > s.w_max
        weight = s.w_max
    end
    if weight < s.w_min
        weight = s.w_min
    end
    return point_temperature(s.F1, c.c1, x, y, z) * (1.0 - weight) + point_temperature(s.F2, c.c2, x, y, z) * weight
end


abstract type AbstractTrenchSlab end

"""
    Trench structure

Structure that defines the geometry of the trench and the slab.

Parameters
===

- `Start`     - Start of the trench (`x`,`y`) coordinates
- `End`       - End of the trench (`x`,`y`) coordinates
- `n_seg`     - The number of segment through which the slab is discretize along the dip
- `Length`    - The length of the slab
- `Thickness` - The thickness of the slab
- `Lb`        - Critical distance through which apply the bending angle functions Lb ∈ [0,Length];
- `θ_max`     - maximum angle of bending ∈ [0°,90°].
- `direction` - the direction of the dip
               The rotation of the coordinate system is done as such that the new X is parallel to the segment. Since the
               rotation is anticlockwise the coordinate y has specific values: direction tells if the subduction is directed along
               the positive or negative direction of the new y coordinate system. In practice, it apply an additional transformation
               to y by multiplying it with -1 or +1;
- `d_decoupling` - depth at which the slab is fully submerged into the mantle.
- `type_bending` - is the type of bending angle of the slab [`:Linear`, `:Ribe`].
    The angle of slab changes as a function of `l` (∈ [0,Length]). `l` is the actual distance along the slab length from
    the trench.
    In case:
        - `:Linear`
            ```math θ(l) = ((θ_max - 0.0)/(Lb-0))*l ```;
        - `:Ribe`
            ```math θ(l) =  θ_max*l^2*((3*Lb-2*l))/(Lb^3) ```;
            which is taken from Ribe 2010 [Bending mechanics and mode selection in free subduction: a thin-sheet analysis]

    For l>Lb, θ(l) = θ_max;
- `WeakzoneThickness` - Thickness of the weakzone [km]
- `WeakzonePhase` - Phase of the weakzone

"""
@with_kw_noshow struct Trench{Nseg} <: AbstractTrenchSlab
    Start::NTuple{Nseg, Float64} = (0.0, 0.0)   # Start (x,y) coordinates of trench (in mapview)
    End::NTuple{Nseg, Float64} = (0.0, 1.0)     # End (x,y) coordinates of trench (in mapview)
    n_seg::Int64 = 50                          # number of segments in downdip direction
    Length::Float64 = 400.0                   # length of the slab
    Thickness::Float64 = 100.0                # thickness of the slab
    Lb::Float64 = 200.0                       # Length at which all the bending is happening (Lb<=Length)
    θ_max::Float64 = 45.0                      # max bending angle, (must be converted into radians)
    direction::Float64 = -1.0                  # Direction of the bending angle (-1= left to right or 1.0=right to left)
    d_decoupling::Float64 = 100               # decoupling depth of the slab
    type_bending::Symbol = :Ribe               # Mode Ribe | Linear | Customize
    WeakzoneThickness::Float64 = 0.0           # Thickness of the weakzone
    WeakzonePhase::Int64 = 5                   # Phase of the weak zone
end

function show(io::IO, g::Trench{Nseg}) where {Nseg}
    println(io, "Trench{$Nseg}, $(g.n_seg) segments")
    println(io, "     Trench [Start/End] : $(g.Start) - $(g.End) [km]")
    println(io, "       Slab [Thickness] : $(g.Thickness) km")
    println(io, "          Slab [Length] : $(g.Length) km")
    println(io, "    Bending length [Lb] : $(g.Lb) km")
    println(io, "     Max. angle [θ_max] : $(g.θ_max)ᵒ")
    if g.direction == -1.0
        println(io, "        Dip [direction] : left to right [$(g.direction)]")
    else
        println(io, "     Dip [direction] : right to left [$(g.direction)]")
    end
    println(io, "    Depth [d_decoupling]: $(g.d_decoupling) km")
    println(io, "  Bending [type_bending]: $(g.type_bending)")
    println(io, "    [WeakzoneThickness] : $(g.WeakzoneThickness) km")
    if g.WeakzoneThickness > 0
        println(io, "   Weakzone phase : $(g.WeakzonePhase)")
    end

    return nothing
end

"""
    Top, Bot = compute_slab_surface(trench::Trench)

Computes the (`x`,`z`) coordinates of the slab top, bottom surface using the mid surface of the slab as reference.

Parameters
===
- `trench`          - `Trench` structure that contains the relevant parameters

Method
===

It computes it by discretizing the slab surface in `n_seg` segments, and computing the average bending angle (which is a function of the current length of the slab).
Next, it compute the coordinates assuming that the trench is at 0.0, and assuming a positive `θ_max` angle.
"""
function compute_slab_surface(trench::Trench)

    @unpack Thickness, Length, n_seg, Lb, θ_max, type_bending, direction, WeakzoneThickness = trench

    # Convert θ_max into radians
    θ_max *= π / 180

    # Allocate the top, mid and bottom surface
    Top = zeros(n_seg + 1, 2)
    Bottom = zeros(n_seg + 1, 2)
    WeakZone = zeros(n_seg + 1, 2)
    Bottom[1, 2] = -Thickness
    WeakZone[1, 2] = WeakzoneThickness
    MidS = zeros(n_seg + 1, 2)
    MidS[1, 2] = -Thickness / 2

    # Initialize the length.
    dl = Length / n_seg  # dl
    for it in 1:n_seg
        l = (it - 1) * dl

        # Compute the mean angle within the segment
        θ = compute_bending_angle(θ_max, Lb, l, type_bending)
        θ_n = compute_bending_angle(θ_max, Lb, l + dl, type_bending)
        θ_mean = (θ + θ_n) / 2

        # Mid surface coordinates (x,z)
        sinθ, cosθ = sincos(θ_mean)

        MidS[it + 1, 1] = MidS[it, 1] + dl * cosθ
        MidS[it + 1, 2] = MidS[it, 2] - dl * sinθ

        # Top surface coordinates (x,z)
        Top[it + 1, 1] = MidS[it + 1, 1] + 0.5 * Thickness * abs(sinθ)
        Top[it + 1, 2] = MidS[it + 1, 2] + 0.5 * Thickness * abs(cosθ)

        # Bottom surface coordinate
        Bottom[it + 1, 1] = MidS[it + 1, 1] - 0.5 * Thickness * abs(sinθ)
        Bottom[it + 1, 2] = MidS[it + 1, 2] - 0.5 * Thickness * abs(cosθ)

        # Compute the top surface for the weak zone
        WeakZone[it + 1, 1] = MidS[it + 1, 1] + (0.5 * Thickness + WeakzoneThickness) * abs(sinθ)
        WeakZone[it + 1, 2] = MidS[it + 1, 2] + (0.5 * Thickness + WeakzoneThickness) * abs(cosθ)

    end
    Top[:, 1] *= direction
    Bottom[:, 1] *= direction
    WeakZone[:, 1] *= direction

    return Top, Bottom, WeakZone
end

"""
    θ = compute_bending_angle(θ_max,Lb,l,type)

function that computes the bending angle `θ` as a function of length along the slab `l`.

Parameters
===
`θ_max` = maximum bending angle
`Lb`    = length at which the function of bending is applied (Lb<=Length)
`l`     = current position within the slab
`type`  = type of bending [`:Ribe`,`:Linear`]

"""
function compute_bending_angle(θ_max::Float64, Lb::Float64, l::Float64, type::Symbol)

    if l > Lb
        θ = θ_max
    elseif type === :Ribe
        # Compute theta
        θ = θ_max * l^2 * ((3 * Lb - 2 * l)) / (Lb^3)
    elseif type === :Linear
        # Compute the actual angle
        θ = l * (θ_max - 0) / (Lb)
    end
    return θ
end

# Bounding box, in the plane perpendicular to the trench, of the slab between the polylines `Top` and `Bottom`
function slab_bounding_box(Top, Bottom)
    (ymin, ymax), (zmin, zmax) = extrema([Top[:, 1]; Bottom[:, 1]]), extrema([Top[:, 2]; Bottom[:, 2]])
    return (; ymin, ymax, zmin, zmax)
end

# Corner points, polygon, bounding box and starting length `l` along the slab of every segment of the slab
# between the polylines `Top` and `Bottom` (n_seg + 1 points each in the plane perpendicular to the trench)
function slab_segments(Top, Bottom, n_seg, dl)
    l = [zero(dl)]  # length at the trench position, then summed segment by segment
    for i in 2:n_seg
        push!(l, l[end] + dl)
    end
    segments = map(i -> slab_segment(Top, Bottom, i, l[i]), 1:n_seg)
    # the bounding boxes separately, as most points are tested only against them
    return (; boxes = map(sg -> (sg.ymin, sg.ymax, sg.zmin, sg.zmax), segments), segments)
end

function slab_segment(Top, Bottom, i, l)
    pa = (Top[i, 1], Top[i, 2])        # D = 0 | L = l
    pb = (Bottom[i, 1], Bottom[i, 2])  # D = -Thickness | L=l
    pc = (Bottom[i + 1, 1], Bottom[i + 1, 2])  # D = -Thickness |L=L+dl
    pd = (Top[i + 1, 1], Top[i + 1, 2]) # D = 0| L = L+dl
    poly_y = @SVector [pa[1], pb[1], pc[1], pd[1]]
    poly_z = @SVector [pa[2], pb[2], pc[2], pd[2]]
    (ymin, ymax), (zmin, zmax) = extrema(poly_y), extrema(poly_z)
    return (; pa, pb, pd, poly_y, poly_z, ymin, ymax, zmin, zmax, l)
end

# Length `ls` along the slab (from the trench) and distance `d` perpendicular to its surface (0 at the
# surface, negative below it) of a point with coordinates (xr, yr, zr) in the frame of the trench, for the
# slab `segments`; NaN for a point in no slab segment. Where segments overlap, the last one applies.
function slab_distances(xr, yr, zr, xb, segments, bbox)
    (0.0 <= xr <= xb && bbox.ymin <= yr <= bbox.ymax && bbox.zmin <= zr <= bbox.zmax) || return NaN, NaN
    for i in reverse(eachindex(segments.boxes))
        ymin, ymax, zmin, zmax = segments.boxes[i]
        (ymin <= yr <= ymax && zmin <= zr <= zmax) || continue
        sg = segments.segments[i]
        if inpoly_either(sg.poly_y, sg.poly_z, yr, zr)
            return distance_to_linesegment((yr, zr), sg.pb, sg.pa) + sg.l, -distance_to_linesegment((yr, zr), sg.pa, sg.pd)
        end
    end
    return NaN, NaN
end


"""
    distance_to_linesegment(p::NTuple{2,_T}, v::NTuple{2,_T}, w::NTuple{2,_T})

Computes the distance normal distance from a point `p` to a line segment defined by the points `v` and `w`.
"""
function distance_to_linesegment(p::NTuple{2, _T}, v::NTuple{2, _T}, w::NTuple{2, _T}) where {_T <: Number}
    dx = w[1] - v[1]
    dy = w[2] - v[2]
    l2 = dx * dx + dy * dy  # i.e. |w-v|^2 -  avoid a sqrt
    if l2 == 0.0
        dx = p[1] - v[1]
        dy = p[2] - v[2]
        return sqrt(dx * dx + dy * dy)   # v == w case
    end
    # Consider the line extending the segment, parameterized as v + t (w - v).
    # We find projection of point p onto the line.
    # It falls where t = [(p-v) . (w-v)] / |w-v|^2
    t = ((p[1] - v[1]) * dx + (p[2] - v[2]) * dy) / l2
    if t < 0.0
        dx = p[1] - v[1]
        dy = p[2] - v[2]
        return sqrt(dx * dx + dy * dy)       # Beyond the 'v' end of the segment
    elseif t > 1.0
        dx = p[1] - w[1]
        dy = p[2] - w[2]
        return sqrt(dx * dx + dy * dy)  # Beyond the 'w' end of the segment
    end
    projection_x = v[1] + t * dx
    projection_y = v[2] + t * dy
    dx = p[1] - projection_x
    dy = p[2] - projection_y
    return sqrt(dx * dx + dy * dy)
end

"""
    add_slab!(Phase, Temp, Grid::AbstractGeneralGrid,  trench::Trench; phase = ConstantPhase(1), T = nothing, cell=false)

Adds a curved slab with phase & temperature structure to a 3D model setup.

Parameters
====
- `Phase`   - Phase array (consistent with Grid)
- `Temp`    - Temperature array (consistent with Grid)
- `Grid`    - grid structure (can be any of the grid types in `GMG`)
- `trench`  - Trench structure
- `phase`   - specifies the phase of the box. See `ConstantPhase()`,`LithosphericPhases()`
- `T`       - specifies the temperature of the box. See `ConstantTemp()`,`LinearTemp()`,`HalfspaceCoolingTemp()`,`SpreadingRateTemp()`,`LithosphericTemp()`
- `cell`    - if true, `Phase` and `Temp` are defined on cells

Examples
========

Example 1) Slab
```julia-repl
julia> x     = LinRange(0.0,1200.0,128);
julia> y     = LinRange(0.0,1200.0,128);
julia> z     = LinRange(-660,50,128);
julia> Cart  = CartData(xyz_grid(x, y, z));
julia> Phase = ones(Int64,size(Cart));
julia> Temp  = fill(1350.0,size(Cart));
# Define the trench:
julia> trench= Trench(Start = (400.0,400.0), End = (800.0,800.0), θ_max = 45.0, direction = 1.0, n_seg = 50, Length = 600.0, Thickness = 80.0, Lb = 500.0, d_decoupling = 100.0, type_bending =:Ribe)
julia> phase = LithosphericPhases(Layers=[5 7 88], Phases = [2 3 4], Tlab=nothing)
julia> TsHC  = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=30, Adiabat=0.4)
julia> add_slab!(Phase, Temp, Cart, trench, phase = phase, T = TsHC)
```

"""
function add_slab!(
        Phase, Temp, Grid::AbstractGeneralGrid, trench::Trench;     # required input
        phase::AbstractPhaseNumber = ConstantPhase(1),                          # Sets the phase number(s) in the slab
        T::Union{AbstractThermalStructure, Nothing} = nothing, cell = false
    )     # Sets the thermal structure (various functions are available),

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Compute top and bottom of the slab
    Top, Bottom, WeakZone = compute_slab_surface(trench)
    @unpack Start, End, n_seg = trench

    # Coordinates in the frame of the trench, rotated along the angle from Start -> End
    StrikeAngle = -atand(End[2] - Start[2], End[1] - Start[1])
    sinStrike, cosStrike = sincosd(StrikeAngle)
    xb = Rot3D(End[1] - Start[1], End[2] - Start[2], 0.0, cosd(StrikeAngle), sind(StrikeAngle), 1.0, 0.0)[1]
    dl = trench.Length / n_seg
    bbox, segments = slab_bounding_box(Top, Bottom), slab_segments(Top, Bottom, n_seg, dl)

    if isa(T, LinearWeightedTemperature)
        l_decouplingind = findall(Top[:, 2] .<= -trench.d_decoupling)
        if !isempty(l_decouplingind)
            l_decoupling = Top[l_decouplingind[1], 1]
            T = LinearWeightedTemperature(T.w_min, T.w_max, abs(l_decoupling), T.dir, T.F1, T.F2)
        end
    end

    # The slab is evaluated in the coordinates (ls, y, d): ls is the length along the slab and d the distance
    # perpendicular to it, with the slab surface at d = 0 and its base at d = -Thickness
    p = (; Start, cosStrike, sinStrike, xb, segments, bbox, Thickness = trench.Thickness)
    frame = (; xmin = 0.0, xmax = trench.Length, ymin = min(Start[2], End[2]), ymax = max(Start[2], End[2]), zbot = -trench.Thickness, ztop = 0.0)
    set_shape!(Phase, Temp, X, Y, Z, slab_coordinates, p, frame, phase, T)

    # Add a weak zone on top of the slab (indicated by a phase number but not by temperature)
    if trench.WeakzoneThickness > 0.0
        bbox_weakzone, segments_weakzone = slab_bounding_box(WeakZone, Top), slab_segments(WeakZone, Top, n_seg, dl)
        pw = (;
            Start, cosStrike, sinStrike, xb, segments = segments_weakzone, bbox = bbox_weakzone,
            WeakzoneThickness = trench.WeakzoneThickness, WeakzonePhase = trench.WeakzonePhase, d_decoupling = trench.d_decoupling,
        )
        P = vec(Phase)   # one wrapper on both sides, so that the broadcast needs no copy
        P .= weakzone_phase.(Ref(pw), P, vec(X), vec(Y), vec(Z))
    end

    return nothing
end

# Coordinates in the frame of the trench, rotated along the angle from Start -> End
rotated(p, x, y, z) = Rot3D(x - p.Start[1], y - p.Start[2], z, p.cosStrike, p.sinStrike, 1.0, -0.0)

function slab_coordinates(p, x, y, z)
    ls, d = slab_distances(rotated(p, x, y, z)..., p.xb, p.segments, p.bbox)
    return -p.Thickness <= d <= 0.0 ? (ls, y, d) : nothing
end

function weakzone_phase(p, ph, x, y, z)
    xr, yr, zr = rotated(p, x, y, z)
    d = slab_distances(xr, yr, zr, p.xb, p.segments, p.bbox)[2]
    return -p.WeakzoneThickness <= d <= 0.0 && zr > -p.d_decoupling ? p.WeakzonePhase : ph
end

"""
    add_fault!(Phase, Temp, Grid::AbstractGeneralGrid;
        Start=(20,100), End=(10,80),
        Fault_thickness=10.0,
        Depth_extent=nothing,
        DipAngle=0e0,
        phase=ConstantPhase(1),
        T=nothing,
        cell=false)

Adds a fault to the given 3D grid by modifying the `Phase` and `Temp` arrays.
For a 2D grid, use `add_box` instead.

# Arguments
- `Phase`: Phase array
- `Temp`: Temp array
- `Grid`: The grid on which the fault is to be added.
- `Start`: Tuple representing the starting coordinates of the fault (X, Y).
- `End`: Tuple representing the ending coordinates of the fault (X, Y).
- `Fault_thickness`: Thickness of the fault.
- `Depth_extent`: Depth extent of the fault. If `nothing`, the fault extends through the entire domain.
- `DipAngle`: Dip angle of the fault.
- `phase`: Phase to be assigned to the fault.
- `T`: Temperature to be assigned to the fault. If `nothing`, the temperature is not modified.


# Example
```julia
add_fault!(Phase, Temp, Grid;
        Start=(20,100), End=(10,80),
        Fault_thickness=10.0,
        Depth_extent=(-25.0, 0.0),
        DipAngle=-10.0,
        phase=ConstantPhase(1)
        )
```
"""
function add_fault!(
        Phase,
        Temp,
        Grid::AbstractGeneralGrid;
        Start = (20, 100),
        End = (10, 80),
        Fault_thickness = 10.0,
        Depth_extent = nothing,
        DipAngle = 0.0e0,
        phase = ConstantPhase(1),
        T = nothing,
        cell = false
    )

    # Extract the coordinates
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Calculate the direction vector from Start to End
    direction = (End[1] - Start[1], End[2] - Start[2])
    length = sqrt(direction[1]^2 + direction[2]^2)
    unit_direction = (direction[1], direction[2]) ./ length

    # Calculate the fault region based on fault thickness and length
    fault_half_thickness = Fault_thickness / 2

    # The fault plane contains the strike line through Start and dips about it, towards the
    # horizontal normal (unit_direction[2], -unit_direction[1]) for positive DipAngle.
    sinDip, cosDip = sincosd(DipAngle)
    p = (; Start, unit_direction, length, fault_half_thickness, cosDip, sinDip, Depth_extent)

    # without a depth extent, the fault spans the whole grid vertically; see #135
    zbot, ztop = isnothing(Depth_extent) ? finite_extrema(Z) : extrema(Depth_extent)
    (xmin, xmax), (ymin, ymax) = extrema((Start[1], End[1])), extrema((Start[2], End[2]))
    frame = (; xmin, xmax, ymin, ymax, zbot, ztop)
    set_shape!(Phase, Temp, X, Y, Z, fault_coordinates, p, frame, phase, T)

    return nothing
end

function in_fault(p, x, y, z)
    dx = x - p.Start[1]
    dy = y - p.Start[2]
    projection_length = dx * p.unit_direction[1] + dy * p.unit_direction[2]
    horizontal_offset = dx * p.unit_direction[2] - dy * p.unit_direction[1]
    # distance to the dipping plane (Z = 0 at the surface)
    perpendicular_distance = abs(horizontal_offset * p.cosDip + z * p.sinDip)
    return 0 ≤ projection_length ≤ p.length && perpendicular_distance ≤ p.fault_half_thickness
end

in_depth(p, z) = isnothing(p.Depth_extent) || p.Depth_extent[1] ≤ z ≤ p.Depth_extent[2]

fault_coordinates(p, x, y, z) = in_fault(p, x, y, z) && in_depth(p, z) ? (x, y, z) : nothing
