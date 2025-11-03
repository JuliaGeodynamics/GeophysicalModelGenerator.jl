using Base: Int64, Float64, NamedTuple
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

"""
    ind2D = flatten_index_dimensions(Phase, ind_vec::Vector{CartesianIndex{3}})

This converts the indices to purely 2D indices if the array `phase` is 2D
"""
function flatten_index_dimensions(Phase, ind_vec::Vector{CartesianIndex{3}})
    if length(size(Phase)) == 2
        ind2D = Vector{CartesianIndex{2}}(undef, length(ind_vec))
        for (num, ind) in enumerate(ind_vec)
            ind2D[num] = CartesianIndex(ind[1], ind[3])
        end
    else
        ind2D = ind_vec
    end

    return ind2D
end

"""
    ind2D = flatten_index_dimensions(Phase, ind_vec::Vector{CartesianIndex{3}})

This converts the indices to purely 2D indices if the array `phase` is 2D
"""
function flatten_index_dimensions(Phase::AbstractArray{T, N}, ind_vec::Array{Bool, 3}) where {T, N}
    if N == 2
        ind2D = Vector{CartesianIndex{2}}(undef, length(ind_vec))
        for (num, ind) in enumerate(ind_vec)
            ind2D[num] = CartesianIndex(ind[1], ind[3])
        end
    else
        ind2D = ind_vec
    end

    return ind2D
end

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
        Origin = (maximum(X), maximum(Y), maximum(Z))  # upper-left corner
    end

    # Perform rotation of 3D coordinates:
    Xrot = X .- Origin[1]
    Yrot = Y .- Origin[2]
    Zrot = Z .- Origin[3]

    Rot3D!(Xrot, Yrot, Zrot, StrikeAngle, DipAngle)

    ph_ind = findall(Phase .== phase.phase)

    ind = Int64[]
    if stripAxes[1] == 1
        indX = findall(abs.(Xrot[ph_ind] .% stripeSpacing) .<= stripeWidth / 2.0)
        ind = vcat(ind, indX)
    end
    if stripAxes[2] == 1
        indY = findall(abs.(Yrot[ph_ind] .% stripeSpacing) .<= stripeWidth / 2.0)
        ind = vcat(ind, indY)
    end
    if stripAxes[3] == 1
        indZ = findall(abs.(Zrot[ph_ind] .% stripeSpacing) .<= stripeWidth / 2.0)
        ind = vcat(ind, indZ)
    end

    Phase[ph_ind[ind]] .= stripePhase.phase

    return nothing
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
    if ylim == nothing
        ylim = (minimum(Y), maximum(Y))
    end

    if Origin == nothing
        Origin = (xlim[1], ylim[1], zlim[2])  # upper-left corner
    end

    if Origin !== nothing && isa(T, McKenzie_subducting_slab)
        @warn  "McKenzie temperature does not require the definition of 'Origin' field; if Origin is defined it must be equal to [xmin,ymin,zmax] of the box that has been defined."
        if Origin[1] != xlim[1] || Origin[2] != ylim[1] || Origin[3] != zlim[2]
            @error  "Origin is not set up correctly. For fixing the problem Origin can be left blank or Origin = [xmin,ymin,zmax] of the box"
        end
    end

    # Perform rotation of 3D coordinates:
    Xrot = X .- Origin[1]
    Yrot = Y .- Origin[2]
    Zrot = Z .- Origin[3]

    Rot3D!(Xrot, Yrot, Zrot, StrikeAngle, DipAngle)

    # Set phase number & thermal structure in the full domain
    ztop = maximum(zlim) - Origin[3]
    zbot = minimum(zlim) - Origin[3]
    ind = findall(
        (Xrot .>= (minimum(xlim) - Origin[1])) .& (Xrot .<= (maximum(xlim) - Origin[1])) .&
            (Yrot .>= (minimum(ylim) - Origin[2])) .& (Yrot .<= (maximum(ylim) - Origin[2])) .&
            (Zrot .>= zbot) .& (Zrot .<= ztop)
    )

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly. See routines below for different options
        if T != nothing
            if isa(T, LithosphericTemp)
                Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], phase)
            end
            if segments !== nothing
                Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], Phase[ind_flat], T, segments)
            else
                Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], Phase[ind_flat], T)
            end
        end
        # Set the phase. Different routines are available for that - see below.
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], phase)
    end

    return nothing
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
    if isnothing(xlim) == isnothing(ylim) == isnothing(zlim)
        error("You need to specify at least one of the limits (xlim, ylim, zlim)")
    end

    if isnothing(xlim)
        xlim = (minimum(X), maximum(X))
    end
    if isnothing(ylim)
        ylim = (minimum(Y), maximum(Y))
    end
    if isnothing(zlim)
        zlim = (minimum(Z), maximum(Z))
    end

    # Set phase number & thermal structure in the full domain
    ind = findall(
        (X .>= (xlim[1])) .& (X .<= (xlim[2])) .&
            (Y .>= (ylim[1])) .& (Y .<= (ylim[2])) .&
            (Z .>= (zlim[1])) .& (Z .<= (zlim[2]))
    )

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly. See routines below for different options
        if !isnothing(T)
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T)
        end

        # Set the phase. Different routines are available for that - see below.
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end


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
    ind = findall(((X .- cen[1]) .^ 2 + (Y .- cen[2]) .^ 2 + (Z .- cen[3]) .^ 2) .^ 0.5 .< radius)

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly. See routines below for different options
        if T != nothing
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T)
        end

        # Set the phase. Different routines are available for that - see below.
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end

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

    if Origin == nothing
        Origin = cen  # center
    end

    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    # Perform rotation of 3D coordinates:
    Xrot = X .- Origin[1]
    Yrot = Y .- Origin[2]
    Zrot = Z .- Origin[3]

    Rot3D!(Xrot, Yrot, Zrot, StrikeAngle, DipAngle)

    # Set phase number & thermal structure in the full domain
    x2 = axes[1]^2
    y2 = axes[2]^2
    z2 = axes[3]^2
    cenRot = cen .- Origin
    ind = findall(
        (
            ((Xrot .- cenRot[1]) .^ 2) ./ x2 + ((Yrot .- cenRot[2]) .^ 2) ./ y2 +
                ((Zrot .- cenRot[3]) .^ 2) ./ z2
        ) .^ 0.5 .<= 1
    )

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly. See routines below for different options
        if T != nothing
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], Phase[ind_flat], T)
        end

        # Set the phase. Different routines are available for that - see below.
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], phase)
    end

    return nothing
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

    # distance between grid points and cylinder base
    dx_b = X .- base[1]
    dy_b = Y .- base[2]
    dz_b = Z .- base[3]

    # find normalized parametric coordinate of a point-axis projection
    t = (axVec[1] .* dx_b .+ axVec[2] .* dy_b .+ axVec[3] .* dz_b) ./ ax2

    # find distance vector between point and axis
    dx = dx_b .- t .* axVec[1]
    dy = dy_b .- t .* axVec[2]
    dz = dz_b .- t .* axVec[3]

    # Set phase number & thermal structure in the full domain
    ind = findall((t .>= 0.0) .& (t .<= 1.0) .& ((dx .^ 2 + dy .^ 2 + dz .^ 2) .^ 0.5 .<= radius))

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly. See routines below for different options
        if !isnothing(T)
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T)
        end

        # Set the phase. Different routines are available for that - see below.
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end

# Internal function that rotates the coordinates
function Rot3D!(X, Y, Z, StrikeAngle, DipAngle)

    # precompute trigonometric functions (expensive!)
    sindStrikeAngle, cosStrikeAngle = sincosd(StrikeAngle)
    sinDipAngle, cosDipAngle = sincosd(-DipAngle)   # note the minus here to be consistent with the earlier version of the code
    for i in eachindex(X)
        X[i], Y[i], Z[i] = Rot3D(X[i], Y[i], Z[i], cosStrikeAngle, sindStrikeAngle, cosDipAngle, sinDipAngle)
    end

    return nothing
end


"""
        add_polygon!(Phase, Temp, Grid::AbstractGeneralGrid; xlim=(), ylim::Tuple = (0.0,0.8), zlim=(), phase = ConstantPhase(1), T=nothing, cell=false )


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
        xlim = (), ylim::Tuple = (0.0, 0.8), zlim = (),           # limits of the box
        phase = ConstantPhase(1),                                   # Sets the phase number(s) in the box
        T = nothing, cell = false
    )                             # Sets the thermal structure (various functions are available)


    xlim_ = Float64.(collect(xlim))
    ylim_ = Float64.(collect(ylim))
    zlim_ = Float64.(collect(zlim))


    # Retrieve 3D data arrays for the grid
    X, Y, Z = coordinate_grids(Grid, cell = cell)

    ind = zeros(Bool, size(X))
    ind_slice = zeros(Bool, size(X[:, 1, :]))

    # find points within the polygon, only in 2D
    for i in 1:size(Y)[2]
        if Y[1, i, 1] >= ylim_[1] && Y[1, i, 1] <= ylim_[2]
            inpolygon!(ind_slice, xlim_, zlim_, X[:, i, :], Z[:, i, :])
            ind[:, i, :] = ind_slice
        else
            ind[:, i, :] = zeros(size(X[:, 1, :]))
        end
    end

    if !isempty(ind)
        # Compute thermal structure accordingly. See routines below for different options
        if T != nothing
            Temp[ind] = compute_thermal_structure(Temp[ind], X[ind], Y[ind], Z[ind], Phase[ind], T)
        end

        # Set the phase. Different routines are available for that - see below.
        Phase[ind] = compute_phase(Phase[ind], Temp[ind], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
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
    ind = zeros(Bool, size(X))
    ind_slice = zeros(Bool, size(X[:, :, 1]))

    for k in 1:size(Z, 3)
        if zlim_[1] <= Z[1, 1, k] <= zlim_[2]
            inpolygon!(ind_slice, xlim_, ylim_, X[:, :, k], Y[:, :, k])
            @views ind[:, :, k] = ind_slice
        else
            @views ind[:, :, k] = zeros(size(X[:, :, 1]))
        end
    end

    if !isempty(ind)
        if T != nothing
            if segments !== nothing
                Temp[ind] = compute_thermal_structure(Temp[ind], X[ind], Y[ind], Z[ind], Phase[ind], T, segments)
            else
                Temp[ind] = compute_thermal_structure(Temp[ind], X[ind], Y[ind], Z[ind], Phase[ind], T)
            end
        end
        Phase[ind] = compute_phase(Phase[ind], Temp[ind], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end

"""
    xrot, yrot, zrot = Rot3D(X::Number,Y::Number,Z::Number, cosStrikeAngle, sindStrikeAngle, cosDipAngle, sinDipAngle)

Perform rotation for a point in 3D space
"""
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
- T - temperature structure of the volcano
- crater - this will create a truncated cone and the option defines the radius of the flat top
- base - this sets the flat topography around the volcano
- background - this allows loading in a topography and only adding the volcano on top (also allows stacking of several cones to get a volcano with different slopes)
"""
function add_volcano!(
        Phases,
        Temp,
        Grid::CartData;
        volcanic_phase = 1,
        center = (0, 0, 0),
        height = 0.0,
        radius = 0.0,
        crater = 0.0,
        base = 0.0,
        background = nothing,
        T = HalfspaceCoolingTemp(Age = 0)
    )
    H = make_volc_topo(
        Grid;
        center = center,
        height = height,
        radius = radius,
        crater = crater,
        base = base,
        background = background
    )

    ni = size(Grid.x)
    ind = fill(false, ni...)
    depth = similar(Grid.z.val)

    for k in axes(ind, 3)
        for j in axes(ind, 2), i in axes(ind, 1)
            depth[i, j, k] = max(H[i, j] - Grid.z.val[i, j, k], 0)

            if Grid.z.val[i, j, k] < H[i, j] && Grid.z.val[i, j, k] ≥ base
                Phases[i, j, k] = volcanic_phase
            end
            if Phases[i, j, k] > 0
                ind[i, j, k] = true
            end
        end
    end

    ind_flat = flatten_index_dimensions(Phases, ind)

    # @views Temp[ind .== false] .= 0.0
    if !isempty(ind_flat)
        if !isnothing(T)
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], Grid.x.val[ind], Grid.y.val[ind], depth[ind], Phases[ind_flat], T)
        end
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
    pos = similar(X)

    for i in eachindex(pos)
        # compute radial distance to volcano center
        DX = X[i] - center[1]
        DY = Y[i] - center[2]
        RD = √(DX^2 + DY^2)

        # get radial distance from crater rim
        RD -= crater

        # find position relative to crater rim
        dr = radius - crater
        pos[i] = -RD / dr + 1
    end

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

    for i in eachindex(pos)
        if 0 ≤ pos[i] < 1
            H[i] = pos[i] * (height - base) + base
        elseif pos[i] ≥ 1
            H[i] = height
        end
    end

    return H
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

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::ConstantTemp)
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

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::LinearTemp)
    @unpack Ttop, Tbot = s

    dz = Z[end] - Z[1]
    dT = Tbot - Ttop

    Temp = abs.(Z ./ dz) .* dT .+ Ttop
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
    Age = 60          # thermal age of plate [in Myrs]
    Adiabat = 0        # Adiabatic gradient in K/km
end

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::HalfspaceCoolingTemp)
    @unpack Tsurface, Tmantle, Age, Adiabat = s

    kappa = 1.0e-6
    SecYear = 3600 * 24 * 365
    dz = Z[end] - Z[1]
    ThermalAge = Age * 1.0e6 * SecYear

    MantleAdiabaticT = Tmantle .+ Adiabat * abs.(Z)    # Adiabatic temperature of mantle

    for i in eachindex(Temp)
        Temp[i] = (Tsurface .- Tmantle) * erfc((abs.(Z[i]) * 1.0e3) ./ (2 * sqrt(kappa * ThermalAge))) + MantleAdiabaticT[i]
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

Note: the thermal age at the mid oceanic ridge is set to 1 year to avoid division by zero

"""
@with_kw_noshow mutable struct SpreadingRateTemp <: AbstractThermalStructure
    Tsurface = 0       # top T
    Tmantle = 1350     # bottom T
    Adiabat = 0        # Adiabatic gradient in K/km
    MORside = "left"   # side of box where the MOR is located
    SpreadingVel = 3   # spreading velocity [cm/yr]
    AgeRidge = 0       # Age of the ridge [Myrs]
    maxAge = 60       # maximum thermal age of plate [Myrs]
end

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::SpreadingRateTemp)
    @unpack Tsurface, Tmantle, Adiabat, MORside, SpreadingVel, AgeRidge, maxAge = s

    kappa = 1.0e-6
    SecYear = 3600 * 24 * 365
    dz = Z[end] - Z[1]

    MantleAdiabaticT = Tmantle .+ Adiabat * abs.(Z)    # Adiabatic temperature of mantle

    if MORside == "left"
        Distance = X .- X[1, 1, 1]
    elseif MORside == "right"
        Distance = X[end, 1, 1] .- X
    elseif MORside == "front"
        Distance = Y .- Y[1, 1, 1]
    elseif MORside == "back"
        Distance = Y[1, end, 1] .- Y

    else
        error("unknown side")
    end

    for i in eachindex(Temp)
        ThermalAge = abs(Distance[i] * 1.0e3 * 1.0e2) / SpreadingVel + AgeRidge * 1.0e6    # Thermal age in years
        if ThermalAge > maxAge * 1.0e6
            ThermalAge = maxAge * 1.0e6
        end

        ThermalAge = ThermalAge * SecYear
        if ThermalAge == 0
            ThermalAge = 1.0e-6   # doesn't like zero
        end

        Temp[i] = (Tsurface .- Tmantle) * erfc((abs.(Z[i]) * 1.0e3) ./ (2 * sqrt(kappa * ThermalAge))) + MantleAdiabaticT[i]

    end

    return Temp
end

"""
    SpreadingRateTemp(Temp, X, Y, Z, Phase, s::SpreadingRateTemp, segments)

Calculates the temperature distribution across the plate considering multiple ridge segments.

This function computes the thermal structure based on the perpendicular distance from each point to its corresponding ridge segment, and applies a thermal model using a spreading velocity and thermal age.

Parameters
==========
- Temp     : Temperature field to be updated (array)
- X, Y, Z  : Coordinates of the points (arrays)
- Phase    : Phase of the material (unused in this version)
- s        : SpreadingRateTemp object containing the thermal and spreading parameters
- segments : List of ridge segments, where each segment is defined by two tuples representing the start and end coordinates (x1, y1) and (x2, y2) for each segment.

Note
====
The temperature at each point is calculated using the thermal age, which is determined by the distance from the point to the nearest ridge segment and the spreading velocity.

The function works in the context of one or more segments. The key difference from the previous function is that the ridge can now be placed at any position within the box, not just at the boundary.

The thermal age is capped at `maxAge` years, and the temperature is adjusted based on the distance to the ridge and the corresponding thermal gradient.

"""

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::SpreadingRateTemp, segments::Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}})
    @unpack Tsurface, Tmantle, Adiabat, SpreadingVel, AgeRidge, maxAge = s
    kappa = 1.0e-6
    SecYear = 3600 * 24 * 365
    dz = Z[end] - Z[1]

    MantleAdiabaticT = Tmantle .+ Adiabat * abs.(Z)

    #Create delimiters
    delimiters = [(segments[i][2], segments[i + 1][1]) for i in 1:(length(segments) - 1)]

    for I in eachindex(X)
        px, py, pz = X[I], Y[I], Z[I]

        # Determine region of point
        region = determine_region(px, py, delimiters, segments)

        # Select the corresponding segment
        x1, y1 = segments[region][1]
        x2, y2 = segments[region][2]

        # Calculate distance to segment
        Distance = perpendicular_distance_to_segment(px, py, x1, y1, x2, y2)

        # Calculate thermal age
        ThermalAge = abs(Distance * 1.0e5) / SpreadingVel + AgeRidge * 1.0e6  # Thermal age in years
        if ThermalAge > maxAge * 1.0e6
            ThermalAge = maxAge * 1.0e6
        end

        ThermalAge = ThermalAge * SecYear  # Convert to seconds
        if ThermalAge == 0
            ThermalAge = 1.0e-6  # Avoid zero
        end

        # Calculate temperature
        Temp[I] = (Tsurface - Tmantle) * erfc(abs(pz) * 1.0e3 / (2 * sqrt(kappa * ThermalAge))) + MantleAdiabaticT[I]
    end

    return Temp

end

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
@with_kw_noshow mutable struct LithosphericTemp <: AbstractThermalStructure
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

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::LithosphericTemp)
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

    interp_linear_T = linear_interpolation(-z ./ 1.0e3, T .- 273.15)      # create interpolation object
    Temp = interp_linear_T(-Z)

    return Temp
end

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

function compute_phase(Phase, Temp, X, Y, Z, s::ConstantPhase)
    Phase .= s.phase
    return Phase
end


"""
    LithosphericPhases(Layers=[10 20 15], Phases=[1 2 3 4], Tlab=nothing )

This allows defining a layered lithosphere. Layering is defined from the top downwards.

Parameters
===
- Layers : The thickness of each layer, ordered from top to bottom. The thickness of the last layer does not have to be specified.
- Phases : The phases of the layers, ordered from top to bottom.
- Tlab   : Temperature of the lithosphere asthenosphere boundary. If specified, the phases at locations with T>Tlab are set to Phases[end].

"""
@with_kw_noshow mutable struct LithosphericPhases <: AbstractPhaseNumber
    Layers = [10.0, 20.0, 15.0]
    Phases = [1, 2, 3, 4]
    Tlab = nothing
end


"""
    Phase = compute_phase(Phase, Temp, X, Y, Z, s::LithosphericPhases, Ztop)

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
function compute_phase(Phase, Temp, X, Y, Z, s::LithosphericPhases; Ztop = 0)
    @unpack Layers, Phases, Tlab = s

    Phase .= Phases[end]

    for i in 1:length(Layers)
        Zbot = Ztop - Layers[i]
        ind = findall((Z .>= Zbot) .& (Z .<= Ztop))
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
@with_kw_noshow mutable struct McKenzie_subducting_slab <: AbstractThermalStructure
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
of McKenzie (1969) ["Speculations on the consequences and causes of plate motions"]. The functions assumes
that the bottom of the slab is the coordinate Z=0. Internally the function shifts the coordinate.

Parameters

=============================
- `Temp`:  Temperature array
- `X`:    X Array
- `Y`:    Y Array
- `Z`:    Z Array
- `Phase`: Phase array
- `s`:    `McKenzie_subducting_slab`
"""
function compute_thermal_structure(Temp, X, Y, Z, Phase, s::McKenzie_subducting_slab)
    @unpack Tsurface, Tmantle, Adiabat, v_cm_yr, κ, it = s

    # Thickness of the layer:
    Thickness = (maximum(Z) - minimum(Z))
    Zshift = Z .- Z[end]       # McKenzie model is defined with Z = 0 at the bottom of the slab

    # Convert subduction velocity from cm/yr -> m/s;
    convert_velocity = 1 / (100.0 * 365.25 * 60.0 * 60.0 * 24.0)
    v_s = v_cm_yr * convert_velocity

    # calculate the thermal Reynolds number
    Re = (v_s * Thickness * 1000) / 2 / κ      # factor 1000 to transfer Thickness from km to m

    # McKenzie model
    sc = 1 / Thickness
    σ = ones(size(Temp))
    # Dividi et impera
    for i in 1:it
        a = (-1.0) .^ (i) ./ (i .* pi)
        b = (Re .- (Re .^ 2 .+ i^2.0 .* pi^2.0) .^ (0.5)) .* X .* sc
        c = sin.(i .* pi .* (1 .- abs.(Zshift .* sc)))
        e = exp.(b)
        σ .+= 2 * a .* e .* c
    end

    Temp .= Tsurface .+ (Tmantle - Tsurface) .* σ
    Temp .= Temp + (Adiabat * abs.(Z))

    return Temp
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
@with_kw_noshow mutable struct LinearWeightedTemperature <: AbstractThermalStructure
    w_min::Float64 = 0.0
    w_max::Float64 = 1.0
    crit_dist::Float64 = 100.0
    dir::Symbol = :X
    F1::AbstractThermalStructure = ConstantTemp()
    F2::AbstractThermalStructure = ConstantTemp()
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
function compute_thermal_structure(Temp, X, Y, Z, Phase, s::LinearWeightedTemperature)
    @unpack w_min, w_max, crit_dist, dir = s
    @unpack F1, F2 = s

    if dir === :X
        dist = X
    elseif dir === :Y
        dist = Y
    else
        dist = Z
    end

    # compute the 1D thermal structures
    Temp1 = zeros(size(Temp))
    Temp2 = zeros(size(Temp))
    Temp1 = compute_thermal_structure(Temp1, X, Y, Z, Phase, F1)
    Temp2 = compute_thermal_structure(Temp2, X, Y, Z, Phase, F2)

    # Compute the weights
    weight = w_min .+ (w_max - w_min) ./ (crit_dist) .* (dist)

    ind_1 = findall(weight .> w_max)
    ind_2 = findall(weight .< w_min)

    # Change the weight
    weight[ind_1] .= w_max
    weight[ind_2] .= w_min

    # Average temperature
    Temp .= Temp1 .* (1.0 .- weight) + Temp2 .* weight

    return Temp
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
@with_kw_noshow mutable struct Trench{Nseg} <: AbstractTrenchSlab
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
    l = 0.0       # initial length
    it = 1         # iteration

    dl = Length / n_seg  # dl
    while l < Length

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

        # update l
        l = l + dl
        it = it + 1
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

"""
    find_slab_distance!(ls, d, X,Y,Z, trench::Trench)

Function that finds the perpendicular distance to the top and bottom of the slab `d`, and the current length of the slab `l`.

"""
function find_slab_distance!(ls, d, X, Y, Z, Top, Bottom, trench::Trench)
    @unpack Thickness, Length, n_seg, Start, End, direction = trench

    # Perform rotation of 3D coordinates along the angle from Start -> End:
    Xrot = X .- Start[1]
    Yrot = Y .- Start[2]

    StrikeAngle = -atand((End[2] - Start[2]) / (End[1] - Start[1]))
    Rot3D!(Xrot, Yrot, Z, StrikeAngle, 0.0)

    xb = Rot3D(End[1] - Start[1], End[2] - Start[2], 0.0, cosd(StrikeAngle), sind(StrikeAngle), 1.0, 0.0)

    # dl
    dl = trench.Length / n_seg
    l = 0  # length at the trench position
    D = @SVector [Top[1, 2], Bottom[1, 2], Bottom[1, 2], Top[1, 2]]

    # Construct the slab
    for i in 1:(n_seg - 1)
        ln = l + dl

        pa = (Top[i, 1], Top[i, 2])        # D = 0 | L = l
        pb = (Bottom[i, 1], Bottom[i, 2])  # D = -Thickness | L=l

        pc = (Bottom[i + 1, 1], Bottom[i + 1, 2])  # D = -Thickness |L=L+dl
        pd = (Top[i + 1, 1], Top[i + 1, 2]) # D = 0| L = L+dl

        # Create the polygon
        poly_y = @SVector [pa[1], pb[1], pc[1], pd[1]]
        poly_z = @SVector [pa[2], pb[2], pc[2], pd[2]]

        # find a sub set of particles
        ymin, ymax = extrema(poly_y)
        zmin, zmax = extrema(poly_z)

        ind_s = findall(0.0 .<= Xrot .<= xb[1] .&& ymin .<= Yrot .<= ymax .&& zmin .<= Z .<= zmax)

        # Find the particles
        yp = Yrot[ind_s]
        zp = Z[ind_s]

        # Initialize the ind that are going to be used by inpoly
        ind = zeros(Bool, size(zp))
        inpolygon!(ind, poly_y, poly_z, yp, zp)         # determine whether points are inside the polygon or not

        # indexes of the segment
        ind_seg = ind_s[ind]

        # Loop over the chosen particles and interpolate the current value of L and D.
        for ip in ind_seg
            point_ = (Yrot[ip], Z[ip])
            d[ip] = -distance_to_linesegment(point_, pa, pd)
            ls[ip] = distance_to_linesegment(point_, pb, pa) + l
        end

        #Update l
        l = ln
    end
    return
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

    # Find the distance to the slab (along & perpendicular)
    d = fill(NaN, size(Grid))        # -> d = distance perpendicular to the slab
    ls = fill(NaN, size(Grid))       # -> l = length from the trench along the slab
    find_slab_distance!(ls, d, X, Y, Z, Top, Bottom, trench)

    # Function to fill up the temperature and the phase.
    ind = findall((-trench.Thickness .<= d .<= 0.0))

    if !isempty(ind)
        if isa(T, LinearWeightedTemperature)
            l_decouplingind = findall(Top[:, 2] .<= -trench.d_decoupling)
            if !isempty(l_decouplingind)
                l_decoupling = Top[l_decouplingind[1], 1]
                T.crit_dist = abs(l_decoupling)
            end
        end

        # Compute thermal structure accordingly. See routines below for different options {Future: introducing the length along the trench for having lateral varying properties along the trench}
        if !isnothing(T)
            Temp[ind] = compute_thermal_structure(Temp[ind], ls[ind], Y[ind], d[ind], Phase[ind], T)
        end

        # Set the phase
        Phase[ind] = compute_phase(Phase[ind], Temp[ind], ls[ind], Y[ind], d[ind], phase)

        # Add a weak zone on top of the slab (indicated by a phase number but not by temperature)
        if trench.WeakzoneThickness > 0.0
            d_weakzone = fill(NaN, size(Grid))        # -> d = distance perpendicular to the slab
            ls_weakzone = fill(NaN, size(Grid))       # -> l = length from the trench along the slab
            find_slab_distance!(ls_weakzone, d_weakzone, X, Y, Z, WeakZone, Top, trench)

            ind = findall((-trench.WeakzoneThickness .<= d_weakzone .<= 0.0) .& (Z .> -trench.d_decoupling))
            Phase[ind] .= trench.WeakzonePhase
        end
    end

    return nothing
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

    # Create a mask for the fault region
    fault_mask = falses(size(X))

    for k in 1:size(Z, 3), j in 1:size(Y, 2), i in 1:size(X, 1)
        # Rotate the point using the dip angle
        x_rot, y_rot, z_rot = Rot3D(X[i, j, k], Y[i, j, k], Z[i, j, k], 1.0, 0.0, cosd(DipAngle), sind(DipAngle))

        # Calculate the projection of the rotated point onto the fault line
        projection_length = (x_rot - Start[1]) * unit_direction[1] + (y_rot - Start[2]) * unit_direction[2]
        if 0 ≤ projection_length ≤ length
            # Calculate the perpendicular distance to the fault line
            perpendicular_distance = abs((x_rot - Start[1]) * unit_direction[2] - (y_rot - Start[2]) * unit_direction[1])
            if perpendicular_distance ≤ fault_half_thickness
                fault_mask[i, j, k] = true
            end
        end
    end

    ind = findall(fault_mask)

    # Apply depth extent if provided
    if !isnothing(Depth_extent)
        ind = ind[Z[ind] .≥ Depth_extent[1] .&& Z[ind] .≤ Depth_extent[2]]
    end

    ind_flat = flatten_index_dimensions(Phase, ind)

    if !isempty(ind_flat)
        # Compute thermal structure accordingly
        if T != nothing
            Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T)
        end

        # Set the phase
        Phase[ind_flat] = compute_phase(Phase[ind_flat], Temp[ind_flat], X[ind], Y[ind], Z[ind], phase)
    end

    return nothing
end
