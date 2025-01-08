using GeophysicalModelGenerator

# Grid parameters
nx, ny, nz = 128, 128, 57
x = range(-1000, 1000, length=nx)
y = range(-1000, 1000, length=ny)
z = range(-660, 0, length=nz)
Grid = CartData(xyz_grid(x, y, z))

# Phases and temperature
Phases = fill(2, nx, ny, nz)
Temp = fill(1350.0, nx, ny, nz)

# Ridge Segments
segments = [
    ((-500.0, -1000.0), (-500.0, 0.0)),  # Segment 1
    ((-250.0, 0.0), (-250.0, 200.0)),    # Segment 2
    ((-750.0, 200.0), (-750.0, 1000.0))  # Segment 3
]

# Create a vector of SpreadingRateTemp parameters (one for each segment)
s = [
    SpreadingRateTemp(Tsurface=0.0, Tmantle=1350.0, Adiabat=0.0, SpreadingVel=3.0, AgeRidge=0.0, maxAge=100.0),  # Segmento 1
    SpreadingRateTemp(Tsurface=0.0, Tmantle=1350.0, Adiabat=0.0, SpreadingVel=1.0, AgeRidge=0.0, maxAge=100.0),  # Segmento 2
    SpreadingRateTemp(Tsurface=0.0, Tmantle=1350.0, Adiabat=0.0, SpreadingVel=6.0, AgeRidge=0.0, maxAge=100.0)   # Segmento 3
]

# Lithospheric phases and temperature setup
lith = LithosphericPhases(Layers=[15, 55], Phases=[1, 2], Tlab=1250)

# Add a box for the test (adjust parameters as needed)
add_box!(Phases, Temp, Grid; xlim=(-1000.0,0.0), ylim=(-500.0, 500.0), zlim=(-80.0, 0.0), phase=lith, 
         T=s, segments=segments)

# Add and save results
Grid = addfield(Grid, (; Phases, Temp))
write_paraview(Grid, "Ridge_Thermal_Structure_test_2")


################3
# Case 3: Multiple segments with different parameters

function compute_thermal_structure(Temp, X, Y, Z, Phase, s::Vector{SpreadingRateTemp}, segments::Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}})
    @show "Using compute_thermal_structure with Vector{SpreadingRateTemp_3}"
    kappa = 1e-6
    SecYear = 3600 * 24 * 365
    dz = Z[end] - Z[1]

    # MantleAdiabaticT debe calcularse fuera del bucle para ser consistente en todas las regiones
    MantleAdiabaticT = [params.Tmantle + params.Adiabat * abs.(Z) for params in s]

    # Crear delimitadores
    delimiters = [(segments[i][2], segments[i + 1][1]) for i in 1:length(segments) - 1]

    # Iterar sobre cada índice de X
    for I in eachindex(X)
        px, py, pz = X[I], Y[I], Z[I]

        # Determinar la región de este punto
        region = determine_region(px, py, delimiters, segments)

        # Seleccionar el segmento correspondiente y sus parámetros
        segment = segments[region]
        params = s[region]  # Parámetros correspondientes al segmento

        # Desestructuración de los parámetros
        @unpack Tsurface, Tmantle, Adiabat, SpreadingVel, AgeRidge, maxAge = params

        # Calcular MantleAdiabaticT para este segmento (dentro del bucle)
        MantleAdiabaticT = Tmantle .+ Adiabat .* abs.(pz)

        x1, y1 = segment[1]
        x2, y2 = segment[2]

        # Calcular la distancia al segmento
        Distance = perpendicular_distance_to_segment(px, py, x1, y1, x2, y2)

        # Calcular la edad térmica
        ThermalAge = abs(Distance * 1e3 * 1e2) / SpreadingVel + AgeRidge * 1e6  # Edad térmica en años
        if ThermalAge > maxAge * 1e6
            ThermalAge = maxAge * 1e6
        end

        ThermalAge = ThermalAge * SecYear  # Convertir a segundos
        if ThermalAge == 0
            ThermalAge = 1e-6  # Evitar cero
        end

        # Calcular temperatura
        Temp[I] = (Tsurface - Tmantle) * erfc(abs(pz) * 1e3 / (2 * sqrt(kappa * ThermalAge))) + MantleAdiabaticT[I]
    end

    return Temp
end




if segments !== nothing
    if isa(T, Vector{SpreadingRateTemp})
        Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T, segments)
    else
        Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], X[ind], Y[ind], Z[ind], Phase[ind_flat], T, segments)
    end
else
    Temp[ind_flat] = compute_thermal_structure(Temp[ind_flat], Xrot[ind], Yrot[ind], Zrot[ind], Phase[ind_flat], T)
end