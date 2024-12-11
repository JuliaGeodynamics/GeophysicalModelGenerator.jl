# The script allows for modeling complex oceanic ridge configurations, supporting multiple ridge segments 
# and enabling the calculation of temperature variations based on ridge spreading, age, and thermal properties.

using GeophysicalModelGenerator   
using SpecialFunctions: erfc

# Grid parameters
nx, ny, nz = 512, 512, 128
x = range(-1000, 1000, length=nx)
y = range(-1000, 1000, length=ny)
z = range(-660, 0, length=nz)
Grid = CartData(xyz_grid(x, y, z))

# Phases and temperature
Phases = fill(2, nx, ny, nz)
Temp = fill(1350.0, nx, ny, nz)

# Ridge segments
segments = [
    ((-500.0, -1000.0), (-500.0, 0.0)),  # Segment 1
    ((-250.0, 0.0), (-250.0, 200.0)),    # Segment 2
    ((-750.0, 200.0), (-750.0, 1000.0))  # Segment 3  #It is possible to add more segments
]

# Create delimiters for each segment
# This ensures that each region is influenced by only one ridge segment
delimiters = [(segments[i][2], segments[i + 1][1]) for i in 1:length(segments) - 1]

# Ridge parameters
Tsurface = 0.0
Tmantle = 1350.0
Adiabat = 0.0
SpreadingVel = 3.0  # cm/yr
AgeRidge = 0.0
maxAge = 100        # Myr
kappa = 1e-6
SecYear = 3600 * 24 * 365

# Function to calculate the perpendicular distance from a point to the ridge segment
function distance_to_line(x, y, x1, y1, x2, y2)
    num = abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    den = sqrt((y2 - y1)^2 + (x2 - x1)^2)
    return num / den
end

# Function to determine the side of a point respect to the delimiter
# Adapted to take into account the direction of the segments
function side_of_line(x, y, x1, y1, x2, y2, direction)
    side = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
    if direction == :left
        return side > 0
    else
        return side < 0
    end
end

# Function to determine in which region a point is, adjusted by address
function determine_region(px, py, delimiters, segments)
    for i in 1:length(delimiters)
        x1, y1 = delimiters[i][1]
        x2, y2 = delimiters[i][2]
        
        # Determine the direction of the segments
        if x2 < x1
            direction = :left  # Shift left
        else
            direction = :right  # Shift to the right
        end

        # Check the side of the line considering the direction
        if side_of_line(px, py, x1, y1, x2, y2, direction)
            return i  # Region corresponding to segment i
        end
    end
    return length(segments)  # Last region
end

# Function to calculate the ridge thermal structure

function compute_ridge_thermal_structure!(Temp, x, y, z, Phases; 
                                          segments, delimiters, Tsurface, Tmantle, 
                                          Adiabat, SpreadingVel, AgeRidge, maxAge, 
                                          kappa, SecYear)
    MantleAdiabaticT = Tmantle .+ Adiabat * abs.(z)

    for ix in 1:length(x), iy in 1:length(y), iz in 1:length(z)
        px, py, pz = x[ix], y[iy], z[iz]

        # Determine region of point
        region = determine_region(px, py, delimiters, segments)

        # Select the corresponding segment
        x1, y1 = segments[region][1]
        x2, y2 = segments[region][2]

        # Calculate distance to segment
        Distance = distance_to_line(px, py, x1, y1, x2, y2)

        # Calculate thermal age
        ThermalAge = abs(Distance * 1e3 * 1e2) / SpreadingVel + AgeRidge * 1e6
        ThermalAge = min(ThermalAge, maxAge * 1e6) * SecYear

        # Calculate temperature
        Temp[ix, iy, iz] = (Tsurface - Tmantle) * erfc(abs(pz) * 1e3 / (2 * sqrt(kappa * ThermalAge))) + MantleAdiabaticT[iz]
    end
end

# Call the function to calculate the thermal structure
compute_ridge_thermal_structure!(Temp, x, y, z, Phases; 
                                 segments=segments, delimiters=delimiters, 
                                 Tsurface=Tsurface, Tmantle=Tmantle, Adiabat=Adiabat, 
                                 SpreadingVel=SpreadingVel, AgeRidge=AgeRidge, maxAge=maxAge, 
                                 kappa=kappa, SecYear=SecYear)

# Add and save results
Grid = addfield(Grid, (; Phases, Temp))
write_paraview(Grid, "Ridge_Thermal_Structure")
