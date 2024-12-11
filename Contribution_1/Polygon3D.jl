# This script creates 3D polygons from "x" and "y" coordinates, and allows you to define a depth range for the polygon.
# It assigns phase and temperature values to the grid within the polygon.

using GeophysicalModelGenerator     
using CSV 
using DataFrames  

# Grid parameters
nx, ny, nz = 512, 512, 128
x = range(-1000, 1000, length=nx)
y = range(-1000, 1000, length=ny)
z = range(-660, 0, length=nz)
Grid = CartData(xyz_grid(x, y, z))

# Phases and temperature
Phases = fill(2, nx, ny, nz)
Temp = fill(1350.0, nx, ny, nz)

# Set the desired depth range for the polygon
z_lower_limit = -80.0
z_upper_limit = 0.0

# Define phase and temperature values for the polygon
phase_poly = 5
temp_poly = 1000.0

# Read polygon using CSV and DataFrame
polygon_file = "Test_block.txt"
df = CSV.read(polygon_file, DataFrame, delim=' ')  # Read CSV file with a space as delimiter
xpoly, ypoly = df[:, 1], df[:, 2]  # Extract coordinates

# Function to verify if a point (px, py) is inside a polygon defined by the coordinates (poly_x, poly_y)
function point_in_polygon(px, py, poly_x, poly_y)
    n = length(poly_x)
    inside = false
    j = n
    for i in 1:n
        if ((poly_y[i] > py) != (poly_y[j] > py)) && 
           (px < (poly_x[j] - poly_x[i]) * (py - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])
            inside = !inside
        end
        j = i
    end
    return inside
end

# Assign phase and temperature for the polygon within a specific depth range
function assign_polygon_properties!(Phases, Temp, x, y, z, xpoly, ypoly, phase_value, temp_value, z_lower_limit, z_upper_limit)
    for ix in 1:length(x), iy in 1:length(y), iz in 1:length(z)
        px, py, pz = x[ix], y[iy], z[iz]

        # Check if the point is inside the polygon (in 2D projection) and within the z-range
        if point_in_polygon(px, py, xpoly, ypoly) && pz >= z_lower_limit && pz <= z_upper_limit
            Phases[ix, iy, iz] = phase_value
            Temp[ix, iy, iz] = temp_value
        end
    end
end

# Assign properties within the specified depth range
assign_polygon_properties!(Phases, Temp, x, y, z, xpoly, ypoly, phase_poly, temp_poly, z_lower_limit, z_upper_limit)

# Add and save results
Grid = addfield(Grid, (; Phases, Temp))
write_paraview(Grid, "Polygon3D")


