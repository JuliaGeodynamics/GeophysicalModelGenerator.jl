# tests transformations from GeoData <=> Cartesian
using Test
using GeophysicalModelGenerator

# Create 3D volume with some fake data
Lon, Lat, Depth = lonlatdepth_grid(5:25, 20:50, (-1300:100:0)km);
Data_set3D = GeoData(Lon, Lat, Depth, (Depthdata = Depth * 2 + Lon * km, LonData = Lon))

proj = ProjectionPoint(Lon = 20, Lat = 35)

# Convert this 3D dataset to a Cartesian dataset (the grid will not be orthogonal)
Data_set3D_Cart = convert2CartData(Data_set3D, proj)
@test sum(abs.(Value(Data_set3D_Cart.x))) ≈ 5.293469089428514e6km

# Create Cartesian grid
X, Y, Z = xyz_grid(-400:100:400, -500:200:500, (-1300:100:0)km);
Data_Cart = CartData(X, Y, Z, (Z = Z,))

# Project values of Data_set3D to the cartesian data
Data_Cart = project_CartData(Data_Cart, Data_set3D, proj)
@test sum(Data_Cart.fields.Depthdata) ≈ -967680.9136292854km
#@test sum(Data_Cart.fields.Depthdata) ≈ -1.416834287168597e6km
@test sum(Data_Cart.fields.LonData) ≈ 15119.086370714615


# Next, 3D surface (like topography)
Lon, Lat, Depth = lonlatdepth_grid(5:25, 20:50, 0);
Depth = cos.(Lon / 5) .* sin.(Lat) * 10;
Data_surf = GeoData(Lon, Lat, Depth, (Z = Depth,));
Data_surf_Cart = convert2CartData(Data_surf, proj);

# Cartesian surface
X, Y, Z = xyz_grid(-500:10:500, -900:20:900, 0);
Data_Cart = CartData(X, Y, Z, (Z = Z,))

Data_Cart = project_CartData(Data_Cart, Data_surf, proj)
@test sum(Value(Data_Cart.z)) ≈ 1858.2487019158766km
@test sum(Data_Cart.fields.Z) ≈ 1858.2487019158766


# Cartesian surface when UTM data is used
WE, SN, depth = xyz_grid(420000:1000:430000, 4510000:1000:4520000, 0);

Data_surfUTM = UTMData(WE, SN, depth, 33, true, (Depth = WE,));
Data_Cart = CartData(X, Y, Z, (Z = Z,))
Data_Cart = project_CartData(Data_Cart, Data_surfUTM, proj)

@test sum(Value(Data_Cart.z)) ≈ 0.0km
@test sum(Data_Cart.fields.Depth) ≈ 3.9046959539921126e9

# GeoData given in 0-360° longitudes, projected around a point with negative longitude
Lon, Lat, Depth = lonlatdepth_grid(340:355, 30:40, 0);
Data_west = GeoData(Lon, Lat, Depth, (LonData = Lon,));
proj_west = ProjectionPoint(Lon = -12, Lat = 35)
Data_Cart = CartData(xyz_grid(-100:50:100, -100:50:100, 0)..., (Z = zeros(5, 5, 1),))
Data_Cart = project_CartData(Data_Cart, Data_west, proj_west)
@test size(Data_Cart.fields.LonData) == (5, 5, 1)
@test Data_Cart.fields.LonData[3, 3] ≈ 348.0
@test all(346 .< Data_Cart.fields.LonData .< 350)

# point_in_tetrahedron
a, b, c, d = [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
@test GeophysicalModelGenerator.point_in_tetrahedron([0.1, 0.1, 0.1], a, b, c, d)
@test GeophysicalModelGenerator.point_in_tetrahedron([0.0, 0.0, 0.0], a, b, c, d)        # vertex
@test !GeophysicalModelGenerator.point_in_tetrahedron([0.5, 0.5, 0.5], a, b, c, d)       # inside bounding box, outside tetrahedron
@test !GeophysicalModelGenerator.point_in_tetrahedron([2.0, 0.1, 0.1], a, b, c, d)       # outside bounding box

# project_FEData_CartData: two tetrahedra in opposite corners of the unit cube
vertices = [
    0.0 1.0 0.0 0.0 1.0 0.0 1.0 1.0;
    0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0;
    0.0 0.0 0.0 1.0 1.0 1.0 1.0 0.0
]
connectivity = [1 5; 2 6; 3 7; 4 8]
data_fe = FEData(vertices, connectivity, NamedTuple(), (regions = [1, 2],))
data_cart = CartData(xyz_grid(0.1:0.4:0.9, 0.1:0.4:0.9, 0.1:0.4:0.9))
data_cart = project_FEData_CartData(data_cart, data_fe)
regions = data_cart.fields.regions
@test size(regions) == (3, 3, 3)
@test regions[1, 1, 1] == 1
@test regions[3, 3, 3] == 2
@test regions[2, 2, 2] == 0
@test count(==(1), regions) == 4
@test count(==(2), regions) == 4

# 3D UTMData projected onto a Cartesian volume: linear fields are interpolated exactly
EW, NS, Depth_utm = xyz_grid(range(proj.EW - 20.0e3, proj.EW + 20.0e3, 5), range(proj.NS - 20.0e3, proj.NS + 20.0e3, 5), range(-10.0e3, 0, 3))
Data_UTM3D = UTMData(EW, NS, Depth_utm, proj.zone, proj.isnorth, (EW2 = EW, D2 = Depth_utm))
Data_Cart = CartData(xyz_grid(-10:5:10, -10:5:10, -8:4:0)..., (Z = zeros(5, 5, 3),))
Data_Cart = project_CartData(Data_Cart, Data_UTM3D, proj)
@test size(Data_Cart.fields.EW2) == (5, 5, 3)
@test Data_Cart.fields.EW2 ≈ proj.EW .+ 1000 .* Data_Cart.x.val
@test Data_Cart.fields.D2 ≈ 1000 .* ustrip.(Data_Cart.z.val)

# CartData projected onto another CartData grid (volume and surface)
X, Y, Z = xyz_grid(0:4, 0:3, -2:0)
Data_src = CartData(X, Y, Z, (T = X .+ 2 .* Y .+ 3 .* Z,))
X, Y, Z = xyz_grid(0.5:1:3.5, 0.5:1:2.5, -1.5:1:-0.5)
Data_Cart = project_CartData(CartData(X, Y, Z, (Z = Z,)), Data_src)
@test size(Data_Cart.fields.T) == (4, 3, 2)
@test Data_Cart.fields.T ≈ X .+ 2 .* Y .+ 3 .* Z
@test keys(Data_Cart.fields) == (:T,)

X, Y, _ = xyz_grid(0:4, 0:3, 0)
Data_src = CartData(X, Y, X .+ Y, (T = 2 .* X,))
X, Y, Z = xyz_grid(0.5:1:3.5, 0.5:1:2.5, 0)
Data_Cart = project_CartData(CartData(X, Y, Z, (Z = Z,)), Data_src)
@test ustrip.(Data_Cart.z.val) ≈ X .+ Y
@test Data_Cart.fields.T ≈ 2 .* X
