using Test

# Create a CartGrid
Grid_cart = CartData(xyz_grid(-20:20, -20:0.1:20, -30:30))

# create 2D GeoData struct
Lon, Lat, Depth = lonlatdepth_grid(10:20, 30:40, 0);
CM = zeros(size(Depth)); CM[1:5, 1:5] .= 1.0
Data_set2D = GeoData(Lon, Lat, Depth, (Count = CM,))

using StableRNGs

rng = StableRNG(123)

# Create some random point data
pt = rand(rng, 10_000, 3) .- 0.5
pt[:, 1] .*= 20
pt[:, 2] .*= 20
pt[:, 3] .*= 30

EQ_cart = CartData(pt[:, 1], pt[:, 2], pt[:, 3], (z = pt[:, 3],))

# test the basic routine
counts = point_to_nearest_grid(pt[:, 1], pt[:, 2], pt[:, 3], NumValue(Grid_cart.x), NumValue(Grid_cart.y), NumValue(Grid_cart.z); radius_factor = 2)
@test extrema(counts) == (0, 114)

# the search radius is a length: the same points and grid in m instead of km give the same counts
counts_m = point_to_nearest_grid(1000 .* pt[:, 1], 1000 .* pt[:, 2], 1000 .* pt[:, 3], 1000 .* NumValue(Grid_cart.x), 1000 .* NumValue(Grid_cart.y), 1000 .* NumValue(Grid_cart.z); radius_factor = 2)
@test counts_m == counts

# Test if the grid is on a CartData grid
Grid_Count = point_to_nearest_grid(pt[:, 1], pt[:, 2], pt[:, 3], Grid_cart; radius_factor = 2)
@test extrema(Grid_Count.fields.Count) == (0, 114)

# Test in case the EQ data is also specified as CartData
Grid_Count = point_to_nearest_grid(EQ_cart, Grid_cart; radius_factor = 2)
@test extrema(Grid_Count.fields.Count) == (0, 114)

# GeoData: the Cartesian grid and points, mapped to lon/lat around a projection point, give nearly the same counts,
# as distances are measured in km; Earth curvature causes the small differences
proj = ProjectionPoint(Lon = 10, Lat = 45)
Grid_geo = convert(GeoData, convert2UTMzone(Grid_cart, proj))
EQ_geo = convert(GeoData, convert2UTMzone(EQ_cart, proj))
Count_cart = point_to_nearest_grid(EQ_cart, Grid_cart; radius_factor = 2).fields.Count

# Test if the grid is on a GeoData grid
Grid_Count = point_to_nearest_grid(NumValue(EQ_geo.lon), NumValue(EQ_geo.lat), NumValue(EQ_geo.depth), Grid_geo; radius_factor = 2)
@test isapprox(sum(Grid_Count.fields.Count), sum(Count_cart); rtol = 0.01)
@test maximum(abs, Grid_Count.fields.Count .- Count_cart) <= 7

# Test in case the EQ data is also specified as GeoData
@test point_to_nearest_grid(EQ_geo, Grid_geo; radius_factor = 2).fields.Count == Grid_Count.fields.Count

# Test countmap
Data_countMap = countmap(Data_set2D, "Count", 5, 5)
@test Data_countMap.fields.countmap[1, 1] == 1.0
@test Data_countMap.fields.countmap[2, 2] == 0.4444444444444444
@test Data_countMap.fields.countmap[4, 4] == 0.0
