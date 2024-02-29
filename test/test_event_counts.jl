using Test

# Create a CartGrid
Grid_cart = CartData(XYZGrid(-20:20,-20:.1:20,-30:30))

Lon,Lat,Depth = LonLatDepthGrid(-20:20,-20:.1:20,-30:30);
Grid_geo = GeoData(Lon,Lat,Depth,(;Depth))

using Random
rng = MersenneTwister(123)

# Create some random point data
pt = rand(rng, 10_000,3) .- 0.5
pt[:,1] .*= 20
pt[:,2] .*= 20
pt[:,3] .*= 30

EQ_cart = CartData(pt[:,1],pt[:,2], pt[:,3], (z=pt[:,3],))
EQ_geo  = GeoData(pt[:,1],pt[:,2], pt[:,3], (z=pt[:,3],))

R = (sum(pt.^2,dims=2)).^(1/3)
ind = findall(R.< 5);

# test the basic routine
counts = PointData2NearestGrid(pt[:,1],pt[:,2],pt[:,3], NumValue(Grid.x),NumValue(Grid.y), NumValue(Grid.z); radius_factor=2)
@test extrema(counts) == (0, 86)

# Test if the grid is on a CartData grid
Grid_Count = PointData2NearestGrid(pt[:,1],pt[:,2],pt[:,3], Grid_cart; radius_factor=2)
@test extrema(Grid_Count.fields.Count) == (0, 86)

# Test in case the EQ data is also specified as CartData
Grid_Count = PointData2NearestGrid(EQ_cart, Grid_cart; radius_factor=2)
@test extrema(Grid_Count.fields.Count) == (0, 86)

# Test if the grid is on a GeoData grid
Grid_Count = PointData2NearestGrid(pt[:,1],pt[:,2],pt[:,3], Grid_geo; radius_factor=2)
@test extrema(Grid_Count.fields.Count) == (0, 86)

# Test in case the EQ data is also specified as GeoData
Grid_Count = PointData2NearestGrid(EQ_geo, Grid_geo; radius_factor=2)
@test extrema(Grid_Count.fields.Count) == (0, 86)
