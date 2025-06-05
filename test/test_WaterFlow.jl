using Test, GMT


# Download some topographic data
Topo = import_topo([6.5, 7.3, 50.2, 50.6], file = "@earth_relief_03s");

# Flow the water through the area:
Topo_water, sinks, pits, bnds = waterflows(Topo)

@test maximum(Topo_water.fields.area) ≈ 9.309204547276944e8
@test sum(Topo_water.fields.c) == 834501044
@test sum(Topo_water.fields.nin) == 459361
@test sum(Topo_water.fields.dir) == 2412566

# With rain in m3/s per cell
rainfall = ones(size(Topo.lon.val[:, :, 1])) * 1.0e-3 # 2D array with rainfall per cell area
Topo_water1, sinks, pits, bnds = waterflows(Topo, rainfall = rainfall)

@test maximum(Topo_water1.fields.area) ≈ 169.79800000000208
