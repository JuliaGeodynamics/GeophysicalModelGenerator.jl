using Test
# test various surface routines

# Create surfaces
cartdata1 = CartData(xyz_grid(1:4, 1:5, 0))
cartdata2 = CartData(xyz_grid(1:4, 1:5, 2))
cartdata3 = CartData(xyz_grid(1:4, 1:5, 2:5))
cartdata2 = addfield(cartdata2, "Z2", cartdata2.x.val)

@test is_surface(cartdata1)
@test is_surface(cartdata2)
@test is_surface(cartdata3) == false

geodata1 = GeoData(lonlatdepth_grid(1:4, 1:5, 0))
geodata2 = GeoData(lonlatdepth_grid(1:4, 1:5, 2))
geodata3 = GeoData(lonlatdepth_grid(1:4, 1:5, 2:5))

@test is_surface(geodata1)
@test is_surface(geodata2)
@test is_surface(geodata3) == false

# Test add & subtraction of surfaces
cartdata4 = cartdata1 + cartdata2
@test length(cartdata4.fields) == 2
@test cartdata4.z.val[2] == 2.0

cartdata5 = cartdata1 - cartdata2
@test length(cartdata5.fields) == 2
@test cartdata5.z.val[2] == -2.0

geodata4 = geodata1 + geodata2
@test length(geodata4.fields) == 1
@test geodata4.depth.val[2] == 2.0

geodata5 = geodata1 - geodata2
@test length(geodata5.fields) == 1
@test geodata5.depth.val[2] == -2.0

# Test removing NaN;
Z = NumValue(cartdata5.z)
Z[2, 2] = NaN;
remove_NaN_surface!(Z, NumValue(cartdata5.x), NumValue(cartdata5.y))
@test any(isnan.(Z)) == false

# Test draping values on topography
X, Y, Z = xyz_grid(1:0.14:4, 1:0.02:5, 0);
v = X .^ 2 .+ Y .^ 2;
values1 = CartData(X, Y, Z, (; v))
values2 = CartData(X, Y, Z, (; colors = (v, v, v)))

cart_drape1 = drape_on_topo(cartdata2, values1)
@test  sum(cart_drape1.fields.v) ≈ 366.02799999999996

cart_drape2 = drape_on_topo(cartdata2, values2)
@test  cart_drape2.fields.colors[1][10] ≈ 12.9204

values1 = GeoData(X, Y, Z, (; v))
values2 = GeoData(X, Y, Z, (; colors = (v, v, v)))

geo_drape1 = drape_on_topo(geodata2, values1)
@test  sum(geo_drape1.fields.v) ≈ 366.02799999999996

geo_drape2 = drape_on_topo(geodata2, values2)
@test  geo_drape2.fields.colors[1][10] ≈ 12.9204

# test fit_surface_to_points
cartdata2b = fit_surface_to_points(cartdata2, X[:], Y[:], v[:])
@test sum(NumValue(cartdata2b.z)) ≈ 366.02799999999996


#-------------
# test above_surface with the Grid object
Grid = create_CartGrid(size = (10, 20, 30), x = (0.0, 10), y = (0.0, 10), z = (-10.0, 2.0))
@test Grid.Δ[2] ≈ 0.5263157894736842

Temp = ones(Float64, Grid.N...) * 1350;
Phases = zeros(Int32, Grid.N...);

Topo_cart = CartData(xyz_grid(-1:0.2:20, -12:0.2:13, 0));
ind = above_surface(Grid, Topo_cart);
@test sum(ind[1, 1, :]) == 5

ind = below_surface(Grid, Topo_cart);
@test sum(ind[1, 1, :]) == 25


#-------------
# test above_surface with the Q1Data object
q1data = Q1Data(xyz_grid(1:4, 1:5, -5:5))
ind = above_surface(q1data, cartdata2);
@test sum(ind) == 60

ind = below_surface(q1data, cartdata2);
@test sum(ind) == 140

#-------------
# Add & subtract ParaviewData surfaces
pvdata1 = ParaviewData(xyz_grid(1:4, 1:5, 0)..., (a = ones(4, 5, 1),))
pvdata2 = ParaviewData(xyz_grid(1:4, 1:5, 2)..., (b = ones(4, 5, 1),))
pvdata3 = pvdata1 + pvdata2
@test keys(pvdata3.fields) == (:a, :b)
@test all(pvdata3.z.val .== 2.0)
pvdata4 = pvdata1 - pvdata2
@test all(pvdata4.z.val .== -2.0)

# above_surface & below_surface with ParaviewData
pvvol = ParaviewData(xyz_grid(1:4, 1:5, -5:5)..., (z = ones(4, 5, 11),))
ind = above_surface(pvvol, pvdata2)
@test size(ind) == (4, 5, 11)
@test sum(ind) == 60
@test all(ind[:, :, end])
ind = below_surface(pvvol, pvdata2)
@test sum(ind) == 140
@test all(ind[:, :, 1])

# interpolate_data_surface with GeoData: a linear field is interpolated exactly on a dipping surface
Lon, Lat, Depth = lonlatdepth_grid(1:4, 1:5, -5:5)
geovol = GeoData(Lon, Lat, Depth, (v = 2 .* ustrip.(Depth),))
Lon, Lat, Depth = lonlatdepth_grid(1.5:1:3.5, 1.5:1:4.5, 0)
Depth = -0.2 .* Lon .- 0.5 .* Lat
geosurf = GeoData(Lon, Lat, Depth, (a = zeros(size(Lon)),))
geosurf_interp = interpolate_data_surface(geovol, geosurf)
@test size(geosurf_interp.fields.v) == (3, 4, 1)
@test geosurf_interp.fields.v ≈ 2 .* Depth

# Add & subtract UTMData surfaces
utmdata1 = UTMData(xyz_grid(1:4, 1:5, 0)..., 33, true, (a = ones(4, 5, 1),))
utmdata2 = UTMData(xyz_grid(1:4, 1:5, 2)..., 33, true, (b = ones(4, 5, 1),))
utmdata3 = utmdata1 + utmdata2
@test keys(utmdata3.fields) == (:a, :b)
@test all(utmdata3.depth.val .== 2.0)
@test all(utmdata3.zone .== 33)
@test all(utmdata3.northern)
utmdata4 = utmdata1 - utmdata2
@test all(utmdata4.depth.val .== -2.0)

# fit_surface_to_points with GeoData: every surface point takes the depth of the closest point
Lon, Lat, Depth = lonlatdepth_grid(1:4, 1:5, 0)
geosurf0 = GeoData(Lon, Lat, Depth, (a = zeros(size(Lon)),))
geosurf_fit = fit_surface_to_points(geosurf0, Lon[:], Lat[:], -Lon[:])
@test geosurf_fit.depth.val ≈ -Lon
@test all(geosurf0.depth.val .== 0)        # input is not modified
geosurf_fit = fit_surface_to_points(geosurf0, [1.2, 3.9], [1.1, 4.8], [-3.0, -7.0])
@test geosurf_fit.depth.val[1, 1] == -3.0
@test geosurf_fit.depth.val[4, 5] == -7.0
@test sort(unique(geosurf_fit.depth.val)) == [-7.0, -3.0]

# interpolate_data_surface with CartData and ParaviewData
X, Y, Z = xyz_grid(1:4, 1:5, -5:5)
Xs, Ys, _ = xyz_grid(1.5:1:3.5, 1.5:1:4.5, 0)
Zs = -0.2 .* Xs .- 0.5 .* Ys
cartsurf_interp = interpolate_data_surface(CartData(X, Y, Z, (v = 2 .* Z,)), CartData(Xs, Ys, Zs, (a = zeros(size(Xs)),)))
@test cartsurf_interp isa CartData
@test cartsurf_interp.fields.v ≈ 2 .* Zs
@test ustrip.(cartsurf_interp.z.val) ≈ Zs
pvsurf_interp = interpolate_data_surface(ParaviewData(X, Y, Z, (v = 2 .* Z,)), ParaviewData(Xs, Ys, Zs, (a = zeros(size(Xs)),)))
@test pvsurf_interp isa ParaviewData
@test pvsurf_interp.fields.v ≈ 2 .* Zs
