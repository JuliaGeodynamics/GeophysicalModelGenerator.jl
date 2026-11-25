# Normal gravity and Bouguer / terrain corrections, to turn observed gravity into anomalies that can be compared with
# voxel_grav

export normal_gravity, bouguer_correction, topographic_effect

# GRS80: semi-major axis a [m], flattening f, normal gravity at equator and pole [m/s²], m = ω²a²b/(GM)
const GRS80 = (a = 6378137.0, f = 1 / 298.257222101, γe = 9.7803267715, γp = 9.8321863685, m = 0.00344978600308)

"""
    γ = normal_gravity(lat, h = 0; ellipsoid = GRS80)

Normal gravity [mGal] of the reference ellipsoid (GRS80 by default) at geodetic latitude `lat` [°] and height `h` [km] above the ellipsoid: Somigliana's formula with the second-order height correction (Hofmann-Wellenhof & Moritz, 2006, eq. 2-215).
Subtracting it from observed (absolute) gravity gives the gravity disturbance (with ellipsoidal heights) or, approximately, the free-air anomaly (with heights above sea level). `ellipsoid` is a NamedTuple with fields `a, f, γe, γp, m`.
"""
function normal_gravity(lat, h = 0; ellipsoid = GRS80)
    (; a, f, γe, γp, m) = ellipsoid
    b = a * (1 - f)
    e2 = f * (2 - f)
    s2 = sind(lat)^2
    γ0 = γe * (1 + (b * γp / (a * γe) - 1) * s2) / sqrt(1 - e2 * s2)
    H = 1.0e3 * h
    return 1.0e5 * γ0 * (1 - 2 / a * (1 + f + m - 2f * s2) * H + 3 * H^2 / a^2)
end

"""
    bc = bouguer_correction(topo; density = 2670, density_water = 1040)

Simple Bouguer correction [mGal]: the attraction `2πGρh` of an infinite slab with the thickness of the topography `topo` [km] (a number or an array). For negative `topo` (oceans), the water is replaced by rock, i.e. the correction is `2πG(density_water - density)|topo|`. Subtract it from the free-air anomaly to get the simple Bouguer anomaly. [`topographic_effect`](@ref) includes the terrain.
"""
function bouguer_correction(topo; density = 2670.0, density_water = 1040.0)
    slab(h) = 2π * G_grav * (h >= 0 ? density * h : (density - density_water) * h) * 1.0e3 * mGal
    return slab.(topo)
end

"""
    te = topographic_effect(Topo; density = 2670, density_water = 1040, dz = 0.1, cutoff = 166.7, stations = nothing, kwargs...)

Complete Bouguer correction [mGal]: the attraction of the topography between sea level and the surface `Topo` with the density `density` (and of the oceans, as `density_water - density` between the sea floor and sea level), computed with [`voxel_grav`](@ref) at stations on the topography (or at sea level for the oceans), or at `stations`. Unlike [`bouguer_correction`](@ref), this includes the terrain correction. Subtract it from the free-air anomaly (or gravity disturbance) to get the complete Bouguer anomaly.

`Topo` is a `GeoData` surface (computed on the ellipsoid) or a `CartData` surface (flat Earth), with the topography in km as depth/z, e.g. from `import_topo`. The masses are discretized in layers of thickness `dz` [km]. Cells beyond the horizontal distance `cutoff` [km] are ignored; 166.7 km is the standard outer radius of terrain corrections (Hayford-Bowie zone O). Further keyword arguments (e.g. `backend`, `float_type`, `gradients`) are passed on to `voxel_grav`. Returns the result of `voxel_grav`, i.e. `(te, dtedx, dtedy)`.
"""
function topographic_effect(Topo::Union{GeoData, CartData}; density = 2670.0, density_water = 1040.0, dz = 0.1, cutoff = 166.7, stations = nothing, kwargs...)
    a, b, topo = topo_surface(Topo)
    zmin, zmax = min(0.0, minimum(topo)), max(0.0, maximum(topo))
    n = max(2, ceil(Int, (zmax - zmin) / dz) + 1)
    A, B, Z = lonlatdepth_grid(a, b, range(zmin, zmax, n))
    Grid = Topo isa GeoData ? GeoData(A, B, Z, (z = Z,)) : CartData(A, B, Z, (z = Z,))
    t = reshape(topo, size(topo, 1), size(topo, 2), 1)
    RHO = @. ifelse((Z >= 0) & (Z <= t), density, ifelse((Z < 0) & (Z >= t), density_water - density, 0.0))
    if isnothing(stations)
        return voxel_grav(Grid, RHO; refMod = zeros(n), zRef = max.(topo, 0.0), cutoff, kwargs...)
    end
    return voxel_grav(Grid, RHO; refMod = zeros(n), stations, cutoff, kwargs...)
end

function topo_surface(Topo::GeoData)
    size(Topo.lon.val, 3) == 1 || error("Topo must be a surface (size(Topo.lon, 3) == 1).")
    return Topo.lon.val[:, 1, 1], Topo.lat.val[1, :, 1], Topo.depth.val[:, :, 1]
end
function topo_surface(Topo::CartData)
    size(Topo.x.val, 3) == 1 || error("Topo must be a surface (size(Topo.x, 3) == 1).")
    return Topo.x.val[:, 1, 1], Topo.y.val[1, :, 1], Topo.z.val[:, :, 1]
end
