# Parallel (multithreaded CPU or GPU) versions of the gravity kernels of voxel_grav and voxel_grav_adjoint
module KernelAbstractions_utils

using GeophysicalModelGenerator, KernelAbstractions
import GeophysicalModelGenerator: grav_forward_ka, grav_adjoint_ka, grav_points_forward_ka, grav_points_adjoint_ka, G_grav, mGal, box_field

# anomaly and horizontal gradients (if grad) at station (i, j): sum over the cells k within the horizontal distance sqrt(c2)
@kernel function forward_kernel!(dg, gx, gy, @Const(x), @Const(y), @Const(zo), @Const(xc), @Const(yc), @Const(zc), @Const(c), @Const(m1), @Const(m2), @Const(m3), c2, grad)
    i, j = @index(Global, NTuple)
    s = sx = sy = zero(eltype(dg))
    for k in 1:length(c)
        dx = x[i] - xc[k]
        dy = y[j] - yc[k]
        h = zo[i, j] - zc[k]
        gz, ge, gn = box_field(ifelse(dx * dx + dy * dy <= c2, c[k], zero(c2)), dx, dy, h, m1[k], m2[k], m3[k])
        s += gz
        sx += ge
        sy += gn
    end
    dg[i, j] = s
    if grad
        gx[i, j] = sx
        gy[i, j] = sy
    end
end

# sensitivity of cell k: sum over the stations (i, j)
@kernel function adjoint_kernel!(S, @Const(x), @Const(y), @Const(zo), @Const(R), @Const(xc), @Const(yc), @Const(zc), @Const(m1), @Const(m2), @Const(m3), c2)
    k = @index(Global)
    s = zero(eltype(S))
    for j in 1:length(y), i in 1:length(x)
        dx = x[i] - xc[k]
        dy = y[j] - yc[k]
        h = zo[i, j] - zc[k]
        s += box_field(ifelse(dx * dx + dy * dy <= c2, R[i, j], zero(c2)), dx, dy, h, m1[k], m2[k], m3[k])[1]
    end
    S[k] = s
end

function device(backend, T, v)
    a = KernelAbstractions.allocate(backend, T, size(v))
    copyto!(a, T.(v))
    return a
end

# stations, and centres and second moments of the cells I, relative to the first grid point (keeps Float32 accurate)
function device_geometry(backend, T, geo, I)
    (; x, y, xc, yc, zc, zo, mx, my, mz) = geo
    x0, y0 = x[1], y[1]
    d(v) = device(backend, T, v)
    return d(x .- x0), d(y .- y0), d(zo), d([xc[i[1]] - x0 for i in I]), d([yc[i[2]] - y0 for i in I]), d([zc[i[3]] for i in I]),
        d([mx[i[1]] for i in I]), d([my[i[2]] for i in I]), d([mz[i[3]] for i in I])
end

# squared cutoff in the number type of the device (Inf stays Inf)
cutoff2(T, cutoff) = T(min(cutoff, floatmax(Float64))^2)

# anomaly and gradients (or nothing) as Float64 arrays on the host; the gradient arrays are dummies if !grad
function run_forward(backend, T, kernel!, n, args, grad)
    dg = KernelAbstractions.zeros(backend, T, n...)
    gx = KernelAbstractions.zeros(backend, T, (grad ? n : map(one, n))...)
    gy = KernelAbstractions.zeros(backend, T, (grad ? n : map(one, n))...)
    kernel!(backend)(dg, gx, gy, args..., grad; ndrange = n)
    KernelAbstractions.synchronize(backend)
    host(a) = Float64.(Array(a))
    return host(dg), (grad ? host(gx) : nothing), (grad ? host(gy) : nothing)
end

function grav_forward_ka(backend::Backend, T, geo, DRHO, rhoTol, cutoff, grad)
    I = findall(a -> abs(a) > rhoTol, DRHO)
    x, y, zo, xc, yc, zc, m1, m2, m3 = device_geometry(backend, T, geo, I)
    c = device(backend, T, G_grav * mGal .* geo.vol[I] .* DRHO[I])
    return run_forward(backend, T, forward_kernel!, (length(geo.x), length(geo.y)), (x, y, zo, xc, yc, zc, c, m1, m2, m3, cutoff2(T, cutoff)), grad)
end

function grav_adjoint_ka(backend::Backend, T, geo, R, cutoff)
    I = CartesianIndices(geo.vol)
    x, y, zo, xc, yc, zc, m1, m2, m3 = device_geometry(backend, T, geo, I)
    S = KernelAbstractions.zeros(backend, T, length(I))
    adjoint_kernel!(backend)(S, x, y, zo, device(backend, T, R), xc, yc, zc, m1, m2, m3, cutoff2(T, cutoff); ndrange = length(I))
    KernelAbstractions.synchronize(backend)
    return G_grav * mGal .* geo.vol .* reshape(Float64.(Array(S)), size(geo.vol))
end

# point masses with stations on the ellipsoid (GeoData): component along the station normal
@kernel function points_forward_kernel!(
        dg, gx, gy, @Const(sx), @Const(sy), @Const(sz), @Const(nx), @Const(ny), @Const(nz), @Const(ex), @Const(ey), @Const(qx), @Const(qy), @Const(qz),
        @Const(px), @Const(py), @Const(pz), @Const(c), @Const(m1), @Const(m2), @Const(m3), c2, grad
    )
    i = @index(Global)
    s = se = sn = zero(eltype(dg))
    for k in 1:length(c)
        dx, dy, dz = sx[i] - px[k], sy[i] - py[k], sz[i] - pz[k]
        h = dx * nx[i] + dy * ny[i] + dz * nz[i]
        t1 = dx * ex[i] + dy * ey[i]
        t2 = dx * qx[i] + dy * qy[i] + dz * qz[i]
        gz, ge, gn = box_field(ifelse(t1 * t1 + t2 * t2 <= c2, c[k], zero(c2)), t1, t2, h, m1[k], m2[k], m3[k])
        s += gz
        se += ge
        sn += gn
    end
    dg[i] = s
    if grad
        gx[i] = se
        gy[i] = sn
    end
end

@kernel function points_adjoint_kernel!(
        S, @Const(sx), @Const(sy), @Const(sz), @Const(nx), @Const(ny), @Const(nz), @Const(ex), @Const(ey), @Const(qx), @Const(qy), @Const(qz),
        @Const(R), @Const(px), @Const(py), @Const(pz), @Const(m1), @Const(m2), @Const(m3), c2
    )
    k = @index(Global)
    s = zero(eltype(S))
    for i in 1:length(sx)
        dx, dy, dz = sx[i] - px[k], sy[i] - py[k], sz[i] - pz[k]
        h = dx * nx[i] + dy * ny[i] + dz * nz[i]
        t1 = dx * ex[i] + dy * ey[i]
        t2 = dx * qx[i] + dy * qy[i] + dz * qz[i]
        s += box_field(ifelse(t1 * t1 + t2 * t2 <= c2, R[i], zero(c2)), t1, t2, h, m1[k], m2[k], m3[k])[1]
    end
    S[k] = s
end

stations(backend, T, geo) = map(v -> device(backend, T, v), (geo.sx, geo.sy, geo.sz, geo.nx, geo.ny, geo.nz))
directions(backend, T, geo) = map(v -> device(backend, T, v), (geo.ex, geo.ey, geo.qx, geo.qy, geo.qz))

function grav_points_forward_ka(backend::Backend, T, geo, DRHO, rhoTol, cutoff, grad)
    I = findall(a -> abs(a) > rhoTol, DRHO)
    d(v) = device(backend, T, v)
    c = d(G_grav * mGal .* geo.vol[I] .* DRHO[I])
    args = (stations(backend, T, geo)..., directions(backend, T, geo)..., d(geo.px[I]), d(geo.py[I]), d(geo.pz[I]), c, d(geo.m1[I]), d(geo.m2[I]), d(geo.m3[I]), cutoff2(T, cutoff))
    return run_forward(backend, T, points_forward_kernel!, (length(geo.sx),), args, grad)
end

function grav_points_adjoint_ka(backend::Backend, T, geo, R, cutoff)
    d(v) = device(backend, T, v)
    S = KernelAbstractions.zeros(backend, T, length(geo.vol))
    points_adjoint_kernel!(backend)(
        S, stations(backend, T, geo)..., directions(backend, T, geo)..., d(R), d(vec(geo.px)), d(vec(geo.py)), d(vec(geo.pz)),
        d(vec(geo.m1)), d(vec(geo.m2)), d(vec(geo.m3)), cutoff2(T, cutoff); ndrange = length(S)
    )
    KernelAbstractions.synchronize(backend)
    return G_grav * mGal .* vec(geo.vol) .* Float64.(Array(S))
end

end
