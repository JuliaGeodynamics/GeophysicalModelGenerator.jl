module GLMakie_Visualisation
# This contains visualisation widgets which are optionally made available when GLMakie is loaded along with GMG

using Statistics
using GeophysicalModelGenerator: lonlatdepth_grid, GeoData, CartData, km, AbstractGeneralGrid
import GeophysicalModelGenerator: visualise, ustrip

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
using GLMakie

import GLMakie: heatmap!, heatmap

export visualise, heatmap, heatmap!

println("Loading GLMakie extensions for GMG")

function visualise(Data::AbstractGeneralGrid; Topography = nothing, Topo_range = nothing)


    axis_equal = false   # in case we use x/y/z data in km, this is useful

    if isa(Data, GeoData)
        x = Data.lon.val[:, 1, 1];    xlab = "lon"
        y = Data.lat.val[1, :, 1];    ylab = "lat"
        z = Data.depth.val[1, 1, :];  zlab = "depth [km]"
    elseif isa(Data, CartData)
        # Determine
        x = Data.x.val[:, 1, 1];    xlab = "X [km]"
        y = Data.y.val[1, :, 1];    ylab = "Y [km]"
        z = Data.z.val[1, 1, :];  zlab = "Z [km]"
        axis_equal = true
    else
        error("not yet implemented ")
    end

    length(z) > 1 || error("visualise needs a 3D data set with more than one point in the vertical direction")
    if isa(Data, CartData)
        X, Y, Z = Data.x.val, Data.y.val, Data.z.val
        all(X .== X[:, 1:1, 1:1]) && all(Y .== Y[1:1, :, 1:1]) && all(Z .== Z[1:1, 1:1, :]) ||
            error("visualise needs an orthogonal grid; use project_CartData to project this CartData onto one")
    end

    if !axis_equal
        x_vec = 1:length(x)
        y_vec = 1:length(y)
        z_vec = 1:length(z)
    else
        x_vec = range(x[1], x[end], length(x))
        y_vec = range(y[1], y[end], length(y))
        z_vec = range(z[1], z[end], length(z))
    end
    # determine width of axis
    dx = (x[end] - x[1]) / (length(x) - 1)
    dy = (y[end] - y[1]) / (length(y) - 1)
    dz = (z[end] - z[1]) / (length(z) - 1)

    if !isnothing(Topography)
        if isa(Topography, GeoData)
            x_topo = (Topography.lon.val .- x[1]) / dx .+ 1
            y_topo = (Topography.lat.val .- y[1]) / dy .+ 1
            z_topo = (Topography.depth.val .- z[1]) / dz .+ 1
        elseif isa(Topography, CartData)
            x_topo = Topography.x.val
            y_topo = Topography.y.val
            z_topo = Topography.z.val
        end

        if isnothing(Topo_range)
            # base
            topo_max = maximum(abs, ustrip.(Topography.fields.Topography))
            iszero(topo_max) && (topo_max = 1.0)
            Topo_range = (-topo_max, topo_max)
        end

    end

    # Names of the scalar fields (vector fields stored as Tuples cannot be sliced)
    data_names = filter(k -> !(Data.fields[k] isa Tuple), collect(keys(Data.fields)))
    isempty(data_names) && error("visualise needs at least one scalar field in the data set")
    data_selected = Observable(Symbol(data_names[1]))
    data_string = @lift String($data_selected)
    get_vol(f_name) = ustrip.(Data.fields[f_name])
    vol = lift(get_vol, data_selected)

    fig = Figure(size = (2000, 2000), fontsize = 20)
    ax = Axis3(
        fig[1, 1:2], xlabel = xlab, ylabel = ylab, zlabel = zlab,
        aspect = axis_equal ? :data : (1, 1, 1)
    )

    # Create sliders
    slider_specs = Any[
        (
            label = xlab, range = 1:length(x_vec),
            format = v -> string(round((v - 1) * dx + x[1], digits = 2)),
        ),
        (
            label = ylab, range = 1:length(y_vec),
            format = v -> string(round((v - 1) * dy + y[1], digits = 2)),
        ),
        (
            label = zlab, range = 1:length(z_vec),
            format = v -> string(round((v - 1) * dz + z[1], digits = 2)),
        ),
    ]
    isnothing(Topography) || push!(slider_specs, (label = "Transparency topo", range = 0:0.01:1))
    sgrid = SliderGrid(fig[2, 2], slider_specs...)

    # Create dropdown menus
    menu_dataset = Menu(fig, options = [String.(data_names)...], default = String(data_selected[]))
    menu_colormap = Menu(
        fig, options = ["roma", "romaO", "vik", "turku", "davos", "batlow", "tab10", "tab20", "bone"],
        default = "roma"
    )

    # Colorbar limits
    # Equal extrema are widened, since a degenerate color range cannot be resolved
    crange = @lift begin
        vmin, vmax = extrema($vol)
        vmin == vmax && (vmax = vmin + (iszero(vmin) ? 1 : abs(vmin)))
        (round(vmin, RoundDown; sigdigits = 4), round(vmax, RoundUp; sigdigits = 4))
    end
    cmin = @lift $crange[1]
    cmax = @lift $crange[2]

    cmin_box = Textbox(fig, stored_string = string(cmin[]), width = nothing)
    cmax_box = Textbox(fig, stored_string = string(cmax[]), width = nothing)

    iso_level = Observable([1.7])
    iso_alpha = Observable(0.5)
    iso_box = Textbox(fig, stored_string = "$(iso_level[][1])", validator = Float64, width = nothing)
    iso_toggle = Toggle(fig, active = false)
    iso_slide = Slider(fig, range = 0:0.01:1)
    set_close_to!(iso_slide, iso_alpha[])

    Label(fig[3, 1:4], " ", width = nothing)

    fig[2, 1] = vgrid!(
        hgrid!(Label(fig, "Dataset", width = nothing), menu_dataset),
        hgrid!(Label(fig, "Colormap", width = nothing), menu_colormap),
        hgrid!(Label(fig, "Color axis limits", width = nothing), hgrid!(cmin_box, Label(fig, "-", width = 20), cmax_box)),
        hgrid!(hgrid!(Label(fig, "Isovalue", width = nothing), iso_toggle), hgrid!(iso_box, Label(fig, "α: ", width = nothing), iso_slide));
        tellheight = false
    )


    lo = sgrid.layout
    nc = ncols(lo)

    # GeoData is plotted in index coordinates, with tick labels showing lon/lat/depth: the contour
    # isosurface needs a uniformly spaced grid, which index coordinates provide.
    plt = volumeslices!(ax, x_vec, y_vec, z_vec, vol, colorrange = (cmin[], cmax[]), colormap = :roma)
    iso = contour!(
        ax, extrema(x_vec), extrema(y_vec), extrema(z_vec), vol, levels = iso_level, alpha = iso_alpha,
        visible = iso_toggle.active[], colormap = plt.colormap, colorrange = plt.colorrange
    )

    topo_alpha = Observable(0.5)
    if !isnothing(Topography)

        # in case topography is supplied
        topo_surf = surface!(ax, x_topo[:, :, 1], y_topo[:, :, 1], z_topo[:, :, 1], colormap = (:oleron, topo_alpha[]), color = ustrip.(Topography.fields.Topography[:, :, 1]), colorrange = Topo_range, transparency = true)
        cb_surf = Colorbar(fig[1, 4], topo_surf, vertical = true, label = "Topography", height = Relative(0.6))

    end

    if !axis_equal
        ax.xticks = ([1, length(x)], string.(round.([x[1], x[end]], digits = 2)))
        ax.yticks = ([1, length(y)], string.(round.([y[1], y[end]], digits = 2)))
        ax.zticks = ([1, length(z)], string.(round.([z[1], z[end]], digits = 2)))
    end

    #
    cb = Colorbar(fig[1, 3], plt, vertical = true, label = data_string, height = Relative(0.6))

    # connect sliders to `volumeslices` update methods
    sl_yz, sl_xz, sl_xy = sgrid.sliders[1:3]

    on(sl_yz.value) do v
        plt.yz_index = v
    end
    on(sl_xz.value) do v
        plt.xz_index = v
    end
    on(sl_xy.value) do v
        plt.xy_index = v
    end

    set_close_to!(sl_yz, 0.5length(x_vec))
    set_close_to!(sl_xz, 0.5length(y_vec))
    set_close_to!(sl_xy, 0.5length(z_vec))
    isnothing(Topography) || set_close_to!(sgrid.sliders[4], 0.5)

    # change color limits; the textboxes reject an inverted range, which stalls the Colorbar tick computation
    cmin_box.validator = s -> something(tryparse(Float64, s), Inf) < plt.colorrange[][2]
    cmax_box.validator = s -> something(tryparse(Float64, s), -Inf) > plt.colorrange[][1]
    on(cmin_box.stored_string) do s
        plt.colorrange = (parse(Float64, s), plt.colorrange[][2])
    end
    on(cmax_box.stored_string) do s
        plt.colorrange = (plt.colorrange[][1], parse(Float64, s))
    end

    # Change data
    on(menu_dataset.selection) do s
        data_selected[] = Symbol(s)
        plt.colorrange = (cmin[], cmax[])
        for (box, c) in ((cmin_box, cmin[]), (cmax_box, cmax[]))
            box.stored_string = box.displayed_string = string(c)
        end
    end

    # Change colormap
    on(menu_colormap.selection) do s
        plt.colormap = s
    end

    # Create isosurface?
    on(iso_toggle.active) do s
        iso.visible = s
    end
    on(iso_box.stored_string) do s
        iso_level[] = [parse(Float64, s)]
    end
    on(iso_slide.value) do v
        iso_alpha[] = v
    end

    if !isnothing(Topography)
        on(sgrid.sliders[4].value) do v
            topo_alpha[] = v
            topo_surf.colormap = (:oleron, topo_alpha[])
        end
    end

    # add toggles to show/hide slices
    hmaps = [plt[Symbol(:heatmap_, s)][] for s in (:yz, :xz, :xy)]
    toggles = [Toggle(lo[i, nc + 1], active = true) for i in 1:(isnothing(Topography) ? 3 : 4)]

    shown = isnothing(Topography) ? hmaps : [hmaps; topo_surf]
    for (h, t) in zip(shown, toggles)
        on(a -> h.visible = a, t.active)
    end

    display(fig)

    return fig
end

"""
    heatmap(x::GeoData, args...; field=:Topography, kwargs...)
heatmap for a 2D GeoData object (surface)
"""
function heatmap(x::GeoData, args...; field = :Topography, kwargs...)
    @assert size(x.depth.val, 3) == 1
    @assert hasfield(typeof(x.fields), field)

    return heatmap(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(x.fields[field][:, :, 1]), args...; axis = (xlabel = "lon", ylabel = "lat"), kwargs...)

end

"""
    heatmap(x::GeoData, a::Array{_T,N}, args...; kwargs...)
in-place heatmap for a 2D GeoData object (surface) with an array `a`. 
"""
function heatmap(x::GeoData, a::Array{_T, N}, args...; kwargs...) where {_T, N}
    @assert size(x.depth.val, 3) == 1

    return if N == 3
        heatmap(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(a[:, :, 1]), args...; axis = (xlabel = "lon", ylabel = "lat"), kwargs...)
    elseif N == 2
        heatmap(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(a), args...; axis = (xlabel = "lon", ylabel = "lat"), kwargs...)
    else
        error("heatmap needs a 2D or 3D array, got a $(N)D array")
    end

end

"""
    heatmap(x::CartData, args...; field=:Topography, kwargs...)
heatmap for a 2D CartData object (surface)
"""
function heatmap(x::CartData, args...; field = :Topography, kwargs...)
    @assert size(x.z.val, 3) == 1
    @assert hasfield(typeof(x.fields), field)

    return heatmap(x.x.val[:, 1], x.y.val[1, :], ustrip.(x.fields[field][:, :, 1]), args...; axis = (xlabel = "x [km]", ylabel = "y [km]"), kwargs...)

end

"""
    heatmap(x::CartData, a::Array{_T,N}, args...; kwargs...)
in-place heatmap for a 2D CartData object (surface) with an array `a`. 
"""
function heatmap(x::CartData, a::Array{_T, N}, args...; kwargs...) where {_T, N}
    @assert size(x.z.val, 3) == 1

    return if N == 3
        heatmap(x.x.val[:, 1], x.y.val[1, :], ustrip.(a[:, :, 1]), args...; axis = (xlabel = "x [km]", ylabel = "y [km]"), kwargs...)
    elseif N == 2
        heatmap(x.x.val[:, 1], x.y.val[1, :], ustrip.(a), args...; axis = (xlabel = "x [km]", ylabel = "y [km]"), kwargs...)
    else
        error("heatmap needs a 2D or 3D array, got a $(N)D array")
    end

end


"""
    heatmap!(x::GeoData, args...; field=:Topography, kwargs...)
in-place heatmap for a 2D GeoData object (surface), 
"""
function heatmap!(x::GeoData, args...; field = :Topography, kwargs...)
    @assert size(x.depth.val, 3) == 1
    @assert hasfield(typeof(x.fields), field)

    return heatmap!(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(x.fields[field][:, :, 1]), args...; kwargs...)

end

"""
    heatmap!(x::GeoData, a::Array{_T,N}, args...; kwargs...)
in-place heatmap for a 2D GeoData object (surface) with an array `a`. 
"""
function heatmap!(x::GeoData, a::Array{_T, N}, args...; kwargs...) where {_T, N}
    @assert size(x.depth.val, 3) == 1

    return if N == 3
        heatmap!(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(a[:, :, 1]), args...; kwargs...)
    elseif N == 2
        heatmap!(x.lon.val[:, 1], x.lat.val[1, :], ustrip.(a), args...; kwargs...)
    else
        error("heatmap needs a 2D or 3D array, got a $(N)D array")
    end

end

"""
    heatmap!(x::CartData, args...; field=:Topography, colorbar=false, kwargs...)
in-place heatmap for a 2D CartData object (surface). With `colorbar=true`, a colorbar is added
to the right of the plot, at position `[1, 2]` of the current figure. Returns the heatmap plot.
"""
function heatmap!(x::CartData, args...; field = :Topography, colorbar = false, kwargs...)
    @assert size(x.z.val, 3) == 1
    @assert hasfield(typeof(x.fields), field)

    data = ustrip.(x.fields[field][:, :, 1])

    hm = heatmap!(x.x.val[:, 1], x.y.val[1, :], data, args...; kwargs...)
    colorbar && Colorbar(current_figure()[1, 2], hm, label = string(field))

    return hm
end

"""
    heatmap!(x::CartData, a::Array{_T,N}, args...; kwargs...)
in-place heatmap for a 2D CartData object (surface) with an array `a`. 
"""
function heatmap!(x::CartData, a::Array{_T, N}, args...; kwargs...) where {_T, N}
    @assert size(x.z.val, 3) == 1

    return if N == 3
        heatmap!(x.x.val[:, 1], x.y.val[1, :], ustrip.(a[:, :, 1]), args...; kwargs...)
    elseif N == 2
        heatmap!(x.x.val[:, 1], x.y.val[1, :], ustrip.(a[:, :]), args...; kwargs...)
    else
        error("heatmap needs a 2D or 3D array, got a $(N)D array")
    end

end

end
