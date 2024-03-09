module GLMakie_Visualisation
# This contains visualisation widgets which are optionally made available when GLMakie is loaded along with GMG

using Statistics
using GeophysicalModelGenerator: LonLatDepthGrid, GeoData, CartData, km, AbstractGeneralGrid
import GeophysicalModelGenerator: Visualise

# We do not check `isdefined(Base, :get_extension)` as recommended since
# Julia v1.9.0 does not load package extensions when their dependency is
# loaded from the main environment.
if VERSION >= v"1.9.1"
    using GLMakie
else
    using ..GLMakie
end

export Visualise

println("Loading GLMakie extensions for GMG")

"""
    Visualise(DataSet; Topography=Topo_Data, Topo_range=nothing)

This starts an interactive widget that allows you to explore a 3D data set `DataSet` in an interactive manner.
All fields in the dataset can be explored, and if the optional parameter `Topography` is provided, the topography will be drawn on top.

Note that this requires orthogonal grids, so it will work with a `GeoData` set, or with an orthogonal `CartData` set.
Note that you may have to use `ProjectCartData` to project it to orthogonal cartesian coordinates.
"""
function Visualise(Data::AbstractGeneralGrid; Topography=nothing, Topo_range=nothing)


    axis_equal = false;  # in case we use x/y/z data in km, this is useful 

    if isa(Data,GeoData)
        x = Data.lon.val[:,1,1];    xlab = "lon"
        y = Data.lat.val[1,:,1];    ylab = "lat"
        z = Data.depth.val[1,1,:];  zlab = "depth [km]"
        orthogonal = true;
    elseif isa(Data,CartData)
        # Determine 
        x = Data.x.val[:,1,1];    xlab = "X [km]"
        y = Data.y.val[1,:,1];    ylab = "Y [km]"
        z = Data.z.val[1,1,:];  zlab = "Z [km]"
        axis_equal = true

        if sum(abs.(Data.x.val[:,1,2] - Data.x.val[:,1,1]))>1e-10
            orthogonal = false
            warning("Non-orthogonal CartData - can only show topography")
        else
            orthogonal = true
        end

    else
        error("not yet implemented ")
    end
    
    if !axis_equal
        x_vec = 1:length(x)
        y_vec = 1:length(y)
        z_vec = 1:length(z)
    else
        x_vec = range(x[1],x[end],length(x));
        y_vec = range(y[1],y[end],length(y));
        z_vec = range(z[1],z[end],length(z));
    end
    # determine width of axis
    dx = (maximum(x) - minimum(x))/(length(x)-1);
    dy = (maximum(y) - minimum(y))/(length(y)-1);
    dz = (maximum(z) - minimum(z))/(length(z)-1);
    
    if !isnothing(Topography)
        if isa(Topography,GeoData)
            x_topo = (Topography.lon.val .- x[1])/dx; 
            y_topo = (Topography.lat.val .- y[1])/dy; 
            z_topo = (Topography.depth.val .- z[1])/dz; 
        elseif isa(Topography,CartData)
            x_topo = Topography.x.val; 
            y_topo = Topography.y.val; 
            z_topo = Topography.z.val; 
        end

        if isnothing(Topo_range)
            # base 
            topo_max = round(maximum(ustrip.(Topography.fields.Topography)),digits=1)
            Topo_range = (-topo_max, topo_max)
        end

    end
    
    data_names  = keys(Data.fields)             # Names of the fields
    data_selected = Observable(Symbol(data_names[1]))
    data_string   = @lift String($data_selected)
    get_vol(f_name) = Data.fields[f_name]
    vol = lift(get_vol, data_selected)              

    fig = Figure(resolution = (2000,2000), fontsize=20)
    ax = LScene(fig[1, 1:2], scenekw = (camera = cam3d!, raw = false))

    # Create sliders
    sgrid = SliderGrid(
    fig[2, 2],
        (label = xlab, range = 1:length(x_vec),
            format = v -> string(round((v-1)*dx + x[1], digits = 2))),
        (label = ylab, range = 1:length(y_vec),
            format = v -> string(round((v-1)*dy + y[1], digits = 2))),
        (label = zlab, range = 1:length(z_vec),
            format = v -> string(round((v-1)*dz + z[1], digits = 2))),
        (label = "Transparency topo", range = 0:.01:1),
    )
    
    # Create dropdown menus
    menu_dataset  = Menu(fig, options = [String.(data_names)...], default=String(data_selected[]))
    menu_colormap = Menu(fig, options = ["roma","romaO","vik","turku","davos","batlow","tab10","tab20","bone"], 
                               default="roma")

    # Colorbar limits
    cmin = @lift round(minimum($vol), digits=2)
    cmax = @lift round(maximum($vol), digits=2)
    cmin_str = @lift string($cmin)
    cmax_str = @lift string($cmax)
    
    cmin_box    = Textbox(fig, stored_string = cmin_str,width = nothing)
    cmax_box    = Textbox(fig, stored_string = cmax_str,width = nothing)

    iso_level = Observable([1.7])
    iso_alpha = Observable(0.5);
    iso_box     = Textbox(fig, stored_string ="$(iso_level[][1])",width = nothing)
    iso_toggle  = Toggle(fig, active = true); 
    iso_slide   = Slider(fig, range = 0:.01:1)
    set_close_to!(iso_slide, iso_alpha[])

    Label(fig[3,1:4], " ",  width = nothing)

    fig[2, 1] = vgrid!(
        hgrid!(Label(fig, "Dataset",  width = nothing),menu_dataset),
        hgrid!(Label(fig, "Colormap", width = nothing),menu_colormap),
        hgrid!(Label(fig, "Color axis limits", width = nothing), hgrid!(cmin_box, Label(fig, "-", width = 20), cmax_box)),
        hgrid!(hgrid!(Label(fig, "Isovalue", width = nothing), iso_toggle), hgrid!(iso_box,Label(fig, "α: ", width = nothing),iso_slide));
        tellheight = false)


    lo = sgrid.layout
    nc = ncols(lo)
 
    # Note: volumeslices & GLMakie in general seems to have a bit of an issue with 
    # using real coordinates. In many cases the numerical values of lon/lat are much smaller than the depth values,
    # & not centered around zero.
    #
    # The "trick" we do here is to create the axis based on the number of points in the 3D volume
    # and simply overwrite the names of the labels

    if orthogonal
        plt = volumeslices!(ax, x_vec,y_vec,z_vec,vol, colorrange=(cmin,cmax), colormap=:roma)
        iso = GLMakie.contour!(plt, x_vec, y_vec, z_vec, vol, levels = iso_level, alpha=iso_alpha, colormap=plt.attributes.colormap, colorrange=plt.attributes.colorrange)

#        plt = volumeslices!(ax, reverse(x_vec),reverse(y_vec),z_vec,vol, colorrange=(cmin,cmax), colormap=:roma)
#        iso = GLMakie.contour!(plt, reverse(x_vec), reverse(y_vec), z_vec, vol, levels = iso_level, alpha=iso_alpha, colormap=plt.attributes.colormap, colorrange=plt.attributes.colorrange)
    end

    topo_alpha = Observable(0.5)
    if !isnothing(Topography)

         # in case topography is supplied
        topo_surf = surface!(ax, x_topo[:,:,1], y_topo[:,:,1], z_topo[:,:,1], colormap=(:oleron, topo_alpha[]), color=ustrip.(Topography.fields.Topography[:,:,1]), colorrange=Topo_range, transparency = true)
        cb_surf = Colorbar(fig[1, 4], topo_surf, vertical = true, label="Topography", height = Relative(0.6))

    end

    if !axis_equal
        xticks!(ax.scene, xtickrange=[0.;length(x)],xticklabels=["$(x[1])", "$(x[end])"])
        yticks!(ax.scene, ytickrange=[0.;length(y)],yticklabels=["$(y[1])", "$(y[end])"])
        zticks!(ax.scene, ztickrange=[0.;length(z)],zticklabels=["$(z[1])", "$(z[end])"])
    else
        xticks!(ax.scene, xtickrange=[x[1];x[end]],xticklabels=["$(x[1])", "$(x[end])"])
        yticks!(ax.scene, ytickrange=[y[1];y[end]],yticklabels=["$(y[1])", "$(y[end])"])
        zticks!(ax.scene, ztickrange=[z[1];z[end]],zticklabels=["$(z[1])", "$(z[end])"])
    end
    xlabel!(ax.scene, xlab)
    ylabel!(ax.scene, ylab)
    zlabel!(ax.scene, zlab)

    # 
    cb = Colorbar(fig[1, 3], plt, vertical = true, label=data_string, height = Relative(0.6))

    # connect sliders to `volumeslices` update methods
    sl_yz, sl_xz, sl_xy, sl_alpha_topo = sgrid.sliders

    on(sl_yz.value) do v; plt[:update_yz][](v) end
    on(sl_xz.value) do v; plt[:update_xz][](v) end
    on(sl_xy.value) do v; plt[:update_xy][](v) end

    if orthogonal
        set_close_to!(sl_yz, .5length(x_vec))
        set_close_to!(sl_xz, .5length(y_vec))
        set_close_to!(sl_xy, .5length(z_vec))
    end
    set_close_to!(sl_alpha_topo, .5)
    
    # change color limits 
    on(cmin_box.stored_string) do s
        ra = plt[:colorrange]
        plt[:colorrange] = (parse(Float64,s), ra.val[2])
    end
    on(cmax_box.stored_string) do s
        ra = plt[:colorrange]
        plt[:colorrange] = (ra.val[1], parse(Float64,s))
    end

    # Change data 
    on(menu_dataset.selection) do s
        data_selected[] = Symbol(s)
        plt[:colorrange] = (cmin[], cmax[])
        cmin_box.displayed_string = cmin_str[]
        cmax_box.displayed_string = cmax_str[]

        # update values
        set_close_to!(sl_yz,  sl_yz.value[])
        set_close_to!(sl_xz,  sl_xz.value[])
        set_close_to!(sl_xy,  sl_xy.value[])
        
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
        iso_level[] = [parse(Float64,s)] 
    end
    on(iso_slide.value) do v
        iso_alpha[] = v
    end

    on(sl_alpha_topo.value) do v
        topo_alpha[] = v
        if !isnothing(Topography)
            # in case topography is supplied
            topo_surf.attributes.colormap=(:oleron,topo_alpha[])
        end

    end

    # add toggles to show/hide slices
    if orthogonal
        hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
        toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:4]

        if !isnothing(Topography)
            map(zip([hmaps; topo_surf], toggles)) do (h, t)
                connect!(h.visible, t.active)
            end
        else
            map(zip(hmaps, toggles[1:3])) do (h, t)
                connect!(h.visible, t.active)
            end
        end
    end

    # axis data
    ax3 = ax.scene.plots[1]


    display(fig)
    
    return nothing
end


end