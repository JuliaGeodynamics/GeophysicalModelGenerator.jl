using GeophysicalModelGenerator, Shapefile, Plots, Rasters, GeoDatasets, Interpolations

# Data from "Active Faults of Eurasia Database AFEAD v2022" DOI:10.13140/RG.2.2.25509.58084
File   = "AFEAD_v2022/AFEAD_v2022/AFEAD_v2022.shp"

# load data from shapefiles
table  = Shapefile.Table(File)
geoms  = Shapefile.shapes(table)
CONF   = table.CONF

# select data and rasterize it
ind     = findall((table.CONF .== "A") .| (table.CONF .== "B") .| (table.CONF .== "C"))
faults  = Shapefile.Handle(File).shapes[ind]
faults  = rasterize(last,faults; res=(0.12,0.12), missingval=0, fill=1, atol = 0.4, shape=:line)

# get lon lat dimensions from data rasterized data
lon    = faults.dims[1]
lat    = faults.dims[2]

# download coastlines
lonC,latC,dataC = GeoDatasets.landseamask(;resolution='l',grid=10);
# interpolate to fault grid
itp        = LinearInterpolation((lonC, latC), dataC)
coastlines = itp[lon.val,lat.val]
coastlines = map(y -> y > 1 ? 1 : y, coastlines)

# plot the fault data
p = heatmap(lon.val,lat.val,coastlines',legend=false,colormap=cgrad(:gray1,rev=true),alpha=0.4);
p = plot!(faults; color=:red,legend = false,title="Fault Map World",ylabel="Lat",xlabel="Lon")

# restrict lon lat area
indlat = findall((lat .> 35) .& (lat .< 60))
Lat    = lat[indlat]
indlon = findall((lon .> -10) .& (lon .< 35))
Lon    = lon[indlon]
data   = faults.data[indlon,indlat]

ind_lon = findall( (lonC .> 0) .& (lonC .< 30 ) );
ind_lat = findall( (latC .> 35) .& (latC .< 50 ) );

# create density map
stepsize = 2
lonstep = Lon[1:stepsize:end]
latstep = Lat[1:stepsize:end]
countmap = zeros(length(latstep),length(lonstep))
for i = 1:length(lonstep)-1
    for j = 1:length(latstep)-1
        indi    = findall((Lon .>= lonstep[i]) .& (Lon .< lonstep[i+1]))
        indj    = findall((Lat .>= latstep[j]) .& (Lat .< latstep[j+1]))
        dataint = data[indi,indj]
        count   = sum(dataint)
        countmap[j,i] = count
    end
end
maxcount = maximum(countmap)
countmap = countmap ./ maxcount

coastlinesEurope = itp[lonstep.val,latstep.val]
coastlinesEurope = map(y -> y > 1 ? 1 : y, coastlinesEurope)
p2 = heatmap(lonstep.val,latstep.val,coastlinesEurope',colormap=cgrad(:gray1,rev=true),alpha=1.0);
p2 = heatmap!(lonstep.val,latstep.val,countmap,colormap=cgrad(:batlowW,rev=true),alpha = 0.7,legend=true,title="Fault Density Map Europe")

# save as GeoData
Lon3D,Lat3D, Dens3D = LonLatDepthGrid(lonstep,latstep,0);
Dens3D[:,:,1] = countmap'
Data_Faults = GeoData(Lon3D,Lat3D,Dens3D,(Faults=Dens3D,))
Write_Paraview(Data_Faults, "Faults_Spherical")



