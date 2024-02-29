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
heatmap(lon.val,lat.val,coastlines',legend=false,colormap=cgrad(:gray1,rev=true),alpha=0.4);
plot!(faults; color=:red,legend = false,title="Fault Map World",ylabel="Lat",xlabel="Lon")

# restrict area
indlat = findall((lat .> 35) .& (lat .< 60))
Lat    = lat[indlat]
indlon = findall((lon .> -10) .& (lon .< 35))
Lon    = lon[indlon]
data   = faults.data[indlon,indlat]

# create GeoData from restricted data
Lon3D,Lat3D, Faults = LonLatDepthGrid(Lon,Lat,0);
Faults[:,:,1]       = data
Data_Faults         = GeoData(Lon3D,Lat3D,Faults,(Faults=Faults,))

steplon  = Int(length(Lon)/3)
steplat  = Int(length(Lon)/3)
countmap = CountMap(Data_Faults,"Faults",steplon,steplat)

lon = unique(countmap.lon.val)
lat = unique(countmap.lat.val)
coastlinesEurope = itp[lon,lat]
coastlinesEurope = map(y -> y > 1 ? 1 : y, coastlinesEurope)
heatmap(lon,lat,coastlinesEurope',colormap=cgrad(:gray1,rev=true),alpha=1.0);
heatmap!(lon,lat,countmap.fields.CountMap[:,:,1]',colormap=cgrad(:batlowW,rev=true),alpha = 0.8,legend=true,title="Fault Density Map Europe",ylabel="Lat",xlabel="Lon")

Write_Paraview(countmap, "Faults")



