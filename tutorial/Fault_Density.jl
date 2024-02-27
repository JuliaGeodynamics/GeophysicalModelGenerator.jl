using GeophysicalModelGenerator, Shapefile, Plots, Rasters

# Data from "Active Faults of Eurasia Database AFEAD v2022" DOI:10.13140/RG.2.2.25509.58084
File   = "AFEAD_v2022/AFEAD_v2022/AFEAD_v2022.shp"

table  = Shapefile.Table(File)
geoms  = Shapefile.shapes(table)
CONF   = table.CONF

function selectshapes(table)
    geoms = empty(Shapefile.shapes(table))
    for row in table
        if row.CONF != "D"
            push!(geoms, Shapefile.shape(row))
        end
    end
    return geoms
end

ind = findall((table.CONF .== "A") .| (table.CONF .== "B"))
#fau = table.geometry[ind]

faults  = Shapefile.Handle(File).shapes
faultsR = rasterize(last,faults; res=(0.12,0.12), missingval=0, fill=1, atol = 0.4, shape=:line)

# plotting
p = plot(faultsR; color=:red)
plot!(p, faultsR; fillalpha=0, linewidth=0.6)

lon  = faultsR.dims[1]
lat  = faultsR.dims[2]

indlat = findall((lat .> 35) .& (lat .< 60))
Lat    = lat[indlat]
indlon = findall((lon .> -10) .& (lon .< 35))
Lon    = lon[indlon]

data = faultsR.data[indlon,indlat]

# plotting
p = plot(data; color=:red)
plot!(p, data; fillalpha=0, linewidth=0.6)

# create density map
stepsize = 4

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

heatmap(lonstep.val,latstep.val,countmap)

# save as GeoData
Lon3D,Lat3D, Dens3D = LonLatDepthGrid(lonstep,latstep,0);
Dens3D[:,:,1] = countmap'
Data_Faults = GeoData(Lon3D,Lat3D,Dens3D,(Faults=Dens3D,))
Write_Paraview(Data_Faults, "Faults_Spherical")

