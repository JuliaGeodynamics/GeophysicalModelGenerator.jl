"""
This is Tutorial_SeismicityTime.jl. It shows how you can create a movie
to represent your seismicity through time for Campi Flegrei (Italy).

You will need to download the zipped folder containing all files from:
(https://seafile.rlp.net/f/ff2c8424274c4d56b1f7/](https://ngdc.noaa.gov/mgg/global/global.html)
Remember to work in the downloaded directory.

"""


using DelimitedFiles, GeophysicalModelGenerator, Dates

# Seismicity movie
p2                  = @__FILE__;
p2last              = findlast("/",p2);
p3=chop(p2,head=0,tail = length(p2)-p2last[1]+1);
output_path         = string(p3,"/");
movie               = Movie_Paraview(name=string(p3,"/TemporalSeismicity"), Initialize=true);
if isdir(string(p3,"/TemporalSeismicity"))==0
    mkdir(string(p3,"/TemporalSeismicity"));
end

# 1. define where the two files are located on your computer
data                = readdlm("SeismicLocations/Seismicity_UTM_1983_1984.txt", '\t', skipstart=0, header=false);
l = length(data[:,1]);
dates = Date.(zeros(l,1))
for i = 1:l
    df                     = DateFormat("yyyymmdd");
    t1                     = string("19",@sprintf("%5.o", data[i,1]));
    t2                     = Date(t1[1:8],df);
    dates[i,1]              = t2;
end
WE                  = data[:,2];
SN                  = data[:,3];
depth               = data[:,4];
dates_num           = dates - dates[1];

# time steps
nt                  = 50;
dt                  = Dates.Day(nt);
t                   = minimum(dates):dt:maximum(dates);

for itime = 1:length(t)-1
    name            = string(p3,"/TemporalSeismicity/", string(itime));
    tt=findall(x->(x>=t[itime]) & (x<=t[itime+1]),dates);
    we             = WE[tt];
    sn             = SN[tt];
    Depth1          = depth[tt];
    DN              = dates_num[tt];
    label_time      = Dates.value(DN[end]);
    if size(tt,1)>1
        Data_set    = UTMData(we, sn, Depth1, 33, true, (Depth=Depth1*km,Timedata=DN));
        movie       = Write_Paraview(Data_set, name,pvd=movie,time=label_time,PointsData=true);
    end
end
Movie_Paraview(pvd=movie, Finalize=true)
