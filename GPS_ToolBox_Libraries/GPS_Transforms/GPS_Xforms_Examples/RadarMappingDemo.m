%Demo coordinaten manipulations
clear
close all
%Map radar coordinates to surface using precomputed maps
dtr=pi/180;
if ~exist('R','var')
    %%%% Radar Location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rlat= 36.0;              
    Rlon=-75.0;            
    Rheight=30;              
    RADAR.Rlat=Rlat*dtr;
    RADAR.Rlon=Rlon*dtr;
    RADAR.Rheight=Rheight;
    nptsx=500;
    nptsy=200;
    %Lat-lon boundaries of region to be mapped
    lat_min= 35.9; lat_max= 36.1;
    lon_min=-74.9; lon_max=-75.1;
    %Define rectangular lat-long grid & mesh
    latitude  =linspace(lat_min,lat_max,nptsx)*dtr;
    longitude =linspace(lon_min,lon_max,nptsy)*dtr;
    [R,PHI,x,y,Surf]=SurfaceTCSmap(latitude,longitude,RADAR);
    if 1
    figure
    imagesc(x/1000,y/1000,R/1000)
    axis xy
    colorbar
    title('Range(km) to surface at x-y')
    xlabel('x--km')
    ylabel('y--km')

    figure
    imagesc(x/1000,y/1000,PHI/dtr)
    axis xy
    colorbar
    title('Bearing(deg) to surface at x-y')
    xlabel('x--km')
    ylabel('y--km')
    end
end