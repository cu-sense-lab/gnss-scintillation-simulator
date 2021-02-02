function [R,PHI,x,y,Surf]=SurfaceTCSmap(latitude,longitude,RADAR)
%
%  Generate uniformly spaced topocentric x-y maps of surface height, range, and bearing
%  The WGS84 ellipsoid defines surface height
%
%  INPUTS:
%  latitude = uniformly spaced grid of latitudes   (rad +North)
%  longitude= uniformaly spaced grid of longitudes (rad +East)
%  RADAR    = structure containing Rlat, Rlon, Rheight (rad, rad, m)
%
%  RETURN:
%  R(x,y)   = Radar range to surface at (x,y)   (m)
%  PHI(x,y) = Radar azimuth to surface at (x,y) (rad)
%  Surf(x,y)= Surface height from topocentric tangent plane at radar origin (m)
%
dtr=pi/180; 
Rlat=RADAR.Rlat;
Rlon=RADAR.Rlon;
Rheight=RADAR.Rheight;
origin=[Rlat,Rlon,0]';   %TCL reference on ellipsoid at radar location
[LAT,LON]=meshgrid(latitude,longitude);
nptsx=length(longitude);
nptsy=length(latitude);
%Generate llh 3-by-nptsx*nptsy vector of surface (zero height) points
llh=[LAT(:) LON(:) zeros(nptsx*nptsy,1)]';
%Convert llh to tcs
tcs= llh2tcsT(llh,origin); 
clear llh LAT LON

x   = reshape(tcs(1,:),nptsx,nptsy); 
y   = reshape(tcs(2,:),nptsx,nptsy);
z   = reshape(tcs(3,:),nptsx,nptsy); clear tcs
%Compute range & azimuth(bearing) to each grid point
R   = sqrt(x.^2+y.^2+(z-Rheight).^2);
PHI = atan2(y,x); 
%Resample to uniform grid (Trim grid a bit to avoid NaN's in output
xs = linspace(ceil(min(x(:)/1000)),floor(max(x(:)/1000)-1),nptsx)*1000;
ys = linspace(ceil(min(y(:)/1000)),floor(max(y(:)/1000)-1),nptsy)*1000;
[xm,ym]=meshgrid(xs,ys);
Surf    = griddata(x,y,z,  xm,ym); clear z
R       = griddata(x,y,R,  xm,ym);
PHI     = griddata(x,y,PHI,xm,ym);
clear xm ym
x=xs; y=ys;
return