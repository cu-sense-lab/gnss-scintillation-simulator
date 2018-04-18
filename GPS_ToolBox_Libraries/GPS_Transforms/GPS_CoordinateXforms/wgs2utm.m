function  [x,y,utmzone] = wgs2utm(Lat,Lon)
% -------------------------------------------------------------------------
% [x,y,utmzone] = wgs2utm(Lat,Lon)
%
% Description:
%    Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
%    (northing, easting).
%
% Input:
%    Lat: WGS84 Latitude scalar or vector in decimal degrees  
%    Lon: WGS84 Longitude scalar or vector in decimal degrees  
%
% Output:
%    x: UTM easting in meters
%    y: UTM northing in meters
%    utmzone: UTM longitudinal zone
%
% Author notes:
%    I downloaded and tried deg2utm.m from Rafael Palacios but found
%    differences of up to 1m with my reference converters in southern
%    hemisphere so I wrote my own code based on "Map Projections - A
%    Working Manual" by J.P. Snyder (1987). Quick quality control performed
%    only by comparing with LINZ converter
%    (www.linz.govt.nz/apps/coordinateconversions/) and Chuck Taylor's 
%    (http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html) on a 
%    few test points, so use results with caution. Equations not adapted
%    for a latitude of +/- 90deg.
%
% Example:
%    Lat=[48.866667; 34.05;   -36.85 ];
%    Lon=[2.333056;  -118.25; 174.783333];
%    [x,y,utmzone] = wgs2utm(Lat,Lon)
%
% Author: 
%   Alexandre Schimel
%   Coastal Marine Group - University of Waikato
%   Hamilton, New Zealand
%
% Version:
%   April 2007
%-------------------------------------------------------------------------

%% Argument checking
error(nargchk(2, 2, nargin));  %2 arguments required
n1=size(Lat);
n2=size(Lon);
if (n1~=n2)
   error('Lat and Lon should have same size');return
end

%% coordinates in radians
lat = Lat.*pi./180;
lon = Lon.*pi./180;

%% WGS84 parameters
a = 6378137;           %semi-major axis
b = 6356752.314245;    %semi-minor axis
% b = 6356752.314140;  %GRS80 value, originally used for WGS84 before refinements
e = sqrt(1-(b./a).^2); % eccentricity

%% UTM parameters
% lat0 = 0;                % reference latitude, not used here
Lon0 = floor(Lon./6).*6+3; % reference longitude in degrees
lon0 = Lon0.*pi./180;      % in radians
k0 = 0.9996;               % scale on central meridian

FE = 500000;              % false easting
FN = (Lat < 0).*10000000; % false northing 

%% Equations parameters
eps = e.^2./(1-e.^2);  % e prime square
% N: radius of curvature of the earth perpendicular to meridian plane
% Also, distance from point to polar axis
N = a./sqrt(1-e.^2.*sin(lat).^2); 
T = tan(lat).^2;                
C = ((e.^2)./(1-e.^2)).*(cos(lat)).^2;
A = (lon-lon0).*cos(lat);                            
% M: true distance along the central meridian from the equator to lat
M = a.*(  ( 1 - e.^2./4 - 3.*e.^4./64 - 5.*e.^6./256 )  .* lat         ...
         -( 3.*e.^2./8 + 3.*e.^4./32 + 45.*e.^6./1024 ) .* sin(2.*lat) ...
         +( 15.*e.^4./256 + 45.*e.^6./1024 )            .* sin(4.*lat) ...
         -(35.*e.^6./3072 )                             .* sin(6.*lat) );

%% easting
x = FE + k0.*N.*(                                  A       ...
                 + (1-T+C)                      .* A.^3./6 ...
                 + (5-18.*T+T.^2+72.*C-58.*eps) .* A.^5./120 );
               
%% northing 
% M(lat0) = 0 so not used in following formula
y = FN + k0.*M + k0.*N.*tan(lat).*(                                     A.^2./2  ...
                                   + (5-T+9.*C+4.*C.^2)              .* A.^4./24 ...
                                   + (61-58.*T+T.^2+600.*C-330.*eps) .* A.^6./720 );
                                 
%% UTM zone
utmzone = floor(Lon0./6)+31;