function ecf=llh2ecfT(llh)
%
% Description:
% 		This function returns the x, y, and z earth centered fixed (ECF)
%		coordinates of the point specified by llh [lat lon hgt].  
%		Note that longitude is positive east of the Greenwich meridian.
%
% Usage:
%   ecf=llh2ecfT(llh)
%
% Inputs:
%   llh   - 3xN Matrix of Vectors
%           where
%              llh(1,:) is latitude positive (radians)
%              llh(2,:) is longitude positive East (radians)
%              llh(3,:) is height (meters)
%
% Outputs
%   ecf   - 3xN Matrix of Vectors in ECF coordinates (meters)
%           where
%              x(1,:) is Greenwich meridan; (0 lon, 0 lat)
%              y(2,:) is Easterly
%              z(3,:) is North Pole
%

if isvector(llh)
    llh = llh(:);
end

lat=llh(1,:);
lon=llh(2,:);
hgt=llh(3,:);

%  Set up WGS-84 constants.
[a,f] = EarthModel;

%  Convert lat,lon,hgt geographic coords to XYZ Earth Centered Earth Fixed coords.
%		N = a/sqrt( 1 - f*(2-f)*sin(lat)*sin(lat) )
%		X = (N + h)*cos(lat)*cos(lon)
%		Y = (N + h)*cos(lat)*sin(lon)
%		Z = ((1-f)^2 * N + h)*sin(lat)

%  Store some commonly used values.

slat = sin(lat);
N = a./sqrt(1 - f*(2-f) * slat.^2);
Nplushgtclat = (N + hgt) .* cos(lat);

x = Nplushgtclat .* cos(lon);
y = Nplushgtclat .* sin(lon);
z = ((1-f)^2 * N + hgt) .* slat;

ecf = [x; y; z];

return
