function llh=ecf2llhT(ecf)
%
% Description:
%   This function returns the [lat lon hgt] LLH coordinates
%   of the given earth center fixed (ECF) coordinates.
%
% Usage:
%   llh=ecf2llhT(ecf)
%
% Inputs
%   ecf   - 3xN Matrix of Vectors in ECF coordinates (meters)
%           where
%              x(1,:) is Greenwich meridan; (0 lon, 0 lat)
%              y(2,:) is Easterly
%              z(3,:) is North Pole
%
% Outputs:
%   llh   - 3xN Matrix of Vectors
%           where
%              llh(1,:) is latitude positive (radians)
%              llh(2,:) is longitude positive East (radians)
%              llh(3,:) is height (meters)
%

% Need an origin to go from ECF to UVW and UVW to LLH.
origin = [0 0 0];

uvw = ecf2uvwT(ecf,origin);

llh=zeros(size(uvw));
LLH = uvw2llhT(uvw,origin);
llh(1,:)=LLH.lat;
llh(2,:)=LLH.lon;
llh(3,:)=LLH.hgt;

return
