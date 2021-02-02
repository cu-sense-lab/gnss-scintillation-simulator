function [a,f] = EarthModel
%
% Define the constants for the WGS-84 ellipsoidal Earth model.
%
% Outputs:
%   a - semi-major axis of the Earth ellipsoid model
%   f - flattening

a = 6378137.0; %meters
f = 1.0/298.257223563;

return