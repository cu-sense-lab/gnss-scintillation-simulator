function [tcs]=llh2tcsT(llh,origin);
%
% Description:
% 	This function returns the x, y, and z topocentric (TCS)
%	coordinates of the point specified by llh [lat lon hgt],
%   relative to the input origin [lat lon alt].
%
% Usage:
%   tcs=llh2tcs_mks(llh,origin);
%
% Inputs:
%   llh   - 3xN Matrix of Vectors in LLH coordinates
%           where
%              llh(1,:) is latitude   (radians)
%              llh(2,:) is longitude  (radians)
%              llh(3,:) is height     (meters)
%   origin - 3xN Matrix of Vectors in LLH coordinates
%           where
%              origin(1,:) is latitude   (radians)
%              origin(2,:) is longitude  (radians)
%              origin(3,:) is height     (meters)
%
% Outputs
%   tcs   - 3XN Matrix of Vectors in TCS coordinates (meters)
%           where
%              tcs(1,:) is x
%              tcs(2,:) is y
%              tcs(3,:) is z
%


ecf = llh2ecfT(llh);

uvw = ecf2uvwT(ecf,origin);

tcs = uvw2tcsT(uvw,origin);

return
