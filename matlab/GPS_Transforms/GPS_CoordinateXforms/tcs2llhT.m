function [llh]=tcs2llhT(tcs,origin);
%
% Description:
% 	This function returns the LLH coordinates of the point specified by 
%   x,y,z topocentric (TCS) coordinates relative to the input origin.
%
% Usage:
%   llh=tcs2llhT(tcs,origin);
%
% Inputs:
%   tcs   - Nx3 Matrix of Vectors in TCS coordinates (meters)
%           where
%              tcs(1,:) is x
%              tcs(2,:) is y
%              tcs(3,:) is z
%
% Outputs
%   llh   - Nx3 Matrix of Vectors in LLH coordinates
%           where
%              llh(1,:) is latitude   (radians)
%              llh(2,:) is longitude  (radians)
%              llh(3,:) is height     (meters)
%

if isvector(tcs)
    tcs = tcs(:);
end
if isvector(origin)
    origin = origin(:);
end

uvw=tcs2uvwT(tcs,origin);

llh=zeros(size(uvw));
LLH = uvw2llhT(uvw,origin);
llh(1,:)=LLH.lat;
llh(2,:)=LLH.lon;
llh(3,:)=LLH.hgt;

return
