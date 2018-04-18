function uvw = tcs2uvwT( tcs, origin)
%
% Description:
% 		This function will rotate TCS coordinates into UVW coordinates:
%       X axis (U axis) is colinear with longitude of origin.
%
% Usage:
%   uvw = tcs2uvwT(tcs, origin)
%
% Inputs:
%   tcs	3xN array of vectors
%   origin 3x1 or 3XN CLR August, 2005
%
% Outputs:
%   uvw	3xN array of vectors
%
if isvector(tcs)
    tcs = tcs(:);
end
if isvector(origin)
    origin = origin(:);
end

YAW_TYPE = 1;
PITCH_TYPE = 2;
ROLL_TYPE = 3;

origin_ECF = llh2ecfT(origin);             %origin_ECF is Nx3
origin_UVW = ecf2uvwT(origin_ECF,origin);

DC1 = Get_DirCos_ForwardT (pi/2, ROLL_TYPE);
DC2 = Get_DirCos_ForwardT (pi/2, PITCH_TYPE);
DC3 = Get_DirCos_ForwardT (-origin(1,:), ROLL_TYPE);
DC21=DC2*DC1;
sz=size(DC3);
if length(sz)==2
    DC = DC3*DC21;
    DCINV=DC^(-1);
    uvw = DCINV * tcs;
else
    DCINV=zeros(3,3,sz(3));
    uvw  =zeros(3,sz(3));
    for n=1:sz(3)
        DCINV(:,:,n)=DC3(:,:,n)*DC21;
        DCINV(:,:,n)=DCINV(:,:,n)^(-1);
        uvw(:,n) = DCINV(:,:,n) * tcs(:,n);
    end
end

uvw(1,:) = uvw(1,:) + origin_UVW(1,:);
uvw(2,:) = uvw(2,:) + origin_UVW(2,:);
uvw(3,:) = uvw(3,:) + origin_UVW(3,:);

return

% Notes:
% In the file uvw2tcs.m, we have the following transform for UVW --> TCS
%
% TCS = DC * ( UVW - origin_UVW )
%
% To get the inverse transform for TCS --> UVW,
%	DC^(-1) * TCS = DC^(-1)*DC * ( UVW - origin_UVW )
%	DC^(-1)*TCS = UVW - origin_UVW
%	UVW = DC^(-1)*TCS + origin_UVW
