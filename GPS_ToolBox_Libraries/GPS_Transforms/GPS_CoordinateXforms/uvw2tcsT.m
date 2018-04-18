function tcs = uvw2tcs(uvw, origin)
%
% Description:
%	This function will convert a position vector from
%	UVW to TCS coordinates relative to origin.
%
% Usage:
%  tcs = uvw2tcs ( uvw , origin)
%
% Inputs:
%  uvw      3xN array of vectors in geocentric coordinates
%
%           x1 x2 x3 ... xN
%     uvw = y1 y2 y3 ... yN
%           z1 z2 z3 ... zN
%
%     origin 3x1 or 3xN vector array (CLR August 2005)
%
% Outputs:
%  tcs      3xN array of vectors in topocentric coordinates 
%
if isvector(origin)
    origin = origin(:);
end
if isvector(uvw)
    uvw = uvw(:);
end

origin_ECF = llh2ecfT(origin);
origin_UVW = ecf2uvwT(origin_ECF,origin);

YAW_TYPE = 1;
PITCH_TYPE = 2;
ROLL_TYPE = 3;

DC1 = Get_DirCos_ForwardT (pi/2, ROLL_TYPE);
DC2 = Get_DirCos_ForwardT (pi/2, PITCH_TYPE);
DC3 = Get_DirCos_ForwardT (-origin(1,:), ROLL_TYPE);
DC21=DC2*DC1;
sz=size(DC3);
if length(sz)==2
    DC = DC3*DC21;
else
    DC=zeros(3,3,sz(3));
    for n=1:sz(3)
        DC(:,:,n)=DC3(:,:,n)*DC21;
    end
end
    
szuvw = size(uvw);
tcs = zeros(szuvw(1),szuvw(2));
tcs(1,:) = uvw(1,:) - origin_UVW(1,:);
tcs(2,:) = uvw(2,:) - origin_UVW(2,:);
tcs(3,:) = uvw(3,:) - origin_UVW(3,:);
if length(sz)==2
    tcs = DC * tcs;
else
    for n=1:sz(3)
        tcs(:,n)=DC(:,:,n)*tcs(:,n);
    end
end

return