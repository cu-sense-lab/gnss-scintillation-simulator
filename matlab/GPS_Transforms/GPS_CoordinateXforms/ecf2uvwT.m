function uvw = ecf2uvwT( ecf, origin )
%
% Description: 
%    This function will rotate ECF coordinates into UVW coordinates:
%    X axis (U axis) is colinear with longitude of origin
%
% Usage:
%  uvw = ecf2uvw( ecf, origin )
%
% Inputs:
%  ecf         3xN array of vectors in ECF coordinates
%
%                   x1 x2 x3 ... xN
%              ecf= y1 y2 y3 ... yN
%                   z1 z2 z3 ... zN
%
% Outputs:
%  uvw         3xN array of vectors in UVW coordinates 
%

if isvector(origin)
    origin = origin(:);
end
if isvector(ecf)
    ecf = ecf(:);
end

YAW_TYPE = 1;
DC = Get_DirCos_ForwardT (origin(2,:), YAW_TYPE);

sz  = size(origin);
sze = size(ecf);

if sz(2)==1
    uvw = DC*ecf;
elseif sze(2)==1
   uvw = zeros(3,sz(2));
   for n=1:sz(2)
       uvw(:,n)=DC(:,:,n)*ecf;
   end
elseif sz(2)==sze(2)
   uvw = zeros(3,sz(2));
   for n=1:sz(2)
       uvw(:,n)=DC(:,:,n)*ecf(:,n);
   end
else
    error('Inconsistent dimensions');
end

return