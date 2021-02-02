function [ecf] = uvw2ecfT(UVW,origin)
%
% Description:
%   Rotate UVW into ECF: X axis (U axis) is colinear with longitude of origin.
%

if isvector(UVW)
    UVW = UVW(:);
end
if isvector(origin)
    origin = origin(:);
end

YAW_TYPE = 1;

DC = Get_DirCos_ForwardT(-origin(2,:),YAW_TYPE);

sz=size(DC);
if length(sz)==3
    ecf = zeros(3,sz(3));
    for n=1:sz(3)
        ecf(:,n)=DC(:,:,n)*UVW(:,n);
    end
else
    ecf = DC*UVW;
end

return