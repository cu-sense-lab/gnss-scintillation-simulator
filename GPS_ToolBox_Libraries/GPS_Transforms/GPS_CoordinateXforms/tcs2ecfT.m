function [ecf] = tcs2ecfT(tcs,origin);
%
% Description:
%   Transform from TCS into ECF.

uvw = tcs2uvwT(tcs,origin);
ecf = uvw2ecfT(uvw,origin);

return