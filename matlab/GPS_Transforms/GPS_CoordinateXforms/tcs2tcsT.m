function [tcs2] = tcs2tcs(tcs1,origin1,origin2)
%
% Description:
%   Transform TCS relative to origin1 into ECF, 
%   then from ECF to TCS relative to origin2.
%

uvw1 = tcs2uvwT(tcs1,origin1);
ecf  = uvw2ecfT(uvw1,origin1);
uvw2 = ecf2uvwT(ecf,origin2);
tcs2 = uvw2tcsT(uvw2,origin2);

return
