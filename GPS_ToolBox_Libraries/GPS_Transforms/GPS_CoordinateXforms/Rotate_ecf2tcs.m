function D=Rotate_ecf2tcs(origin)
%
%  Rotate free ecf vector into alignment with tcs
%       xyz=DXYZ
%
clat=cos(origin(1)); slat=sin(origin(1));
clon=cos(origin(2)); slon=sin(origin(2));
D=[-slon, clon, 0; -slat*clon, -slat*slon, clat; clat*clon, clat*slon, slat];
return
