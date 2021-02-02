function [theta,phi,R_phi,R_theta]=Rmatrix(u_tcs)
%  Extract polar angles and rotation matrices from unit tcs vector
%
theta=atan2(u_tcs(3),sqrt(u_tcs(1)^2+u_tcs(2)^2));
phi  =atan2(u_tcs(2),u_tcs(1));
cth=cos(theta); sth=sin(theta);
cph=cos(phi); sph=sin(phi);
R_phi  =[cph, sph,      0; -sph, cph, 0,;       0, 0,      1];
R_theta=[cth,      0, sth;       0,      1, 0,; -sth, 0, cth];
end