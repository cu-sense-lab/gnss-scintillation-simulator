function   [sat_ecf]=satposYuma(t_sec,almanac)
%
u = 3.986005e+014;   %  WGS 84 value of the earth's gravitational  constant for GPS user
%  Compute mean motion n & mean anomaly M
n = sqrt(u/almanac.sqrtA^6);              %  Mean motion
toe = almanac.time + almanac.Af0;  %  Reference time of ephemeris

sat_ecf=zeros(3,length(t_sec));
M = almanac.meanAnom + n*(t_sec - toe);   % mean anomaly
%  Compute eccentricity anomaly E (iteratively)
E = M; % initialization
err = 10;
while (err > 1e-5)
    temp_E = M + almanac.e*sin(E);
    err = abs(temp_E - E);
    E = temp_E;
end
%  Compute orbit radius r and true anomaly v
r = almanac.sqrtA^2*(1 - almanac.e*cos(E));  %  Orbit radius
v = atan2(sqrt(1 - almanac.e^2)*sin(E)./(1 - almanac.e*cos(E)),(cos(E) - almanac.e)./(1 - almanac.e*cos(E)));
    
%  Compute SV coordinate in rotated orbit (in meters)
xp = r.*cos(almanac.argPG + v);
yp = r.*sin(almanac.argPG + v);

%   Convert SV position to ECEF coordinate (xt,yt,zt)
%   xt: SV ECEF x coordinate in meters
%   yt: SV ECEF y coordinate in meters
%   zt: SV ECEF z coordinate in meters
OmegaRate_earth = 7.2921151467e-005;  %  The earth's rotation rate
Omega_p = almanac.Omega + almanac.OmegaRate*(t_sec - toe) - OmegaRate_earth*t_sec;
xt = xp.*cos(Omega_p) - yp.*sin(Omega_p).*cos(almanac.incAng);
yt = xp.*sin(Omega_p) + yp.*cos(Omega_p).*cos(almanac.incAng);
zt = yp.*sin(almanac.incAng);
sat_ecf=[xt; yt; zt];
return
    