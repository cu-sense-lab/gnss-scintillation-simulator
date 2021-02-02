function [satXYZ,satXYZdot]= satposvel(t,eph)
%USAGE:  [satXYZ,satXYZdot]= satposvel(t,eph)
%GPS satellite ECF position and velocity versus time (gps time from epoch)
%INPUTS:  
%eph=21 element RINEX format GPS ephemeris
%tk=time since epoch seconds (1XN)
%OUTPUTS:
%satXYZ    = ECS coordinates (3XN)
%satXYZdot = ECS velocity    (3XN)
%
%Adapted from cited MatLaB software
%Charles Rino
%Rino Consulting
%January 20, 2014
%
nsamp=length(t);
satXYZ=zeros(3,nsamp);
satXYZdot=zeros(3,nsamp);

GM = 3.986005e14;             %Universal gravitational parameter m^3/s^2
Omegae_dot = 7.2921151467e-5; %Earth rotation rate, rad/s
%  Units are seconds, meters, and radians
%  Assigning the local variables to eph

% svprn   =   eph(1);           %Not used
% af2	    =   eph(2);           %Not used
M0	    =   eph(3);
roota   =   eph(4);
deltan  =   eph(5);
ecc	    =   eph(6);
omega   =   eph(7);
cuc	    =   eph(8);
cus	    =   eph(9);
crc	    =   eph(10);
crs	    =   eph(11); 
i0	    =   eph(12);
idot    =   eph(13);
cic	    =   eph(14);
cis	    =   eph(15);
Omega0  =   eph(16);
Omegadot=   eph(17);
toe	    =   eph(18);     
%af0	=   eph(19);     %Not used
%af1	=   eph(20);     %Not used
%toc	=   eph(21);     %Not used

%Position calculation from Strang & Borre Ch 15.2
A = roota*roota;
tk = check_t(t-toe);
n0 = sqrt(GM/A^3);
n = n0+deltan;
M = M0+n*tk;
mkdot = n;
E = M;

error = 1;  
nct=0;
while error > 1e-12
    tempE = M+ecc*sin(E);
    error = abs(tempE-E);
    E = tempE;
    nct=nct+1;
    if nct>20
        break
    end
end
ekdot = mkdot./(1.0 - ecc*cos(E));
tak = atan2( sqrt(1.0-ecc*ecc)*sin(E), cos(E)-ecc);
takdot = sin(E).*ekdot.*(1.0+ecc*cos(tak))./(sin(tak).*(1.0-ecc*cos(E)));

phik = tak+omega;
corr_u = cus*sin(2.0*phik) + cuc*cos(2.0*phik);
corr_r = crs*sin(2.0*phik) + crc*cos(2.0*phik);
corr_i = cis*sin(2.0*phik) + cic*cos(2.0*phik);
u = phik + corr_u;
r = A*(1-ecc*cos(E)) + corr_r;
ik = i0+idot*tk + corr_i;

ukdot = takdot +2.0*(cus*cos(2.0*u)-cuc*sin(2.0*u)).*takdot;
rkdot = A*ecc*sin(E)*n./(1.0-ecc*cos(E))+2.0*(crs*cos(2.0*u)-crc*sin(2.0*u)).*takdot;
ikdot = idot + (cis*cos(2.0*u)-cic*sin(2.0*u))*2.0.*takdot;

x1 = cos(u).*r;
y1 = sin(u).*r;

xpkdot = rkdot.*cos(u) - y1.*ukdot;
ypkdot = rkdot.*sin(u) + x1.*ukdot;

Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
omegakdot = (Omegadot-Omegae_dot);

satXYZ(1,:) = (x1.*cos(Omega)-y1.*cos(ik).*sin(Omega))';
satXYZ(2,:) = (x1.*sin(Omega)+y1.*cos(ik).*cos(Omega))';
satXYZ(3,:) = (y1.*sin(ik))';

satXYZdot(1,:)= (( xpkdot-y1.*cos(ik)*omegakdot ).*cos(Omega)...
    - ( x1.*omegakdot+ypkdot.*cos(ik)-y1.*sin(ik).*ikdot ).*sin(Omega))';
satXYZdot(2,:) = (( xpkdot-y1.*cos(ik)*omegakdot ).*sin(Omega)...
    + ( x1.*omegakdot+ypkdot.*cos(ik)-y1.*sin(ik).*ikdot ).*cos(Omega))';
satXYZdot(3,:) = (ypkdot.*sin(ik) + y1.*cos(ik).*ikdot)';


return