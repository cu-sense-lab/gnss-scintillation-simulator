function [Coords, position,e,a]=ellipse(Foci,majorax,npt)
%
%  Generate a set of uniformly spaced points around the circumference of an ellipse
%
%   Foci         = each column contains (x,y) coordinates of a focus [2,2]
%   majorax      = length of major axis
%   npt          = number of points around the circumference (odd=> identical start & end points)
%
%   Coords       = (x,y) cooordinates of points [2,npt];
%   position     = circumferential measure of sample points
%                  0 is on the major axis near the first focus
%
% Chris Wilson 1 February 1999
%
focalvec = Foci(:,2)-Foci(:,1);
center   = (Foci(:,2)+Foci(:,1))/2;
rotang   = atan2(focalvec(2), focalvec(1));
Rot      = [cos(rotang), -sin(rotang); sin(rotang), cos(rotang)];
rmin     = sqrt(sum(focalvec.^2));
e        = rmin/majorax;
a        = majorax/2;

ntab =1001;
phi  = pi*linspace(0,1,ntab);
dphi = phi(2)-phi(1);
dldp = a*(1-e^2)*sqrt(1+e^2+2*e*cos(phi))./(1+e*cos(phi)).^2;
arclen = 0.5*dphi*[0,cumsum(dldp(1:end-1))+cumsum(dldp(2:end))];

N=round(npt);
position = arclen(end)*linspace(-1,1,2*floor(N/2)+1);
position = position(end-N+1:end);
intphi = interp1(arclen,phi,abs(position),'cubic').*sign(position);

rho = a*(1-e^2)./(1+e*cos(intphi));
Cords1 = [rho.*cos(intphi)+a*e; rho.*sin(intphi)];
Coords = Rot*Cords1+center*ones(1,N);
return

