% Returns llh coordinates of h_intercept along ray to satellite defined by u_sat                                          
%                                          
%                                          
% Input(s)                                 
%   h_intercept: Desired intercept height (e.g. middle of ionosphere) (m)
%   u_sat:       Unit vector in direction of satellite (from observer)
%   rng_sat:     Range from observer to satellite (m)
%   origin_llh:  Location of observer [lat,lon,alt], (deg, deg, m)
%                                          
% Output(s)                                
%   llh_test:    llh coordinates of h_intercept along ray to satellite 
%                defined by u_sat                                       
%                                          
% Example:  
%{
h_intercept = 300e3;  % 300 km intercept altitude
u_sat = [0.4501, -0.8930, -0.0032];
rng_sat = 3.0980e+06;
origin_llh = [0.7438   -1.2478  417.0000];
llh_test = findIntercept1(h_intercept,u_sat,rng_sat,origin_llh)
llh_test = 
  [4.6546e-01
  -1.0967e+00
   3.0012e+05]
%}
%                                          
%                                          
% See also: satGEOM, findIntercept                               
%                                          
%                                          
% Dependencies: llh2tcsT, tcs2llhT, ecf2uvwT, llh2ecfT, uvw2tcsT, tcs2uvwT, 
%               Get_DirCos_ForwardT, EarthModel
                                           
%                                          
%                                          
% Written by: John Peach 26-Sep-2014       
% MIT Lincoln Laboratory                   
%                                          
% Revisions:                               
                                           
                                           
function llh_test = findIntercept(h_intercept,u_sat,rng_sat,origin_llh)

% Earth radius
%refEllipsoid = referenceEllipsoid('earth','m');
%R = rsphere('curve',refEllipsoid,origin_llh(1));
[a,f] = EarthModel;
b=a-a*f;
R=(a+b)/2;

% Satellite position in TCS
sat_tcs = repmat(rng_sat,3,1).*u_sat;

% Satellite position in lat, long, alt
sat_llh = tcs2llhT(sat_tcs,origin_llh);
h_sat   = sat_llh(3,:);

% Angle between origin-to-satellite and origin-to-Earth center using law of
% cosines: (R + h_sat)^2 = R^2 + rng_sat^2 - 2*R^2*rng_sat^2*cos(gamma)
cos_gamma = (rng_sat.^2 - 2*R*h_sat - h_sat.^2)./(2*R*rng_sat);

% Use law of cosines again to find intercept range, r_test
% (R + h_intercept)^2 = R^2 + r_test^2 - 2*R*r_test*cos(gamma)
% r_test^2 - 2*R*cos(gamma)*r_test - (2*R*h_intercept + h_intercept^2) = 0
nsamps=length(cos_gamma);
r_test=zeros(size(cos_gamma));
for nsamp=1:nsamps
  r_test(nsamp) = max(roots([1, -2*R*cos_gamma(nsamp), -(2*R*h_intercept + h_intercept^2)]));  
end
% Location of intercept point in TCS coordinates
r_tcs =repmat(r_test,3,1).*u_sat;
% Convert to lat, long, alt
llh_test = tcs2llhT(r_tcs,origin_llh);
