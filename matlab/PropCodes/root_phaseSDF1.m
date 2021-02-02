function          [rootSDF]=root_phaseSDF1(Cp,p1,p2,q0,ky,SF)
%USAGE:     [rootSDF]=root_phaseSDF1(Cp,p1,p2,q0,ky,SF)
%
%Generate root_phaseSDF1 for 1D SDF
%
% Input parameters:
%   Cp       = Turbulent strength 
%   p1,p2= Spectral index parameter for each power-law segment (Kolmogorov=4/3)
%   q0             = Spatial wavenumber at transition
%   SF            = Scale factor for normalization (e.g. SF=dKy*dKz/(2*pi)^2)
%   ky              
%

rootSDF=ones(size(ky));
q=abs(ky);
if p1~=p2
    Cp2=Cp*q0^(-(p1-p2));
    nSeg1=find(q<=q0);
    nSeg2=find(q>q0);
    rootSDF(nSeg1)=sqrt( Cp*SF*q(nSeg1).^(-p1));
    rootSDF(nSeg2)=sqrt(Cp2*SF*q(nSeg2).^(-p2));
else
    rootSDF=rootSDF.*sqrt(Cp*SF*q.^(-p1));
end
locINF= rootSDF==inf;
rootSDF(locINF)=0;
return


