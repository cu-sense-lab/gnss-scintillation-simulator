function          [rootSDF]=root_phaseSDF2(Cp,gnu1,gnu2,q0,Ky,Kz,SF,A,B,C)
%USAGE:     [rootSDF]=root_phaseSDF2(Cp,gnu1,gnu2,q0,Ky,Kz,SF,A,B,C)
%
%Generate root_rfn for 2D SDF
%
% Input parameters:
%   Cs       = Turbulent strength Q~Cs*q^(-2*gnu-1)  NOTE: Cs<=a*b*sec(theta)^2*Cs
%                                                          for consistent scaling
%   gnu1,gnu2= Spectral index parameter for each power-law segment (Kolmogorov=4/3)
%   q0             = Spatial wavenumber at transition
%   SF            = Scale factor for normalization (e.g. SF=dKy*dKz/(2*pi)^2)
%   ky,kz         = mesgrid values of ky, kz
%   A,B,C       = anisotropy/oblique-propagation coefficients
%

K    =sqrt(A*Ky.^2+B*Ky.*Kz+C*Kz.^2);
rootSDF=ones(size(K));
if gnu1~=gnu2
    Cp2=Cp*q0^(-2*(gnu1-gnu2));
    nSeg1=find(K<=q0);
    nSeg2=find(K>q0);
    rootSDF(nSeg1)=sqrt( Cp*SF*K(nSeg1).^(-2*gnu1-1));
    rootSDF(nSeg2)=sqrt(Cp2*SF*K(nSeg2).^(-2*gnu2-1));
else
    rootSDF=rootSDF.*sqrt(Cp*SF*K.^(-2*gnu1-1));
end
locINF= rootSDF==inf;
rootSDF(locINF)=0;
return


