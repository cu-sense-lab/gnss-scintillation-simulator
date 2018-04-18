function   [range,stheta,ctheta]=GHorizon(bearing,source,height)
%
%  Compute surface ellipsoid intercept of ray to surface
%
%  bearing  = bearing angle to surface                        [1xN] or 1                      
%  source   = source location llh                             [3xN] or 3x1                           
%  height   = height of surface above ellipsoid               m
% RETURNS:
%  range    = range to Geometric Horizon                      [1xN]              
%  s/ctheta = sine/cosine of depression angle                 [1xN]            
%  NOTE: v_tcs=[stheta'.*sin(bearing); stheta'.*cos(bearing); ctheta'];   3xN
%        Multiply by range to get tcs coordinates of horizon
%        Use tcs2llh to get geodetic coordinates of horizon
%
%WGS8 Ellipsoid
[a,f] = EarthModel;
a=(a+height);
b=a*(1-f);
aob2=(a/b)^2;
sz_b=size(bearing);
sz_s=size(source);
nsamp=max(sz_b(2),sz_s(2));
if sz_b(2)==1 
    bearing=repmat(bearing,1,nsamp);
end
if sz_s(2)==1
    source=repmat(source,3,nsamp);
end
    
OMGA=tcs2ecfR(source);                 %Rotation matrix tcs to ecf  
source_ecf=llh2ecfT(source);
cphi=cos(bearing(:));
sphi=sin(bearing(:));

A1=squeeze(OMGA(1,1,:)).*sphi+squeeze(OMGA(1,2,:)).*cphi;
A2=squeeze(OMGA(2,1,:)).*sphi+squeeze(OMGA(2,2,:)).*cphi;
A3=squeeze(OMGA(3,1,:)).*sphi+squeeze(OMGA(3,2,:)).*cphi;

alph= (A1.^2-squeeze(OMGA(1,3,:)).^2)+ ...
    (A2.^2-squeeze(OMGA(2,3,:).^2))+(A3.^2-squeeze(OMGA(3,3,:).^2))*aob2;
btea= 2*(A1.*squeeze(OMGA(1,3,:))+A2.*squeeze(OMGA(2,3,:))+A3.*squeeze(OMGA(3,3,:))*aob2);
gamm= squeeze(OMGA(1,3,:)).^2+squeeze(OMGA(2,3,:)).^2+squeeze(OMGA(3,3,:)).^2*aob2;
xi  = source_ecf(1,:)'.*A1+source_ecf(2,:)'.*A2+source_ecf(3,:)'.*A3*aob2;
eta = source_ecf(1,:)'.*squeeze(OMGA(1,3,:)) ...
    +source_ecf(2,:)'.*squeeze(OMGA(2,3,:))+source_ecf(3,:)'.*squeeze(OMGA(3,3,:))*aob2;
D2  = source_ecf(1,:).^2+source_ecf(2,:).^2+source_ecf(3,:).^2*aob2;
c0  = (D2-a^2);

nsamp=length(bearing);
ceps=zeros(nsamp,1);
for isamp=1:nsamp
    ceps(isamp)=fzero(@fq,[-1,0],optimset,alph(isamp),btea(isamp),...
                                     xi(isamp),eta(isamp),gamm(isamp),c0(isamp));
end
seps=sqrt(1-ceps.^2);
a0  = (alph.*ceps.^2-btea.*seps.*ceps+gamm);
range = -(xi.*ceps-eta.*seps)./a0;
stheta =-ceps;
ctheta = seps;

return

function y=fq(ceps,alph,btea,xi,eta,gamm,c0)
seps=sqrt(1-ceps.^2);
a0  = (alph.*ceps.^2-btea.*seps.*ceps+gamm);
boa  = (xi.*ceps-eta.*seps)./a0;
y=boa.^2-c0./a0;
return                          