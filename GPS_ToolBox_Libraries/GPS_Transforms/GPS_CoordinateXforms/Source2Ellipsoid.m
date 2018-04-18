function   [stheta, ctheta]=Source2Ellipsoid(range,bearing,source,height)
%
%  Compute direction of ray from source to surface intercept at range
%  Reference coordinate system is tcs with origin at the source location
%  The surface is parallel to the WGS84 ellipsoid at height, which may be negative                     
%INPUTS
%  range   = range to surface                                [NX1] meters
%  bearing = bearing angle to surface                        [NX1] radians 
%  source  = source location llh                             [NX1] or [1X3]
%  height  = height of surface above ellipsoid               m
%RETURN
%  stheta, ctheta= sine and cosine of polar angle
%  NOTES: v_tcs=[stheta'.*sin(bearing); stheta'.*cos(bearing); ctheta'];   3xN
%        Multiply by range to get tcs coordinates of surface intercept
%        Use tcs2llh to get geodetic coordinates of surface intercept
%        In general there are 2 solutions the smaller of which is choosen:
%        Solutions beyond the geometric horizon where the two solutions converge
%        are legetimate.  Use GHorizon to find the range to the geometric horison
%

%WGS84 Ellipsoid
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
    source=repmat(source,1,nsamp);
end
    
OMGA=tcs2ecfR(source);                 %Rotation matrix tcs to ecf  
source_ecf=llh2ecfT(source);
cphi=cos(bearing(:));
sphi=sin(bearing(:));

A1=squeeze(OMGA(1,1,:)).*sphi+squeeze(OMGA(1,2,:)).*cphi;
A2=squeeze(OMGA(2,1,:)).*sphi+squeeze(OMGA(2,2,:)).*cphi;
A3=squeeze(OMGA(3,1,:)).*sphi+squeeze(OMGA(3,2,:)).*cphi;

alph= (A1.^2-squeeze(OMGA(1,3,:)).^2)+...
      (A2.^2-squeeze(OMGA(2,3,:).^2))+(A3.^2-squeeze(OMGA(3,3,:).^2))*aob2;
btea= 2*(A1.*squeeze(OMGA(1,3,:))+A2.*squeeze(OMGA(2,3,:))+A3.*squeeze(OMGA(3,3,:))*aob2);
gamm= squeeze(OMGA(1,3,:)).^2+squeeze(OMGA(2,3,:)).^2+squeeze(OMGA(3,3,:)).^2*aob2;
xi  = source_ecf(1,:)'.*A1+source_ecf(2,:)'.*A2+source_ecf(3,:)'.*A3*aob2;
eta = source_ecf(1,:)'.*squeeze(OMGA(1,3,:))...
     +source_ecf(2,:)'.*squeeze(OMGA(2,3,:))+source_ecf(3,:)'.*squeeze(OMGA(3,3,:))*aob2;
D2  = source_ecf(1,:).^2+source_ecf(2,:).^2+source_ecf(3,:).^2*aob2;
c0  = (D2-a^2);

nsamp=length(bearing);
ctheta=zeros(nsamp,1);
for isamp=1:nsamp
    ctheta(isamp)=fzero(@fq1,[-1,0],optimset,range(isamp),alph(isamp),btea(isamp),...
                        gamm(isamp),xi(isamp),eta(isamp),c0(isamp));
end
stheta=sqrt(1-ctheta.^2);

return
%fq2 is more accurate for theta angles near pi/2
function y=fq2(ceps,range,alph,btea,gamm,xi,eta,c0)
seps=sqrt(1-ceps.^2);
y=range.^2.*(alph.*ceps.^2-btea.*seps.*ceps+gamm)+2*range.*(xi.*ceps-eta*seps)+c0;
return

function y=fq1(cthet,range,alph,btea,gamm,xi,eta,c0)
sthet=sqrt(1-cthet.^2);
y=range.^2.*(alph.*sthet.^2+btea.*sthet.*cthet+gamm)+2*range.*(xi.*sthet+eta.*cthet)+c0;
return