function [zsurf, llh_srf]=tcs2heightT(tcs, origin, varargin);
%
% PURPOSE: Compute tcs z distance above surface (WGS84 ellipsoid)
% USAGE:  [height]=tcs2height(tcs, origin) or tcs2height(tcs, origin, error)
%         [height]=tcs2height(tcs)                => Spherical earth
% INPUT VARIABLES:
%          tcs    = xyz coordinate              [3xN]
%          origin = llh location of tcs origin  [3x1]
%          error  = surface error (Defalut=1.e-6 m)
%
% OUTPUT VARIABLES:
%          zsurf = z distance from tcs(3,:) to surface
%                = NaN if tcs(3,:) is below surface
%
%
%UNITS:    All distance measures are meters
%          All angle measures are radians
%
err=1.e-6;
[n1,n2]=size(tcs);
llh_surf=NaN*ones(n1,n2);
if nargin==1
    origin=[];
    llh_srf=[];
elseif ~isempty(varargin) 
    err=varargin{1};
end
if ~isempty(origin)
    for nvec=1:n2;
        if n2==1
            tcs_stp=tcs;
        else
            tcs_stp=tcs(:,nvec);
        end
        llh_stp=tcs2llhT(tcs_stp,origin);
        dnorm=2*err;
        while dnorm>err
            llh_stp(3)=0;
            tcs_surf=llh2tcsT(llh_stp,origin);  %Project point to surface
            tcs_stpNew=tcs_stp;
            tcs_stpNew(3)=tcs_surf(3);          %Move step to z of surface projection
            dnorm=norm(tcs_stp-tcs_stpNew);     %Stop when change is negligible
            tcs_stp=tcs_stpNew;
            llh_stp=tcs2llhT(tcs_stp,origin);
        end
        zsurf(nvec)=tcs(3,nvec)-tcs_stp(3);
        llh_srf(:,nvec)=llh_stp;
    end
else %Spherical earth
    Radius_Earth = 6378166.0;       % radius of Earth in meters
    if n2>1
        range = sqrt(tcs(1,:).^2+tcs(2,:).^2+tcs(3,:).^2);
        theta=zeros(size(range));
        iOK=find(range~=0);
        theta(iOK)=asin(tcs(3,iOK)./range(iOK));
        phi       =asin(range.*cos(theta)/Radius_Earth);
    else
        range = norm(tcs);
        if range==0
            theta=0;
        else
            theta=asin(tcs(3)./range);
            phi  =asin(range.*cos(theta)/Radius_Earth);
        end
    end
    zsurf=tcs(3,:)+Radius_Earth.*(1-cos(phi));
end
return