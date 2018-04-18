function [SI,Ibar,nstep,fracMom]=computeSI(I,nspan,varargin)
%USAGE:  [SI,Ibar]=computeSI(I,nspan,nstep)
%
%Symmetric boxcar average with reflected end corrections
%
if isempty(varargin)
    nstep=floor(nspan/2);
else
    nstep=varargin{1};
end
npts=length(I);
nmid=floor(nspan/2)+1; nupd=(nmid-1);
Ibar=zeros(size(I(1:nstep:npts))); SI=Ibar;
fracMom=zeros(5,length(SI));
nn=0;
for m=1:nstep:npts
    %mspan=[max(1,m-nupd):min(m+nupd,npts)]; ntot=length(mspan);
    mspan=m:min(m+nspan-1,npts); ntot=length(mspan);
    nn=nn+1;
    Ibar(nn)=sum(I(mspan))/ntot;
    SI(nn)=sqrt(sum(I(mspan).^2)/(ntot*Ibar(nn)^2)-1);
    fracMom(:,nn)=generate_fracMom(I(mspan)) ;
end
return