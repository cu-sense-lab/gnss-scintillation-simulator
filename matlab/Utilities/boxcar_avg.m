function   [vbar,nstep]=boxcar_avg(v,nspan,varargin)
%
if isempty(varargin)
    nstep=floor(nspan/2);
else
    nstep=varargin{1};
end

npts=length(v);
nmid=floor(nspan/2)+1; nupd=(nmid-1);
vbar=zeros(size(v(1:nstep:npts)));
nn=0;
for m=1:nstep:npts
    mspan=[max(1,m-nupd):min(m+nupd,npts)]; ntot=length(mspan);
    nn=nn+1;
    vbar(nn)=sum(v(mspan))/ntot;
end
return