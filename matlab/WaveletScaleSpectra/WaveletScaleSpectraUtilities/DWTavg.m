function     DWT_bar=DWTavg(DWT,nseg)
[nscale,nsamps]=size(DWT);
nsegs=length(nseg);
nsamp_seg=nseg(2);
DWT_bar=zeros(nscale,length(nseg));
for ns=1:nsegs
    n1=(ns-1)*nsamp_seg+1; n2=min(n1+nsamp_seg-1,nsamps);
    DWT_bar(:,ns)=DWT_bar(:,ns)+sum(DWT(:,n1:n2),2)/nsamp_seg;
end
return
