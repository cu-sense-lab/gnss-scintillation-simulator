
function [S,nsegs,nscls]=ComputeScaleSpectrum(J,wc,jmax)
%USAGE:  [S,nsegs,nscls]=ComputeScaleSpectrum(J,wc,jmax)
%
%PURPOSE: Segmented computation of scale spectrum
%         Segments are subdivisions of the independent DWT samples versus x
%         Dyadic downsampling determines the independent samples
%
%INPUT:   J and wc from DWTSpectrum
%              jmax determines the number of segments nseg=^2^(jmax-1);
%              2^(J-j) segments, covering the wavelet support at scale 
%              
%                               
%OUTPUT:  S=array of scale spectra size(S)=[nscls,nsegs]
%                NaN=>scale spectra beyond jmin
%
%
%Speclab 850 Utilities: dyad 
%
%Written by Charles Rino
%Rino Consulting
%January 2012
%crino@comcast.com
%http://chuckrino.com/wordpress
%
nsegs=2^(jmax-1);
nscls=J-jmax;
S=NaN(nscls,nsegs);
for j=(J-1):-1:jmax
    V=abs(wc(dyad(j))).^2;
    nsum=length(V)/nsegs;
    for nseg=1:nsegs
        n0=(nseg-1)*nsum+1; n1=nseg*nsum;
        S(j,nseg)=sum(V(n0:n1))/2^j;
    end
end
return