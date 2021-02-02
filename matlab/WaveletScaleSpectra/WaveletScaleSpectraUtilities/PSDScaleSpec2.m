function   [xxScale,yyScale]=PSDScaleSpec2(nsamp,nfft,nsegs,dx,S_data,J,jmax,varargin)
%USAGE:    SPEC=PSDScaleSpec(data,nsamp,nfft,nsegs,dx,S_data,J,jmax))
%
%PURPOSE  Compute Scale Spectra
%
%INPUT:
%          data  =data array
%          nsamp =length(data)
%          nfft  =length of segments
%          nseg  =number of segments
%          S_data,J,jmax outputs from ComputeScaleSpectrum
%
%OUTPUT:
%               xxScale      = log10(DWT scales)
%               yyScale     = 10log10(ScaleSpec)
%
sj=(2*dx)*2.^(J-jmax-1:-1:0);
qScale=1./sj;                 

for nseg=1:nsegs
    n1=(nseg-1)*nfft+1; n2=n1+nfft-1;
    xxScale=log10(qScale);
    yyScale=dB10(S_data(jmax:J-1,nseg))';
end %End of segment
return