function [xxScale,yyScale,sj]=WaveletAnalysis2(data,dy,jmin)
%Compute wavelet scale spectrum
Family='Symmlet';
data=data(:)';
%Compute discrete wavelet transform
[~,wc,J]=ComputeDWT(data,Family);
%Compute scale spectrum
[S_PLP]=ComputeScaleSpectrum(J,wc,2);
S_PLP=S_PLP/sqrt(2);   
sj=(2*dy)*2.^(J-1-jmin:-1:0);
qScale=1./sj;
nseg=1;
%Scale spectrum interpolated to xxPSD
xxScale=log10(qScale);
yyScale=dB10(S_PLP(jmin:J-1,nseg))';
return