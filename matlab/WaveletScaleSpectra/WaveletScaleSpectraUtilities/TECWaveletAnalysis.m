function    [DWT,DWT_bar,fScale,J]=TECWaveletAnalysis(data,t_sec,nseg)

nsamps=length(t_sec);
dt=diff(t_sec(1:2));
Family='Symmlet';
%Compute discrete wavelet transform
[DWT,~,J]=ComputeDWT(data,Family);      %J=nextpow2(length(data))+1  <= folded extended data set
DWT=DWT(:,1:nsamps)/dt;                        %Truncate DWT to data support
%sj(1)=2*dt <= smallest scale that supports wavelet
%sj(J-2)      <= largest scale fully within data
sj=(2*dt)*2.^(J-2:-1:0);
fScale=1./sj;
DWT_bar=DWTavg(DWT,nseg);
return
