function data_int=fftInterp(data,nint)
%t_int=(0:nsamps-1)*dt  => (0:nint*nsamps-1)*dt/nint
nsamps=length(data);
nfft=nicefftnum(nsamps);
data_hat=fft(data,nfft);     %zero pad data to nfft>nsamps
data_int=zeros(1,nint*nfft);
data_int(1:nfft/2)=data_hat(1:nfft/2);
data_int((nint-1/2)*nfft+1:nint*nfft)=data_hat(nfft/2+1:nfft);
data_int=nint*ifft(data_int);
data_int=data_int(1:nint*nsamps);
return