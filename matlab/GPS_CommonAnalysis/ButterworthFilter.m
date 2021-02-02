function  [data_bar]=ButterworthFilter(data,dt,f3dB)
%ButterworthFilter
%
nsamps=length(data);
nfft=nicefftnum(nsamps);
df=1/dt/nfft;
fDop=(-nfft/2:nfft/2-1)*df;
BF=sqrt(1./(1+(fDop/f3dB).^12));
data_hat=fft(data,nfft).*fftshift(BF);
data_bar=ifft(data_hat);
data_bar=data_bar(1:nsamps);

DEBUG=0;
if DEBUG==1
    nfig=figure;
    plot(log10(fDop(nfft/2+2:nfft)),dB10(abs(data_hat(2:nfft/2).^2)),'r')
    grid on
    hold on
    plot(log10(fDop(nfft/2+2:nfft)),dB10(abs(BF(nfft/2+2:nfft).^2)),'b')
    OK=input('CR to continue ');
    close(nfig);
end
return

