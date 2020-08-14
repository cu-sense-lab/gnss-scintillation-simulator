function  [data_d,data_bar]=DetrendIntensity(data,t_sec,method,param)
%Intensity detrend utility
%
if method==1
    f3dB=param;
elseif method==2
    Jmax=param;
else
    error('input parameter error')
end
[nfreqs,nsamps]=size(data);
nfft=nicefftnum(nsamps);
dt=diff(t_sec(1:2)); df=1/dt/nfft;
fDop=(-nfft/2:nfft/2-1)*df;

data_d=zeros(size(data));
data_bar=zeros(size(data));
if method==1
   BF=sqrt(1./(1+(fDop/f3dB).^12));
    for nfreq=1:nfreqs
        Data=data(nfreq,:);
        data_hat=fft(Data,nfft);
        data_barT=real(ifft(data_hat.*fftshift(BF)));
        data_bar(nfreq,:)=data_barT(1:nsamps);
        data_d(nfreq,:)=Data./data_bar(nfreq,:);
    end
else
    d=nextpow2(nsamps);
    Jmax=min(Jmax,d);
    for nfreq=1:nfreqs
        Data=data(nfreq,:);
        y=NoiseRemoval(Data,Jmax);
        data_bar(nfreq,:)=y(1:nsamps);
        data_d(nfreq,:)=data(nfreq,:)./y;
    end
end
return