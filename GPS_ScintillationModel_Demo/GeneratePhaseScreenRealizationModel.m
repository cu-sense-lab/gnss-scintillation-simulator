%Generate phase-screen realization using GPS_ScintillationModel parameters
%Written by Charles Rino
%crino@comcast.net
%

clear
close all

%%%%%Fixed Parameters%%%%%%%%%%%%%%%%%%%%%
c = 299792458;            %Speed of light (vacuum)
re=2.819740289e-15;   %Classical electron radius
K=re*c/(2*pi)*1.e16;     %TEC conversion factor
dtr=pi/180;

%GPS frequencies
freqGPS(1)=154*10.23e6;  freqID{1}='L1';
freqGPS(2)=120*10.23e6;  freqID{2}='L2';
freqGPS(3)=115*10.23e6;  freqID{3}='L5';

TN=input('Input dt and # of samples (N) as [dt, N] ');
if isempty(TN)
    dt=0.01;
    N=nicefftnum(600/dt);
else
    dt=TN(1);
    T=TN(2);
end
    
fDop=(-N/2:N/2-1)/(N*dt);  %Full Doppler range
fDop_P=fDop(N/2+2:N);    %Positive Doppler range excluding DC and Nyquist

U0=input('Input universal strength parameter U (0 to > 100) '); 
if isempty(U0);
    U0=1;
end
%Generate representative p1, mu0, and p2 values & model S4
%Theoretical values for Imu and mu are mimal selected for interpolation
[S4,p1,mu0,p2,Imu,mu]=GPS_ScintillationModel(U0);

%Determine the mu range that defines the intensity SDF and S4 
dBrng=50;
dBImu=dB10(Imu);
[dBImu_max,nmax]=max(dBImu);
nOK=find(abs(dBImu-dBImu_max)<dBrng);
mu=mu(nOK); Imu=Imu(nOK);

%Adust mu range to capture as much of the significant range as possible
rhoOveff=mu(1)/(2*pi*fDop_P(1));
mu_p=rhoOveff*2*pi*fDop_P;
logDiff=1;
while logDiff>0
      rhoOveff=2*rhoOveff;
      mu_p=rhoOveff*2*pi*fDop_P;
      log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip');        %Interpolate Imu to mu_p
      logDiff=log10Imu_p(end)-log10Imu_p(1);
end
Imu_p=10.^log10Imu_p;  

%Plot significant model range and interpolated range
figure
plot(log10(mu),dB10(Imu),'b')
hold on
plot(log10(mu_p),dB10(Imu_p),'r')
hold on
grid on
xlabel('log10(mu)')
ylabel('\Phi(\mu)-dB')
title(generateIDstring(U0,p1,p2,mu0,[],S4))
bold_fig

nsamps=length(Imu_p);
%Finer sampling for phase-screen model
nfft=8192;
nfft=max(nicefftnum(2*nsamps+1),nfft);

%Recompute fDop and fDop_P <= nfft
dfDop=1/(nfft*dt);
nfft_span=[-nfft/2:nfft/2-1];
fDop=nfft_span*dfDop;

fDop_P=fDop(nfft/2+2:nfft); dfDop=diff(fDop(1:2));
nfreq_fit=1;
nfreq=nfreq_fit;    
nfreqs=3;

SEED=rng;  %One white noise realization for all frequencies

U_fit=U0; mu0_fit=mu0;rhoOveff_fit=rhoOveff;
%generate realization at nfreq_fig
nfreq_ext=setdiff([1,2,3],nfreq_fit);  %Change if nfreqs<3
Phase=zeros(nfreqs,nfft);
[psi,Phase(nfreq_fit,:),fracMom]=generateSurrogateHk(U_fit,p1,p2,mu0_fit,rhoOveff_fit,dt,nfft,SEED);

S4(nfreq_fit)=sqrt(fracMom(2)-1);
Isim(nfreq_fit,:)=abs(psi).^2;
Psim(nfreq_fit,:)=unwrap(atan2(imag(psi),real(psi)));
for nfreq=nfreq_ext
    [X0]=freqExtrapolate(U_fit,p1,mu0_fit,p2,rhoOveff_fit,freqGPS(nfreq_fit),freqGPS(nfreq));
    U=X0(1); mu0=X0(3); rhoOveff=X0(5);
    [psi,Phase(nfreq,:),fracMom]=generateSurrogateHk(U,p1,p2,mu0,rhoOveff,dt,nfft,SEED);
    S4(nfreq)=sqrt(fracMom(2)-1);
    Isim(nfreq,:)=abs(psi).^2;
    Psim(nfreq,:)=unwrap(atan2(imag(psi),real(psi)));
    Psim(nfreq,:)=Psim(nfreq,:)-linex(1:nfft,1,nfft,Psim(nfreq,1),Psim(nfreq,nfft));
end

t_sec=(0:nfft-1)*dt;
IDstring=generateIDstring(U_fit,p1,p2,mu0,rhoOveff,[]);
fileID=IDstring;

figure
dBmin=min(dB10(Isim(:)));
for nfreq=1:nfreqs
    %intensity_d(nfreq,:)=WaveletFilter(intensity_d(nfreq,:),14);
    subplot(nfreqs,2,2*nfreq-1)
    plot(t_sec/60,dB10(Isim(nfreq,:)),'r')
    text(1,dBmin,['S4= ',num2str(S4(nfreq))])
    grid on
    ylabel(freqID{nfreq})
    axis([0 ceil(max(t_sec/60)) floor((dBmin-5)/10)*10,10])
    if nfreq==1
        title(['                                                      ',fileID])
        legend('I-dB')
    end
end
xlabel('                                                              UT-min')


Pmax=max(Psim(:)); Pmin=min(Psim(:));
for nfreq=1:nfreqs
    subplot(nfreqs,2,2*nfreq)
    plot(t_sec/60,Psim(nfreq,:),'r')
    grid on
    axis([0 ceil(t_sec(end)/60) floor(Pmin) ceil(Pmax)])
    if nfreq==1
        legend('\phi-rad')
    end
end
bold_fig

Display_SpectraModel(Isim,Psim,fDop,nfft,nfreq_fit,U_fit,p1,mu0_fit,p2,rhoOveff_fit,fileID)