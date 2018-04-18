function Display_SpectraModel(Isim,Psim,fDop,nfft,nfreq_fit,U_fit,p1,mu0_fit,p2,rhoOveff_fit,fileID)

%frequencies
freqGPS(1)=154*10.23e6;  freqID{1}='L1';
freqGPS(2)=120*10.23e6;  freqID{2}='L2';
freqGPS(3)=115*10.23e6;  freqID{3}='L5';

fDop_P=fDop(nfft/2+2:nfft); dfDop=diff(fDop(1:2));
PSD_P=zeros(3,nfft/2-1);
PSDphi_P=zeros(3,nfft/2-1);
for nfreq=1:3
    PSD=abs(fft(Isim(nfreq,:),nfft)).^2/nfft;
    PSD_P(nfreq,:)=PSD(2:nfft/2);
    PSD_P(nfreq,:)=PSD_P(nfreq,:)/(nfft*dfDop);                          
    PSD=abs(fft(Psim(nfreq,:),nfft)).^2/nfft;
    PSDphi_P(nfreq,:)=PSD(2:nfft/2);
    PSDphi_P(nfreq,:)=PSDphi_P(nfreq,:)/(nfft*dfDop);          
end

dBmaxData=ceil(max(dB10(PSD_P(3,:)))/10)*10;
logfDmin=floor(min(log10(fDop_P(:))))+0.5;
logfDmax=ceil(max(log10(fDop_P(:))))-0.5;

figure
subplot(2,1,1)
[X0]=freqExtrapolate(U_fit,p1,mu0_fit,p2,rhoOveff_fit,freqGPS(nfreq_fit),freqGPS(1));
U=X0(1); mu0=X0(3); rhoOveff=X0(5);
plot(log10(fDop_P),dB10(PSDphi_P(1,:)/rhoOveff),'r')
hold on
plot(log10(fDop_P),dB10(PSD_P(1,:)/rhoOveff),'b')

[Imu,mu,S4,Cpp]=Ispectrum(U,p1,p2,mu0);
mu_p=2*pi*fDop_P*rhoOveff;                                                        %fDop_P to mu_p
log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip');        %Interpolate Imu to mu_p
Imu_p=10.^log10Imu_p;
hold on
plot(log10(fDop_P),dB10(Imu_p),'m')
Pmu_p=generatePmu(Cpp,p1,p2,mu0,mu_p);
hold on
plot(log10(fDop_P),dB10(Pmu_p),'m')
dBmax=max([dBmaxData,dB10(Pmu_p)]);
grid on
title(fileID)
ylabel('L1 PSD-dB')
text(logfDmin,dBmaxData-50,['S4=',num2str(fix(S4*100)/100)],'col','b')
legend('\phi','Int','theory')
subplot(2,1,2)
[X0]=freqExtrapolate(U_fit,p1,mu0_fit,p2,rhoOveff_fit,freqGPS(nfreq_fit),freqGPS(2));
U=X0(1); mu0=X0(3); rhoOveff=X0(5);
plot(log10(fDop_P),dB10(PSDphi_P(2,:)/rhoOveff),'r')
hold on
plot(log10(fDop_P),dB10(PSD_P(1,:)/rhoOveff),'b')

[Imu,mu,S4,Cpp]=Ispectrum(U,p1,p2,mu0);
mu_p=2*pi*fDop_P*rhoOveff;                
plot(log10(fDop_P),dB10(PSD_P(2,:)/rhoOveff),'b')
mu_p=2*pi*fDop_P*rhoOveff;                                                          %fDop_P to mu_p
log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip');        %Interpolate Imu to mu_p
Imu_p=10.^log10Imu_p;
hold on
plot(log10(fDop_P),dB10(Imu_p),'m')
Pmu_p=generatePmu(Cpp,p1,p2,mu0,mu_p);
hold on
plot(log10(fDop_P),dB10(Pmu_p),'m')
dBmax=max([dBmaxData,dB10(Pmu_p)]);
grid on
ylabel('L2 PSD-dB')
text(logfDmin,dBmaxData-50,['S4=',num2str(fix(S4*100)/100)],'col','b')
xlabel('log10(fDop)-Hz')
bold_fig