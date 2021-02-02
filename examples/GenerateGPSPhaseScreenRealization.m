%Generate GPS phase-screen realization
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
    N=nicefftnum(ceil(300)/dt);   %5 min
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

nP=input('1 for single power law, 2 for 2-component power law ');
if isempty(nP)
    nP=2;
else
    nP=1;
end
switch nP
    case 1
       p1=2; mu0=.5;  p2=2.5; 
    case 2
       p1=2.5; mu0=.5; p2=3.5; 
end
if mu0>=1
    Cpp=U0;
else
    Cpp=U0/mu0^(p2-p1);
end
fprintf('Ispectrum computations \n')
[Imu,mu,S4,~]=Ispectrum(U0,p1,p2,mu0);                                                     %Intensity SDF          (<=U + p1,p2,&,mu0)
[Pmu,nb]  =generatePmu(Cpp,p1,p2,mu0,mu);                                              %Phase Screen SDF  (<= Cpp +p1,p2,&mu0)
[Dmu,mu2]=Dspectrum(U0,p1,p2,mu0);                                                         %Complex signal SDF
plotID=generateIDstring(U0,p1,p2,mu0,[],S4);
plotID=strrep(plotID,'_','-');

[rhoOveff,mu_p,log10Imu_p]=Doppler2mu(mu,Imu,fDop_P,plotID);            
PhaseCorrection=input('Input 0 for no phase unwrapp error correction, else CR ');
if isempty(PhaseCorrection)
    PhaseCorrection=1;
else
    PhaseCorrection=0;
end
%%%%%%%%%%%Display Ispectrum Calculations%%%%%%%%%%%%%
nfig1=figure;
plot(log10(mu),dB10(Imu),'r')
hold on
plot(log10(mu),dB10(Pmu),'b')
grid on
xlabel('log10(mu)')
ylabel('\Phi(\mu)-dB')
title(plotID)
hold on
plot(log10(mu),dB10(Dmu),'g')
hold on
plot(log10(mu_p),10*log10Imu_p,'m')
legend('Ispec','Pspec','Dspectrum','mu-range')
bold_fig

nsamps=length(log10Imu_p);
%Finer sampling for phase-screen model
nfft=8192;
nfft=max(nicefftnum(2*nsamps+1),nfft);

%Recompute fDop and fDop_P <= nfft
dfDop=1/(nfft*dt);
nfft_span=[-nfft/2:nfft/2-1];
fDop=nfft_span*dfDop;
t_sec=(0:nfft-1)*dt;

fDop_P=fDop(nfft/2+2:nfft); dfDop=diff(fDop(1:2));
nfreq_fit=1;  
nfreqs=3;

fprintf('Phase screen realizations \n')
SEED=rng;  %One white noise realization for all frequencies

U_fit=U0; mu0_fit=mu0;rhoOveff_fit=rhoOveff;
%generate realization at nfreq_fig
nfreq_ext=setdiff([1,2,3],nfreq_fit);  %Change if nfreqs<3
Phase=zeros(nfreqs,nfft);
S4=zeros(1,nfreqs);
[psi,Phase(nfreq_fit,:),fracMom]=generateSurrogateHk(U_fit,p1,p2,mu0_fit,rhoOveff_fit,dt,nfft,SEED,1);
S4(nfreq_fit)=sqrt(fracMom(2)-1);
phase=unwrapPhase(atan2(imag(psi),real(psi)));
if PhaseCorrection==1
    [phase_int,nint]=ConstructPhase(psi);
    error=phase_int(1:nint:end)-phase;
    nfitP=figure;
    subplot(2,1,1)
    plot(t_sec,phase,'b')
    hold on
    plot(t_sec,phase_int(1:nint:end),'r')
    grid on
    ylabel('phase-rad')
    legend('unwraped','corrected','location','northwest')
    title([plotID,' nint=',num2str(nint)])
    axis([0, max(t_sec), floor(min(phase)), ceil(max(phase))])
    subplot(2,1,2)
    plot(t_sec,error,'r')
    grid on
    axis([0, max(t_sec), min(-1,min(error)), max(1,max(error))])
    ylabel('error-rad')
    xlabel('t-sec')
    bold_fig
end
Isim(nfreq_fit,:)=abs(psi).^2;
if PhaseCorrection==1
    Psim(nfreq_fit,:)=phase_int(1:nint:end); 
else
    Psim(nfreq_fit,:)=phase; 
end
Psim(nfreq_fit,:)=Psim(nfreq_fit,:)-linex(1:nfft,1,nfft,Psim(nfreq_fit,1),Psim(nfreq_fit,nfft));  %Remove linear trend
for nfreq=nfreq_ext
    %[X0]=freqExtrapolate(U_fit,p1,mu0_fit,p2,rhoOveff_fit,freqGPS(nfreq_fit),freqGPS(nfreq));
    %U=X0(1); mu0=X0(3); rhoOveff=X0(5);
    %[psi,Phase(nfreq,:),fracMom]=generateSurrogateHk(U,p1,p2,mu0,rhoOveff,dt,nfft,SEED);
    
    fscale=freqGPS(nfreq_fit)/freqGPS(nfreq);
    [psiT,Phase,fracMom]=generateSurrogateHk(U_fit,p1,p2,mu0,rhoOveff,dt,nfft,SEED,fscale);  
    S4(nfreq)=sqrt(fracMom(2)-1);
    Isim(nfreq,:)=abs(psi).^2;
    Psim(nfreq,:)=unwrap(atan2(imag(psi),real(psi)));
    if PhaseCorrection~=1
       [phase_int,nint]=ConstructPhase(psi);
       Psim(nfreq,:)=phase_int(1:nint:end);
    else
        Psim(nfreq,:)=phase;
    end
    Psim(nfreq,:)=Psim(nfreq,:)-linex(1:nfft,1,nfft,Psim(nfreq,1),Psim(nfreq,nfft));  %Remove linear trend
end

t_sec=(0:nfft-1)*dt;
IDstring=generateIDstring(U_fit,p1,p2,mu0,rhoOveff,[]);
fileID=IDstring;

nfig3=figure;
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