%PhaseScreenModelDemo for spatial realization
%Written by Charles Rino
%crino@comcast.net
%

clear
close all

N=2^16;
Y=200000;
dy=Y/N;
dq=2*pi/Y;
q=(-N/2:N/2-1)*dq;
U=1;
SAVE_FIGS=input('Input 1 to save figures ');
if ~isempty(SAVE_FIGS)
    figs_dir=[pwd,'\JPGs'];
    if ~exist('figs_dir',dir)
        mkdir(figs_dir)
    end
end
DateStamp=generateDateStamp;

nP=input('1 for single power law, 2 for 2-component power law ');
if isempty(nP)
    nP=2;
else
    nP=1;
end
switch nP
    case 1
       p1=2; mu0=.5;  p2=2; 
    case 2
       p1=2.5; mu0=.5; p2=3.5; 
end

if mu0>=1
    Cpp=U;
else
    Cpp=U/mu0^(p2-p1);
end
fprintf('Ispectrum computations \n')
[Imu,mu,S4]=Ispectrum(U,p1,p2,mu0);                                      %Intensity SDF          (<=U + p1,p2,&,mu0)
[Pmu,nb]  =generatePmu(Cpp,p1,p2,mu0,mu);                           %Phase Screen SDF  (<= Cpp +p1,p2,&mu0)
[Dmu]=Dspectrum(U,p1,p2,mu0,mu);                                         %Complex signal SDF
%[Smu]=structure_function(U,p1,p2,mu0,mu);                             %StructureFunction
qP=q(N/2+2:N);
[rhoF,mu_p,log10Imu_p,~,nfig1]=q2mu(mu,Imu,qP,1);
figname=['\Fig1',DateStamp];
if ~isempty(SAVE_FIGS)
    saveas(nfig1,[figs_dir,figname],'jpg');
end

plotID=generateIDstring(U,p1,p2,mu0,[],S4);
plotID=strrep(plotID,'_','-');
SEED=rng;
nfft=N;
[psi,phase0,fracMom,~,~,~]=...
                      generateSurrogateHk_y(U,p1,p2,mu0,rhoF,dy,nfft,SEED);

Isim=abs(psi).^2;
SDF_I=abs(fft(Isim,nfft)).^2/nfft;
SDF_IP=SDF_I(2:nfft/2);
                 
SDF_P=abs(fft(psi,nfft)).^2/nfft;
SDF_PP=SDF_P(2:nfft/2);

y=(0:nfft-1)*dy;
nfig2=figure;
plot(y/1000,dB10(Isim),'r')
grid on
text(10,15, ['S4=',num2str(sqrt(fracMom(2)-1))])
xlabel('y-km')
ylabel('I-dB')
title(plotID)
axis([0 max(y/1000) -30 20])
bold_fig

figname=['\Fig2',DateStamp];
if ~isempty(SAVE_FIGS)
    saveas(nfig2,[figs_dir,figname],'jpg');
end

CF=rhoF;   
nfig3=figure;
subplot(2,1,1)
plot(log10(qP),dB10(SDF_IP),'r')
grid on
hold on
plot(log10(mu_p/rhoF),10*log10Imu_p+dB10(CF),'g')
axis([-5 0 -60 40])
title(plotID)
ylabel('SDF-I dB')

subplot(2,1,2)
plot(log10(qP),dB10(SDF_PP),'r')
grid on
hold on
plot(log10(mu/rhoF),dB10(Dmu*CF),'g')
axis([-5 0  -60 40])
xlabel('log10(q)')
ylabel('SDF-P dB')
bold_fig

figname=['\Fig3',DateStamp];
if ~isempty(SAVE_FIGS)
    saveas(nfig3,[figs_dir,figname],'jpg');
end