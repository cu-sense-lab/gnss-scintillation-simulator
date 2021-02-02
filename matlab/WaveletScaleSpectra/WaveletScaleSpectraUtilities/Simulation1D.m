%Simulation 1D
clear
close all

Video=1;
if Video==1
    videoFile='Hodo.avi';
    vidObj = VideoWriter(videoFile);
    open(vidObj);
end

%%%%%Fixed Parameters%%%%%%%%%%%%%%%%%%%%%
c = 299792458;            %Speed of light (vacuum)
re=2.819740289e-15;   %Classical electron radius
K=re*c/(2*pi)*1.e16;     %TEC conversion factor
dtr=pi/180;

%GPS frequencies
freqGPS(1)=154*10.23e6;  freqID{1}='L1';
freqGPS(2)=120*10.23e6;  freqID{2}='L2';
freqGPS(3)=115*10.23e6;  freqID{3}='L5';

nfft=16384;  %2^14
dt=0.01;
veff=100000/nfft/dt;
t_sec=(0:nfft-1)*dt;
y=veff*t_sec-50000;
rhoOveff=10;
k=2*pi*freqGPS(1)/c;
load('rngSEED');
%SEED=rng('shuffle');
%save('rngSEED','SEED');
p=3.5;
p1=p; p2=p; mu0=1;
IDp1=num2str(fix(p*100)/100);
IDrhoOveff=num2str(fix(rhoOveff*100)/100);
IDdt=num2str(fix((1/dt)*100)/100);
ID=['GPS Model p=',IDp1,'  rhoOveff=',IDrhoOveff,'  Sample Rate=',IDdt, ' Hz'];
nSteps=200;
U=linspace(1,30,nSteps);
fracMom=zeros(5,nSteps);
psi=zeros(nfft,nSteps)+1i*zeros(nfft,nSteps);
nfig=figure;
Amax=5;
for n=1:nSteps
    [psi(:,n),phase0,fracMom(:,n)]=generateSurrogateHk(U(n),p1,p2,mu0,rhoOveff,dt,nfft,SEED);
    I=abs(psi(:,n).^2);
    [PDF,~,yI]=Distribution (I,50);
    S4=sqrt(fracMom(2,n)-1);
    alphamu=getAlphaMu(S4);
    [S4,PDFam]=AlphaMu(yI,alphamu(1),alphamu(2));
    [psi_mean,psi_res]=LogDetrend(psi(:,n),5);
    
    figure(nfig)
    subplot(4,1,1)
    hold off
    plot(y/1000,dB10(I),'r')
    axis([-50 50 -40 10])
    ylabel('I-dB')
    xlabel('y-km')
    title(['S4=',num2str(S4) '      U=',num2str(U(n))])
    grid on
    subplot(4,1,2)
    hold off
    plot(dB10(yI),PDF/diff(yI(1:2)),'b')
    grid on
    hold on
    plot(dB10(yI),PDFam,'r')
    hold on
    plot(dB10(yI),exp(-yI),'g')
    axis([-40 10 0 2])
    ylabel('PDF')
    xlabel('I-dB')
    subplot(4,1,[3,4])
    hold off
    plot(real((psi_res)),imag((psi_res)),'r.')
    hold on
    plot(real((psi_mean)),imag((psi_mean)),'g.')
    grid on
    axis([-Amax Amax -Amax Amax])
    xlabel('Real')
    ylabel('Imaginary')
    axis square
    legend('<log(psi)>','\Delta log(psi)','location','SouthEastOutside')
    bold_fig
    drawnow
    if Video==1
        currFrame = getframe(gcf);
        for nrep=1:10
            writeVideo(vidObj,currFrame);
        end
    end
end
if Video==1
    close(vidObj)
end

