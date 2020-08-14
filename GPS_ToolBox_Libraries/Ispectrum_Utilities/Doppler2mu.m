function   [rhoOveff,mu_p,log10Imu_p,funF,nfig]=Doppler2mu(mu,Imu,fDop_P,varargin)
%Generates rhoOveff scaling to capture full Ispectrum range
%Defining variables are U,p1,p2, and mu0
%
%Written by Chuck Rino
% 
if isempty(varargin)
    DISPLAY=0;
    nfig=[];
else
    DISPLAY=1;
    plotID=varargin{1};
end

rhoOveff=mu(1)/fDop_P(1);
f1=@(rhoOveff)obj_S4(rhoOveff,mu,Imu,fDop_P);
[rhoOveff,funF]=fminsearch(f1,rhoOveff);
mu_p=rhoOveff*fDop_P;
log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip'); 

if DISPLAY==1
    %Capture range to dBrng from peak
    figure
    plot(log10(mu),dB10(Imu),'c')
    hold on
    plot(log10(mu_p),10*log10Imu_p,'r')
    grid on
    text(log10(mu(1))+1,dB10(Imu(1)), ['rhoOveff=',num2str(rhoOveff)]) 
    title(plotID)
    ylabel('I(\mu)-dB')
    xlabel('\mu')
    bold_fig
end
return