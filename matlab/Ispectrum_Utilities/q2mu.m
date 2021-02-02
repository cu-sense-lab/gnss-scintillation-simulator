function   [rhoF,mu_p,log10Imu_p,funF,nfig]=q2mu(mu,Imu,q_P,varargin)
%Generates rhoF scaling to capture full Ispectrum range
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
%mu range that defines the intensity SDF and S4 
rhoF0=mu(1)/q_P(1);
f1=@(rhoF)obj_S4(rhoF,mu,Imu,q_P);
[rhoF,funF,~,~]=fminsearch(f1,rhoF0);
mu_p=rhoF*q_P;
log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip'); 
if DISPLAY==1
    nfig=figure;
    plot(log10(mu),dB10(Imu),'c')
    hold on
    plot(log10(mu_p),10*log10Imu_p,'r')
    grid on
    text(log10(mu(1))+1,dB10(Imu(1)), ['rhoF=',num2str(rhoF)])
    title(plotID)
    ylabel('I(\mu)-dB')
    xlabel('\mu')
    bold_fig
end
return