function   [rhoOveff,S4]=InitialScale(U,p1,p2,mu0,PSD_Pfit,fDop_Pfit)
%
[Imu,mu,S4,~,~,~]=Ispectrum(U,p1,p2,mu0);           %Generate Ispectrum
[~,nD_max]=max(PSD_Pfit);
[~,nI_max]=max(Imu);
rhoOveff=mu(nI_max)/(2*pi*fDop_Pfit(nD_max));
return