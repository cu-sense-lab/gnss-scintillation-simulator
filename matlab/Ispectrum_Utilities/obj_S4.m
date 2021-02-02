function F=obj_S4(rhoF,mu,Imu,qP)
%object function for rhoF calculation
rhoF=max(rhoF,mu(1)/qP(1));
mu_p=rhoF*qP;
log10Imu_p=interp1(log10(mu),log10(Imu),log10(mu_p),'pchip');  
Imu_p=10.^log10Imu_p;
F=-sum(Imu_p)*mu_p(1)/pi;
return