function   [psi_mean,psi_res]=LogDetrend(psi,Jmax);
%
psi=psi(:).';
logPsi=log(psi);
psi_mean=NoiseRemoval(logPsi,Jmax);
psi_mean=exp(psi_mean);
psi_res=exp(logPsi-psi_mean);
return