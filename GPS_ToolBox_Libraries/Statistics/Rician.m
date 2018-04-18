function pdf=Rician(p,b)
%
%   Rician Distribution   |sum(v+sqrt(SNR))|^2 
%
%INPUTS:
%   p    = power*Neff (vector)  
%   b  =   signal to noise power ratio*Neff 
%   Neff = Effective number of independent samples
%
pdf=zeros(size(p));
pdf=exp(-(p+b)).*besseli(0,2*sqrt(b*p));
return