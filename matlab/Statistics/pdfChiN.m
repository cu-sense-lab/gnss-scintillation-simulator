function pdf=pdfChiN(p,sigc,N)
%
%   Chi disribution for average of N pulses
%
%INPUTS:
%   p    = power (vector)
%   sigc = <p> 
%   N    = Number of pulses
%
c=N/sigc/gamma(N);
P=N*p/sigc;
pdf=c*P.^(N-1).*exp(-P);
return