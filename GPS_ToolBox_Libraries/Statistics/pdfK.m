function pdf=pdfK(I,alph)
%
%INPUTS:
%   p    = power (vector)
% 
pdf=(2/gamma(alph))*(alph*I).^(alph/2-1/2).*besselk(2*sqrt(alph*I),alph-1);
return