function [PDF] = NakagamiPDF(I,SI)
% Generate Nakagami Probability density function
% INPUTS:
% SI   = Modulation Index (scalar)
% I    = Intensity normalized to mean  (vector)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SI==0
    error('SI>0')
else
    m=1/SI;
end
PDF=m^m/gamma(m).*I.^(m-1).*exp(-m*I);
return