function [PDF] = LogNormalPDF(I,Ibar,sig2)
% Generate Log Normal Probability Density Function
% INPUTS:
% SI   = Scintillation Index (scalar)
% I    = intensity normalized to mean (vector)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDF=1/sqrt(2*pi*sig2)*exp(-((log(I)+Ibar).^2)/(2*sig2))./I;
return