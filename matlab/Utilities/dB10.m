function IdB=dB10(I)
%    Compute 10*log10(I)  
%USAGE:  IdB=dB10(I)
%INPUTS:
%    I = input Intensity vector  
%OUTPUTS:
%    dB if I>0 else NaN
%

% Chuck Rino
% Rino Consulting
% July 2010
%
IdB=NaN*ones(size(I));
iOK=find(abs(I)>0);
IdB(iOK)=10*log10(abs(I(iOK)));
return