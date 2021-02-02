function IspecParams=generateIspecParams(U,p1,p2,mub)
%USAGE:  IspecParams=generateIspecParams(U,p1,p2,mub)
%
%Generates character parameter input for call to Spectrumintegration.exe
%
IspecParams=[num2str(U) ,' ', num2str(p1), ' ',num2str(p2), ' ',num2str(mub), ' 0.0 0.0'];
end