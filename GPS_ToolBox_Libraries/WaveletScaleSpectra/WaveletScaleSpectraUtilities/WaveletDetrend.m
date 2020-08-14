function [y,d,Wavelet_Struct,yp]=WaveletDetrend(f,varargin)
%
ndata=length(f);
datap=[f(:)',fliplr(f(:)')];   
d=nextpow2(2*ndata); n=2^d;
if n>2*ndata
    yp=[datap,datap(2*ndata)*ones(1,n-2*ndata)];
end
J=d;
if isempty(varargin)
    Jmin=8;
else
    Jmin=varargin{1};
end
Family='Symmlet';
switch Family
    case 'Haar'
        par=[];
    case 'Daubechies'
        par=4;      %4,6,8,...20
    case 'Coiflet'  
        par=1;      %1,2,3,4,5 
    case 'Symmlet'
        par=10;      %4,5,6,...10
    otherwise
        error('Not Implemented');
end
Type=Family;
qmf = MakeONFilter(Type,par);
wc = FWT_PO(yp,1,qmf);
Wavelet_Struct=struct('d',d,'Type',Type,'qmf',qmf,'wc',wc);
for j=0:Jmin
    wc(dyad(j))=0;
end
y=IWT_PO(wc,1,qmf);
y=y(1:ndata);
return