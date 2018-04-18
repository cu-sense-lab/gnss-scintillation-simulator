function [y,d,Wavelet_Struct,yp]=WaveletFilter(data,varargin)
%
%Filter data by selectively removing DWT coefficients from smallest scale
%j=d toward largest scale j=0   
%Within the filter routine d is next power of 2 greater than twice the data
%extent to accomodate folded periodic extension
%

ndata=length(data);
%Folded periodic extension
datap=[data(:)',fliplr(data(:)')];   
d=nextpow2(2*ndata); n=2^d;
if n>2*ndata
    yp=[datap,datap(2*ndata)*ones(1,n-2*ndata)];
else
    yp=datap;
end
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
%Discrete wavelet transform
wc = FWT_PO(yp,1,qmf);

Wavelet_Struct=struct('d',d,'Type',Type,'qmf',qmf,'wc',wc);

if Jmin<=d-1
    %Denoise => remove coefficients d, d-1, ..., Jmin
    for j=d-1:-1:Jmin
        wc(dyad(j))=0;
    end
end
%Extract ndata subset:
y=IWT_PO(wc,1,qmf);
y=y(1:ndata);
return