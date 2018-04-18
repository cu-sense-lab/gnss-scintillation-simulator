function [DWT_Spec,wc,J,yp,qmf]=ComputeDWT(data,varargin)
%USAGE:  [DWT_Spec,wc,J,yp,qmf]=ComputeDWT(data,varargin)   
%
%PURPOSE: Compute Discret Wavelet Transform using Speclab 850 utilities
%
%INPUT:
%       data    = data array
%       varargin= 'wavelet' family default 'Daubechies'
%OUTPUT:
%       DWT_Spec= discrete wavelet transformation (Jx2^J)
%       wc      = dyadic array of wavelet coefficients 
%       J       = dyadic dimension extended to nearest power of 2
%       yp      = data extension to neasest power of 2 (folded periodic)  
%       qmf     = quadrature mirror filter 
%
%Speclab 850 Utilities: MakeONFilter, FWT_PO, dyadlength, dyad 
%
%Written by Charles Rino
%Rino Consulting
%January 2012
%crino@comcast.com
%http://chuckrino.com/wordpress
%
if ~isempty(varargin)
    Family=varargin{1};
else
    Family='Daubechies';
end
DEBUG=0;
ndata=length(data);
datap=[data(:)',fliplr(data(:)')];   
d=nextpow2(2*ndata); n=2^d;
if n>2*ndata
    yp=[datap,datap(2*ndata)*ones(1,n-2*ndata)];
else
    yp=datap;
end

%d=nextpow2(length(data)); n=2^d;
%datap=[data,fliplr(data)];      %folded periodic extension of data
%yp=datap(1:n);

%Discrete Wavelet Transform
%*****************Wavelet****************************
switch Family
    case 'Haar'
        par=[];
    case 'Daubechies'
        par=4;      %4,6,8,...20
    case 'Coiflet'  
        par=1;      %1,2,3,4,5 
    case 'Symmlet'
        par=4;      %4,5,6,...10
    otherwise
        error('Not Implemented');
end
Type=Family;
if DEBUG==1
    wave=MakeWavelet(0,0,Family,par,'Mother',1024);
    figure
    plot(wave)
    grid on
    xlabel('n')
    title(Family)
    fprintf('n=%i5 %s Wavelet \n\n',2^d,Type)
end

qmf = MakeONFilter(Type,par);
wc = FWT_PO(yp,1,qmf);

if DEBUG==1
    PlotWaveCoeff(wc,1,.5);
end

[n,J] = dyadlength(wc);
%for j=(J-1):-1:1  wc(dyad(j))=> wavelet coefficients at j-th resolution level
DWT_Spec=zeros(J-1,n);
for j=(J-1):-1:1
    V=abs(wc(dyad(j))).^2;
    nd=n/length(V);
    for nn=1:nd
      DWT_Spec(j,nn:nd:n)=V;
    end
end
return