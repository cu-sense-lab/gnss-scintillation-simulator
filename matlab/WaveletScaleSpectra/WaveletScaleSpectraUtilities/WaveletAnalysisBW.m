function [SPEC_Summary,xscale,nskip,wc,qmf,J]=WaveletAnalysisBW(data,dx,jmin,Family,varargin)
%USAGE:  [SPEC_Summary,xscale,nskip,wc,qmf,J]=WaveletAnalysis(data,dx,jmin,Family,PLOTID,SCALEID,scale)
%USAGE:                                          ...,PLOTID,SCALEID,scale,DynRange,MaxOffset)
%PURPOSE: Segmented Wavelet Analysis
%
%INPUT:
%           data     = data array                      (NX1)
%           dx       = spatial resolution
%           jmin     = nextpow2(N)-jmin=largest DWT scale used in scale spectrum computation
%                      Number of segments: nsegs=2^(jmin-1);
%           
%           DynRange = Dynamic range for display       (dB)
%           MaxOffset= Offset from maximum for display (dB)
%           PLOTID   = Plot label                      (char)
%OUTPUT:
%           SPEC_Summary <= PSDScaleSpec                 (structure)
%                      
%
%SUBROUTINES:  ComputeDWT, ComputeScaleSpectrum, PSDScaleSpec
%
%Written by Charles Rino
%Rino Consulting
%January 2012 
%Revised August 2012
%crino@comcast.com
%http://chuckrino.com/wordpress
%
%'Haar', 'Daubechies', 'Coiflet','Symmlet'
if isempty(Family)
   Family='Daubechies';
end
if isempty(varargin)
    DISP=0;
elseif length(varargin)==3
    DISP=1;
    DynRange=100;
    MaxOffset=10;
    PLOTID=varargin{1};
    SCALEID=varargin{2};
    scale=varargin{3};
elseif length(varargin)==5
    DISP=1;
    PLOTID=varargin{1};
    SCALEID=varargin{2};
    scale=varargin{3};
    DynRange=varargin{4}; 
    MaxOffset=varargin{5};
else
    error('varargin error')
end
data=data(:)';
nsamp=length(data);
%Compute discrete wavelet transform
[DWT_Spec,wc,J,~,qmf]=ComputeDWT(data,Family);

%Compute spatial wavenumber scales & extension to nearest power of 2
sj=(2*dx)*2.^(J-2:-1:0);
qScale=1./sj;

%Extended scale for nearest power of 2
[~,nsamp_ext]=size(DWT_Spec);
x_data_ext=dx*(0:nsamp_ext-1);
xscale    =x_data_ext(1:nsamp);
if DISP~=0
    DWT_max=10*ceil(max(dB10(DWT_Spec(:)))/10)-MaxOffset;
    nfig=figure;   %Display DWT segmentation
    colormap('gray')
    imagesc(xscale*scale,log10(qScale),dB10(DWT_Spec(:,1:nsamp)))
    grid on
    xlabel(SCALEID)
    ylabel('log10(1/\lambda)')
    caxis([DWT_max-DynRange, DWT_max])
    colorbar
    title([PLOTID, ' DWT'])
end

%Compute scale spectrum
if isempty(jmin)
    fprintf('J=%5i \n',J)
    jmin=input('Input jmin (1<=jmin <J) ');
end
if jmin<1 || jmin>=J
    error('jmin out of range ')
end
[S_PLP,nsegs,~]=ComputeScaleSpectrum(J,wc,jmin);
nskip=nsamp_ext/nsegs;       %Number of samples per segment

if DISP~=0
    figure(nfig)
    hold on
    nn=0;
    for nseg=nskip/2:nskip:nsamp
        nn=nn+1;
        if nn==1
            ntics_max=length(nskip/2:nskip:nsamp);
            nmod=floor(ntics_max/8);
        end
        if mod(nn,nmod)==0,
            hold on,
            text(x_data_ext(nseg)*scale,log10(qScale(1)),num2str(nn),'color','k')
        end
    end
    bold_fig
end

nfft=nskip;
SPEC_Summary=PSDScaleSpec(data,nsamp,nfft,nsegs,dx,S_PLP,J,jmin);
return