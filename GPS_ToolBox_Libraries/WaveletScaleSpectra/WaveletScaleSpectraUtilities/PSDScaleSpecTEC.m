function     [SPEC_Summary,J,nSmax]=PSDScaleSpecTEC(data,t_sec,t_seg,fileID,nseg,varargin)
%USAGE: SPEC_Summary=PSDScaleSpec
%
%OUTPUT:
%          SPEC_Summary{}.SPEC= nsegs cell array of summary data
%          SPEC fields
%               data_seg  =data samples
%               xxPSD     = log10(frequency scale)
%               yyPSD     = 10log10(PSD)
%               yy0       = least squares polynomial fit or PSD at scale spec samples
%               xxScale   = log10(DWT scales)
%               yyScale   = 10log10(ScaleSpec)
%               offset    = dB offset to align Scale Spec with PSD
%NOCLASSIFICATION
%
%SUBROUTINES: pdgm, dB10
%
if ~isempty(varargin)
    DISPLAY=varargin{1};
else
    DISPLAY=0;
    nfig1=[];
    nfig2=[];
end
[DWT,DWT_bar,fScale,J]=TECWaveletAnalysis(data,t_sec,nseg);
if DISPLAY==1 || DISPLAY==2
    [~,DWT_max,DynRange]=DisplayDWT(DWT,t_sec,fScale,nseg,t_seg,fileID);
end
if DISPLAY==2
    figure
    imagesc(t_seg/60,log10(fScale),dB10(DWT_bar))
    xlabel('t-min')
    ylabel('log10(1/s_j)')
    caxis([DWT_max-DynRange, DWT_max])
    colorbar
    title([fileID,'  ScaleSpec'])
end

[~,nsegs]=size(DWT_bar);
SPEC_Summary=cell(1,nsegs);
SPEC=struct('data_seg',[],'tspan',[],'xxPSD',[],'yyPSD',[],'xxScale',[],'yyScale',[]);
nsamp=length(data);
nsamp_seg=nseg(2)-1;
nfft=nicefftnum(nsamp_seg);
dt=diff(t_sec(1:2));
fspan=(-nfft/2:nfft/2-1)/(nfft*dt);
noff=0;
nSmax=nextpow2(length(data))-nextpow2(nseg(2));
for ns=1:nsegs
    n1=(ns-1)*nsamp_seg+1; n2=min(n1+nsamp_seg-1,nsamp);
    if n2>nsamp
        break
    end
    data_seg=data(n1:n2);
    y=(0:length(data_seg)-1)*dt;
    %Periodogram
    d=(data_seg(end)-data_seg(1))*(y-y(1))/(y(end)-y(1));
    data_seg=data_seg-d;
    psd=pdgm(data_seg,nfft,noff);
    xxPSD=log10(fspan(nfft/2+2:nfft));
    yyPSD=dB10(psd(2:nfft/2));
    %Scale spectrum interpolated to xxPSD
    xxScale=log10(fScale);               xxScale=xxScale(nSmax+1:end);
    yyScale=dB10(DWT_bar(:,ns))';  yyScale=yyScale(nSmax+1:end);
    SPEC.data_seg=data_seg;
    SPEC.tspan=y;
    SPEC.xxPSD   =xxPSD;
    SPEC.yyPSD   =yyPSD;
    SPEC.xxScale =xxScale;
    SPEC.yyScale =yyScale-dB10(nfft*dt);
    SPEC_Summary{ns}=SPEC;  
end %End of segment
return