function [nfig]=DisplayScaleSpec(ScaleSpec,y,nseg_locs,qScale,PlotID,varargin)
%USAGE:      DisplayScaleSpec(ScaleSpec,y,nseg_locs,qScale,PlotID,varargin)
%
if isempty(varargin)
    MaxOffset=10;
    DynRange=100;
else
    MaxOffset=varargin{1};
    DynRange=varargin{2};
end

nsegs=length(nseg_locs);
ss=y(nseg_locs);
SPECmax=dB10(max(max(ScaleSpec(:,1: nsegs))))-MaxOffset;
nfig=figure;
imagesc(ss/1000,log10(qScale),dB10(ScaleSpec))
colorbar
xlabel('y-km')
ylabel('log10(1/\lambda)')
caxis([SPECmax-DynRange,SPECmax])
title([PlotID,' SCALE SPEC'])
bold_fig
return