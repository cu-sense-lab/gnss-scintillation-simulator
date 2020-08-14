function [nfig,DWT_max,DynRange]=DisplayDWT(DWT,t_sec,fScale,nseg,t_seg,fileID)
%
DynRange=80;
DWT_max=dB10(max(DWT(:)));
nfig=figure;   %Display DWT segmentation
imagesc(t_sec/60,log10(fScale),dB10(DWT))
grid on
xlabel('t-min')
ylabel('log10(1/s_j)')
caxis([DWT_max-DynRange, DWT_max])
colorbar
axisARG=[floor(min(t_sec/60)) ceil(max(t_sec/60)) min(log10(fScale)) max(log10(fScale))];
axis(axisARG)
title([fileID, ' DWT'])
hold on
for nn=1:4:length(nseg)
    hold on
    text(t_seg(nn)/60,log10(fScale(1))+0.25,num2str(nn))
end
bold_fig
return