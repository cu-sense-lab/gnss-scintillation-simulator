function [nfig]=DisplayDWT(DWT_Spec,y,y_ext,q,PlotID,nseg_locs,nskip,varargin)
%USAGE:      DisplayDWT(DWT_Spec,y,y_ext,q,PlotID,nseg_locs,nskip,varargin)

if isempty(varargin)
    MaxOffset=10;
    DynRange=100;
else
    MaxOffset=varargin{1};
    DynRange=varargin{2};
end

ny=length(y);
DWT_max=10*ceil(max(dB10(DWT_Spec(:)))/10)-MaxOffset;
nfig=figure;   %Display DWT segmentation
imagesc(y/1000,log10(q),dB10(DWT_Spec(:,1:ny)))
grid on
xlabel('y-km')
ylabel('log10(1/\lambda)')
caxis([DWT_max-DynRange, DWT_max])
colorbar
title([PlotID, ' DWT'])
hold on
nn=0;
for nseg=nseg_locs
    nn=nn+1;
    if nn==1
        ntics_max=length(nskip/2:nskip:ny);
        nmod=floor(ntics_max/8);
    end
    if mod(nn,nmod)==0
        hold on
        text(y_ext(nseg)/1000,log10(q(1)),num2str(nn))
    end
end
return