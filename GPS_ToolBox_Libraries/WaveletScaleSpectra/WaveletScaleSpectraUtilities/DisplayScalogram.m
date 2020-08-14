function [nfig,ss,xxScale]=DisplayScalogram(SPEC_Summary,nskip,scani,PLOTID,SCALEID,scale)
%
nsegs=length(SPEC_Summary);
nseg_locs=(1:nsegs)*nskip-nskip+1;
for nseg=1:nsegs
    if ~isempty(SPEC_Summary{nseg})
        nseg_scan=nseg;
        xxScale=SPEC_Summary{nseg}.xxScale;
        yyScale=SPEC_Summary{nseg}.yyScale;
        if nseg==1
            SCALOGRAM=NaN(length(yyScale),nsegs);
        end
        SCALOGRAM(:,nseg)=yyScale;
    end
end

ss=(scani(nseg_locs(1:nseg_scan))+scani(nskip/2));
SPECmax=max(max(SCALOGRAM(:,1:nseg_scan)));
nfig=figure;
imagesc(ss*abs(scale),xxScale,SCALOGRAM(:,1:nseg_scan))
if scale>0
    colorbar
else
    colorbar('SouthOutside')
end
xlabel(SCALEID)
ylabel('log10(1/\lambda)')
caxis([SPECmax-100,SPECmax])
title([PLOTID,' SCALE SPEC'])
bold_fig