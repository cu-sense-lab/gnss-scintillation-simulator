function [nfig1]=DisplayScaleSpectrum(SPEC_Summary,nskip,scani,PLOTID,SCALEID,scale,varargin)
%
nsegs=length(SPEC_Summary);
nseg_locs=(1:nsegs)*nskip-nskip+1;
ss=zeros(1,nsegs);
for nseg=1:nsegs
    if ~isempty(SPEC_Summary{nseg})
        nseg_scan=nseg;
        ss(nseg)=(scani(nseg_locs(nseg))+scani(nskip/2));
        xxPSD  =SPEC_Summary{nseg}.xxPSD;
        yyPSD  =SPEC_Summary{nseg}.yyPSD;
        xxScale=SPEC_Summary{nseg}.xxScale;
        yyScale=SPEC_Summary{nseg}.yyScale;
        offset =SPEC_Summary{nseg}.offset;
        if nseg==1
            SCALOGRAM=NaN(length(yyScale),nsegs);
            SPECTOGRAM=NaN(length(yyPSD),nsegs);
        end
        SPECTOGRAM(:,nseg)=yyPSD;
        SCALOGRAM(:,nseg)=yyScale;
    end
end

nfig1=figure;
imagesc(ss,scale,SCALOGRAM(:,1:nseg_scan))
colorbar
xlabel(SCALEID)
ylabel('log10(1/\lambda)')
caxis([SPECmax-50,SPECmax])
title([PLOTID,' SCALE SPEC'])
bold_fig
return