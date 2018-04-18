function [nfig]=DisplayDopplerDIFF2(Doppler,DIFF2,t_sec,nfreqs,fileID)
%Display Doppler &Diff2

Col{1}='r'; Col{2}='b'; Col{3}='g';
nfig=figure;
for nfreq=nfreqs:-1:1
    subplot(2,1,1)
    hold on
    plot(t_sec/60,Doppler(nfreq,:)/1000,Col{nfreq})
    subplot(2,1,2)
    hold on
    plot(t_sec/60,DIFF2(nfreq,:),Col{nfreq})
end
subplot(2,1,1)
grid on
title(fileID)
ylabel('\Delta\phi-kHz')
axis([floor(min(t_sec/600))*10  ceil(max(t_sec/600))*10 min(Doppler(:))/1000  max(Doppler(:))/1000])
subplot(2,1,2)
grid on
ylabel('Diff2-cycles')
axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 -0.5,0.5])
if nfreqs==3
    legend('L5','L2','L1')
else
    legend('L2','L1')
end
bold_fig