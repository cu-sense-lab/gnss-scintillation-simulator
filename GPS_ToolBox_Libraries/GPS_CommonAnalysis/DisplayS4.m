function [nfig]=DisplayS4(S4,t_sec,nfreqs,nsecs,nseg,nsegs,fileID)
%Display S4

nfig=figure;
plot(t_sec(nseg)/60+nsecs/60/2,S4(1,1:nsegs),'go-')
hold on
plot(t_sec(nseg)/60+nsecs/60/2,S4(2,1:nsegs),'bo-')
if nfreqs==3
    hold on
    plot(t_sec(nseg)/60+nsecs/60/2,S4(3,1:nsegs),'ro-')
    legend('L1','L2','L5')
else
    legend('L1','L2')
end
grid on
title(fileID)
ylabel('S4')
xlabel('t-min')
axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 0,1])
bold_fig

return