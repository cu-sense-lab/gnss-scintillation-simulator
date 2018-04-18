function [nfig]=DisplaySTEC(phase,t_sec,nfreqs,fileID)

%frequencies
freqL1=154*10.23e6;  freqID{1}='L1';
freqL2=120*10.23e6;  freqID{2}='L2';
freqL5=115*10.23e6;  freqID{3}='L5';
freqs=[freqL1,freqL2,freqL5];

c = 299792458;            %Speed of light (vacuum)
re=2.819740289e-15;   %Classical electron radius m
K=re*c/(2*pi)*1.e16;     %TEC conversion factor

if nfreqs==3
    %Geometry Free Combinations
    STEC12=(phase(1,:)/freqL1-phase(2,:)/freqL2)/(1/freqL2^2-1/freqL1^2)/(2*pi*K);
    STEC15=(phase(1,:)/freqL1-phase(3,:)/freqL5)/(1/freqL5^2-1/freqL1^2)/(2*pi*K);
    STEC12=STEC12-min(STEC12);
    STEC15=STEC15-min(STEC15);
    
    nfig=figure;
    subplot(2,1,1)
    plot(t_sec/60,STEC12,'r')
    grid on
    hold on
    plot(t_sec/60,STEC15,'b')
    title(fileID)
    legend('L12','L15')
    ylabel('\Delta STEC')
    axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 0,ceil(max([STEC12,STEC15])/10)*10])
    
    subplot(2,1,2)
    STECerr=STEC12-STEC15;
    plot(t_sec/60,STECerr,'r')
    grid on
    ylabel('Error')
    xlabel('t-min')
    bold_fig
    axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 floor(min(STECerr)), ceil(max(STECerr))])
else
    STEC12=(phase(1,:)/freqL1-phase(2,:)/freqL2)/(1/freqL2^2-1/freqL1^2)/(2*pi*K);
    STEC12=STEC12-min(STEC12);
    figure
    plot(t_sec/60,STEC12,'r')
    grid on
    title(fileID)
    legend('L12')
    ylabel('\Delta STEC')
    axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 0,ceil(max([STEC12])/10)*10])
end
bold_fig

STEC=STEC12;
Jmin=8;
STEC_bar=WaveletFilter(STEC,Jmin);
deltaTEC=STEC-STEC_bar;
figure
subplot(2,1,1)
plot(t_sec/60,STEC ,'b')
hold on
plot(t_sec/60,STEC_bar,'r')
grid on
ylabel('TEC');
title(fileID)
axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 0,ceil(max([STEC12])/10)*10])

subplot(2,1,2)
plot(t_sec/60,deltaTEC,'r')
grid on
ylabel('\Delta TEC')
xlabel('t-min')
title(['TEC residual ',fileID])
axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 floor(min(deltaTEC)) ,ceil(max(deltaTEC))])
bold_fig
