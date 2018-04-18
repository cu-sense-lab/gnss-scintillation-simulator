function  [nfig]=DisplayIntensity(intensity,intensity_bar,t_sec,nfreqs,freqID,fileID)
%Display Intensity

nfig=figure;
dBmax=max(dB10(intensity(:)));
for nfreq=1:nfreqs
    subplot(nfreqs,1,nfreq)
    plot(t_sec/60,dB10(intensity(nfreq,:)),'r')
    if ~isempty(intensity_bar)
        hold on
        plot(t_sec/60,dB10(intensity_bar(nfreq,:)),'b')
    end
    grid on
    ylabel([freqID{nfreq},' CN0-dB'])
    axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 floor((dBmax-20)/10)*10,ceil((dBmax+10)/10)*10])
    if nfreq==1
        title(fileID)
        if ~isempty(intensity_bar)
            legend('Intensity','Mean')
        end
    end
end
xlabel('UT-min')
bold_fig
return