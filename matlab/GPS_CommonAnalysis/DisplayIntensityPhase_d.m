function  [nfig]=DisplayIntensityPhase_d(intensity_d,phase_d,t_sec,nfreqs,freqID,fileID,nseg,nsecs)
%Display detrended Intensity
nfig=figure;
if isempty(phase_d)
    ncol=1;
else
    ncol=2;
end
figure
dBmin=min(dB10(intensity_d(:)));
for nfreq=1:nfreqs
    if ncol==1
        nloc=nfreq;
    else
        nloc=ncol*nfreq-1;
    end
    subplot(nfreqs,ncol,nloc)
    plot(t_sec/60,dB10(intensity_d(nfreq,:)),'r')
    grid on
    ylabel([freqID{nfreq}])
    axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 floor((dBmin-5)/10)*10,10])
    if nfreq==1
        if ncol==2
            title(['                                           ',fileID])
        else
            title(fileID)
        end
    end
    if nfreq==nfreqs && ~isempty(nseg)
        hold on
        plot(t_sec(nseg)/60+nsecs/60/2, floor((dBmin-5)/10)*10,'bo')
    end
end
if ncol==2
   xlabel('                                                UT-min')
else
    xlabel('UT-min')
end

if ncol==2
    phase_max=max(abs(phase_d(:)));
    phase_max=ceil(phase_max/10)*10;
    
    for nfreq=1:nfreqs
        subplot(nfreqs,2,2*nfreq)
        plot(t_sec/60,phase_d(nfreq,:),'r')
        grid on
        axis([floor(min(t_sec/600))*10 ceil(max(t_sec/600))*10 -phase_max phase_max])
    end
end
bold_fig