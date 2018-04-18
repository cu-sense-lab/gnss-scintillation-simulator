%GPS_Preprocessor

%Select active data_matfile
if exist('data_matfile','var')~=1 || ~exist('OK','var')
    error('Run GPS_DiagnosticProcessorDemo first \n')
end

nsamps=length(t_sec);
fprintf('\n dt=0.01, %5.2f-min segment length\n log2(nsamps)=%5.2f Jmin=%3i \n',nsecs/60,log2(nsamps),Jmin)
nfreqs=3;
        
%Offset phase for alignment with pseudo range (Coarse resolution of phase M ambiguity)
for nfreq=1:nfreqs
    offset=mean(phase(nfreq,1:100:end)*c/(2*pi*freqs(nfreq))-prnge(nfreq,:));
    phase(nfreq,:)=phase(nfreq,:)-offset*2*pi*freqs(nfreq)/c;
end

DISPLAY=input('Input freq # for alignment check or CR ');
if ~isempty(DISPLAY)
    %Verify alignment & and compare to ephemeris range
    nfreq=DISPLAY;
    figure
    plot(t_sec(1:100:nsamps)/60,phase(nfreq,1:100:nsamps)*c/(2*pi*freqs(nfreq))/1000,'r')
    hold on
    plot(t_sec(1:100:nsamps)/60,prnge(nfreq,:)/1000,'b')
    grid on
    legend('phase','prnge')
    ylabel('Range-km')
    xlabel('t-min')
    title([fileID,' ',freqID{nfreq}])
    bold_fig
end

%%%%%%%%Segmentation%%%%%%%%%%%%%%%%%%%%%%%%
d=nextpow2(nsamps);

nsamp_seg=nicefftnum(nsecs/dt);       %Samples per segment  nsecs is user defined
nseg=(1:nsamp_seg:length(t_sec));    %Sample start at each segment
nsegs=length(nseg);
fprintf('Segment length=%8.2f min \n',nsamp_seg*dt/60)

%Check to see if last segment is fully populated
if  (nseg(end)+nsamp_seg)>length(t_sec)
    nover=nseg(end)+nsamp_seg-length(t_sec);
    nseg=1:nsamp_seg:nsamp_seg*(nsegs-1);
    nsegs=length(nseg);
end
if nsegs==2
    nseg=nseg+nsamp_seg-nover;
end

%Detrend Intensity & Phase
[intensity_d,intensity_bar]=Intensity_Detrend(intensity,Jmin);
[phase_d,Doppler,DIFF2]=Phase_Detrend(t_sec,phase);

%Calculate S4
for nfreq=1:nfreqs
    [S4temp,Imean,nstep,fracMom]=computeSI(intensity_d(nfreq,:),nsamp_seg,nsamp_seg);
    if nfreq==1
        S4=zeros(nfreqs,length(S4temp));
    end
    S4(nfreq,:)=S4temp;
end
