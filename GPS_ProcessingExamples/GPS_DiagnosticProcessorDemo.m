%DiagnosticProcessorDemo
%Demonstrate basic preprocessing operations

%Select active data_matfile
if exist('data_matfile','var')~=1 || ~exist('OK','var')
    clear all
    close all
    OK=1;
else
    fprintf('Using current active file %s \n',data_matfile)
    OK=input('Input 1 for new file CR to continue ');
end

if ~isempty(OK)
    clear all
    close all
    %DEFINE PARAMETERS AND VARIABLES
    dtr=pi/180;
    %frequencies
    freqL1=154*10.23e6;  freqID{1}='L1';
    freqL2=120*10.23e6;  freqID{2}='L2';
    freqL5=115*10.23e6;  freqID{3}='L5';
    freqs=[freqL1,freqL2,freqL5];
    
    c = 299792458;            %Speed of light (vacuum)
    re=2.819740289e-15;   %Classical electron radius m
    K=re*c/(2*pi)*1.e16;     %TEC conversion factor
    
    %USER SET PARAMETERS
    nsecs=600;   %number of seconds per segment  (600 =>10 min)
    dt=1/100;      %Phase sample interval
    Jmin=6;         %Minimum wavelet scale for detrend (smaller is smoother detrend Jmin=d => no detrend)
    %                                                          d is maximum power of 2 larger than total number of samples
    %SELECT PRECOMPUTED DATA FILE
    data_dir=pwd;
    data_dir=strrep(data_dir,'GPS_ProcessingExamples','GPS_Data_Examples');
    data_matfile=selectmfile([data_dir,'\*.mat']);
    load(fullfile(data_dir,data_matfile))
    
    fileID=strrep(data_matfile,'.mat',' ');
    fileID=strrep(fileID,'_','-');
    if contains(data_matfile,'Peru')
            StationID='Peru';
    elseif contains(data_matfile,'HongKong')
            StationID='HongKong';
    else
            error('1 or 2 only for Demo')
    end
    OK=1;
    %PERFORM PREPROCESSING OPERATIONS
    GPS_PreprocessorDemo
end
%OUTPUTS:
%intensity, phase, prnge, nfreqs, t_sec, t_secG, fileID,  sat_rnge, sat_rdot, sat_elev, sat_phi, sat_llh, sat_tcs, GPSTime, eph, origin_llh
%nsamp_seg, nseg, nsegs  <=nsecs (# sec in segment
%intensity_d, intensity_bar, phase_d, Doppler,DIFF2
%S4
DISPLAY=input('Input 1 to display raw intensity and mean ');
if ~isempty(DISPLAY)
    DisplayIntensity(intensity,intensity_bar,t_sec,nfreqs,freqID,fileID);
end

DISPLAY=input('Input 1 to display detrended intensity & phase ');
if ~isempty(DISPLAY)
    DisplayIntensityPhase_d(intensity_d,phase_d,t_sec,nfreqs,freqID,fileID,nseg,nsecs);
    %DisplayIntensityPhase_d(intensity_d,[],t_sec,nfreqs,freqID,fileID,[],nsecs);
end

DISPLAY=input('Input 1 to S4 ');
if ~isempty(DISPLAY)
    DisplayS4(S4,t_sec,nfreqs,nsecs,nseg,nsegs,fileID);
end

DISPLAY=input('Input 1 to diplay Doppler, DIFF2 ');
if ~isempty(DISPLAY)
    DisplayDopplerDIFF2(Doppler,DIFF2,t_sec,nfreqs,fileID);
end
        
%DISPLAY=input('Input 1 for TEC analysis ');
if ~isempty(DISPLAY)
    DisplaySTEC(phase,t_sec,nfreqs,fileID);
end

OK=input('Input 1 to generate Geometry ');
if ~isempty(OK)
    h_intercept=350000;
    Drift=[0,0,0];
    [satGEOM]=ExtractGeometry(fullfile(data_dir,data_matfile),StationID,h_intercept,Drift);
    DISPLAY=input('Input 1 to display geometry ');
    if ~isempty(DISPLAY)
        DisplayGeometry(satGEOM,fileID)
    end
end
