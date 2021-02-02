function [nfig,dBmin]=DisplayGPSIntensity(tmin,IntL1d,IntL2d,IntL5d,PlotID,Jmax,varargin)
%
if ~isempty(varargin)
    dBspan=varargin{1};
else
    dBspan=40;
end
if ~isempty(Jmax)
    PlotIDD=[PlotID,' Detrend Jmax=',num2str(Jmax)];
else
    PlotIDD=[PlotID,' Detrend'];
end
%Display Detrended Intensity
if ~isempty(IntL5d);
    dBmax=dB10(ceil(max([max(IntL1d),max(IntL2d),max(IntL5d)]/10))*10);
    dBmin=dBmax-dBspan;
    nfig=figure;
    subplot(3,1,1)
    plot(tmin,dB10(IntL1d),'r')
    grid on
    ylabel('L1 dB')
    title(PlotIDD)
    axis([min(tmin) max(tmin) dBmin dBmax])
    
    subplot(3,1,2)
    plot(tmin,dB10(IntL2d),'r')
    grid on
    ylabel('L2 dB')
    axis([min(tmin) max(tmin) dBmin dBmax])
    
    
    subplot(3,1,3)
    plot(tmin,dB10(IntL5d),'r')
    grid on
    ylabel('L5 dB')
    axis([min(tmin) max(tmin) dBmin dBmax])
 
    xlabel('t-min')
    bold_fig
elseif ~isempty(IntL2d)
    dBmax=dB10(ceil(max([max(IntL1d),max(IntL2d)]/10))*10);
    dBmin=dBmax-dBspan; 
    nfig=figure;
    subplot(2,1,1)
    plot(tmin,dB10(IntL1d),'r')
    grid on
    ylabel('L1 dB')
    title(PlotIDD)
    axis([min(tmin) max(tmin) -50 20 ])
    
    subplot(2,1,2)
    plot(tmin,dB10(IntL2d),'r')
    grid on
    ylabel('L2 dB')
    axis([min(tmin) max(tmin) dBmin  dBmax])
    xlabel('t-min')
    bold_fig
else
    dBmax=dB10(ceil(max([max(IntL1d)]/10))*10);
    dBmin=dBmax-dBspan;
    nfig=figure;
    plot(tmin,dB10(IntL1d),'r')
    grid on
    ylabel('L1 dB')
    title(PlotIDD)
    xlabel('t-min')
    bold_fig
end
return