function nfig=DisplaySpectra(SPEC_Summary,varargin)
%USAGE:  nfig=DisplaySpectra(SPEC_Summary,PLOTID)
%
if isempty(varargin)
  nstart=1; 
  nstep=1;
  nend=length(SPEC_Summary);
  PLOTID=[];
elseif length(varargin)==1
  nstart=1; 
  nstep=1;
  nend=length(SPEC_Summary);
  PLOTID=varargin{1};
else
  %Input start, step end & Add to currently open figure
  nstart=varargin{1};
  nstep=varargin{2};
  nend =min(varargin{3},length(SPEC_Summary));
  PLOTID=varargin{4};
end
nfig=figure;
nsegs=length(nstart:nstep:nend);
for nseg=nstart:nstep:nend;
    if isempty(SPEC_Summary{nseg})
        %Empty cell => end of segmentation
        break
    end
    fprintf('\n Seg %3i of %5i \n',nseg,nsegs)
    OK=input('1 to plot CR to skip -1 to quit ');
    if ~isempty(OK)
        if OK<0
            break
        end
        xxPSD  =SPEC_Summary{nseg}.xxPSD;
        yyPSD  =SPEC_Summary{nseg}.yyPSD;
        xxScale=SPEC_Summary{nseg}.xxScale;
        yyScale=SPEC_Summary{nseg}.yyScale;
        
        figure(nfig);
        hold off
        plot(xxPSD,yyPSD,'c')
        hold on
        plot(xxScale,yyScale,'b')
        
        hold on
        text(xxPSD(end)+0.2,yyPSD(end),num2str(nseg));
        xlabel('log10(1/\lambda)')
        ylabel('SDF-dB')
        legend('PSD','Scale')
        grid on
        title([PLOTID,'  sig-Spec=',num2str(nseg)])
        bold_fig
    end
    fprintf('\n')
end
return