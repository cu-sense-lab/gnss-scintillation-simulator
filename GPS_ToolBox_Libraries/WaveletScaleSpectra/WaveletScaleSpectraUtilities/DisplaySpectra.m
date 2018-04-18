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
        yy0    =SPEC_Summary{nseg}.yy0;
        yyScale=SPEC_Summary{nseg}.yyScale;
        if isfield(SPEC_Summary{nseg},'Type')
            %ID=SPEC_Summary{nseg}.ID;
            Type   =SPEC_Summary{nseg}.Type;
            if ~isempty(Type)
                Cs1=Type.Cs1; ps1=Type.ps1;
                Cs2=Type.Cs2; ps2=Type.ps2;
                ns_max=length(xxScale);
                ns=Type.ns;
                sig_Spec=min(Type.sig);
                yy1 = Cs1-10*ps1*xxScale(1:ns);
                yy2 = Cs2-10*ps2*xxScale(ns+1:ns_max);
                fprintf('\n #1 Cs=%8.2f p=%5.2f ns=%4i #2 Cs=%8.2f p=%5.2f',...
                    Cs1,ps1,ns,Cs2,ps2);
                yySDF_model=[yy1,yy2];
            else
                sig_Spec=[];
            end
        end
        figure(nfig);
        hold off
        plot(xxPSD,yyPSD,'c')
        %hold on
        %plot(xxScale,yy0,'g')
        hold on
        plot(xxScale,yyScale,'b')
        if ~isempty(Type)
            hold on
            plot(xxScale(1:ns),yySDF_model(1:ns),'r')
            hold on
            plot(xxScale(ns+1:ns_max),yySDF_model(ns+1:ns_max),'r')
        end
        hold on
        text(xxPSD(end)+0.2,yyPSD(end),num2str(nseg));
        xlabel('log10(1/\lambda)')
        ylabel('SDF-dB')
        %legend('PSD','PSD fit','Scale','Model')
        legend('PSD','Scale','log-linear')
        grid on
        title([PLOTID,'  sig-Spec=',num2str(sig_Spec)])
        bold_fig
        %drawnow
        %if ~isempty(PLOTID)
        %    ID=[PLOTID,num2str(nseg)];
        %  saveas(nfig,ID,'jpg')
        %end
    end
    fprintf('\n')
end
return