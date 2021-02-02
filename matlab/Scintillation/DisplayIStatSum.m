function [SImax,nSImax]=DisplayIStatSum(IStatSum,r,nx,PlotID,varargin)
%
nPropSteps=length(r);
ndec=floor(nPropSteps/100);
if isempty(varargin)
    xx=r;
else
    k=varargin{1};
    if k>0
        xx=sqrt((r-r(1))/varargin{1})*1000;
    else
         xx=sqrt((r(64:end)-r(64))/abs(varargin{1}))*1000;
         xx=[zeros(1,63),xx];
    end
end
figure
SI=zeros(1,min(nPropSteps,100)); nPlot=0;
for nPropStep=1:ndec:nPropSteps
    fracMom=IStatSum{nPropStep};
    nPlot=nPlot+1;
    SI(nPlot)=sqrt(max(0,fracMom(2)-fracMom(1)));
    subplot(2,1,1)
    hold on
    if nPropStep<=nx
         plot(xx(nPropStep)/1000,SI(nPlot),'ro')
    else
        plot(xx(nPropStep)/1000,SI(nPlot),'bo')
    end
    subplot(2,1,2)
    plot(fracMom(2),log10(fracMom(3)),'r.')
    hold on
    plot(fracMom(2),log10(fracMom(4)),'b.')
    hold on
    plot(fracMom(2),log10(fracMom(5)),'m.') 
    hold on
    plot(2,log10(factorial(3)),'rp')
    hold on
    plot(2,log10(factorial(4)),'bp')
    hold on
    plot(2,log10(factorial(5)),'mp')
end
[SImax,nSImax]=max(SI);
subplot(2,1,1)
grid on
ylabel('SI')
if isempty(varargin)
   xlabel('x-km')
   axis([min(xx)/1000 max(xx)/1000 0 1.5])
else
   xlabel('\rho_F')
   axis([0 max(xx)/1000 0 1.5])
end
title(PlotID)
subplot(2,1,2)
grid on
ylabel('log10(FM(m))')
xlabel('FM(2)')
axis([1 3 0 5])
legend('m=3','m=4','m=5')
bold_fig

if 0
    figure
    plot(xx/1000,SI,'r')
    grid on
    ylabel('SI')
    if isempty(varargin)
        xlabel('x-km')
        axis([min(xx)/1000 max(xx)/1000 0 1.5])
    else
        xlabel('\rho_F')
        axis([0 max(xx)/1000 0 1.5])
    end
    title(PlotID)
    bold_fig
end
return
    
    
  
    