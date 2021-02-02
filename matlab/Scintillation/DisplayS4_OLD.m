function nfig=DisplayS4(IStatSum,x,PlotID,varargin)
%
nx=length(x);
SI=zeros(1,nx);
for nStep=1:nx
    fracMom=IStatSum{nStep};
    SI(nStep)=sqrt(max(0,fracMom(2)-fracMom(1)));
end
if isempty(varargin)
    nfig=figure;
else
    figure(varargin{1})
end
plot(x/1000,SI,'r')
grid on
xlabel('x-km')
ylabel('S4')
title(PlotID)
bold_fig

return


