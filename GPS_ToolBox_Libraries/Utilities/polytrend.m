function [y]=polytrend(data,t,norder)
%
[P,~,MU] = polyfit(t,data,norder);
y = polyval(P,t,[],MU);
return