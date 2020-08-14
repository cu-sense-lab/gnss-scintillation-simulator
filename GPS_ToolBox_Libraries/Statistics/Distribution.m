function [PDF,CDF,yU]=Distribution(Y,varargin)
%USAGE: [PDF,CDF,yU]=Distribution(Y,varargin)
%  Compute CDF & PDF of random sequence Y
%  Repeated values are replaced with midpoints
%  PDF is unscaled => sum(PDF)=1;
%  ymid values are uniformly spaced  
%  varargin=number of PDF estimates averaged


% Written by Chuck Rino
% Rino Consulting
% April 2, 2007

if ~isempty(varargin)
    ns=varargin{1};
else
    ns=1;
end
a=1;
b=ones([1,ns])/ns;

y=sort(Y(:));
CDF=[0:length(y)-1]/length(y); %Eliminate multiple values
iStep=find(diff([0;y])>0);
if iStep==1
    fprintf('Singular Distribution \n')
    CDF=[];
    PDF=[];
    yU=[];
    return
end
ymid=(y(iStep(2:end))+y(iStep(1:end-1)))/2;
CDFm=CDF(iStep(2:end));
%Interpolate to uniform sampling
yU=linspace(ymid(1),ymid(end),length(ymid));
%Test for discrete random variables
[ymid,nU]=unique(ymid);
CDFm=CDFm(nU);
CDF=interp1(ymid,CDFm,yU);
PDF=diff(filter(b,a,[0,CDF]));
return