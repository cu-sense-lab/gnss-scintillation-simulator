function pdf=pdfLogN(p,xbar,sigx)
%
%   Log Normal Distribution 
%
%INPUTS:
%   p    = power (vector)
%   pbar = <p> 
%   sigx = std(log(p))
%
iOK=find(p>0);
pdf=zeros(size(p));
c=1/sqrt(2*pi*sigx^2);
pdf(iOK)=c./p.*exp(-(log(p(iOK))/2-xbar).^2/(2*sigx^2));
return