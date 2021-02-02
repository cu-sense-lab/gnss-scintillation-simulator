function [Exc,ylin,n]=exceedance(I,varargin)
%  [Exc,ylin,n]=exceedance(I,varargin)
%  Computes excedance versus I or 10log10(I) if second argument is nonzero
%  Duplicate Exceedances are removed with I or logI truncated accordingly
%

% Written by Chuck Rino
% Vista Research, Inc.
% October 2006
if isempty(varargin)
    dBout=0;
else
    dBout=varargin{1};
end
if dBout==1
    y=sort(dB10(I(:))); 
else
    y=sort(I(:));
end
%The resort eliminates multiple Exc values
%that cause subsequent interpolation to bomb
[yu,j]=sort(y);
n=length(yu);
ylin=yu;
Exc=1-j/n;
return