function [phase,diffW,Wdiff]=unwrapPhase(Wphase)
%Same a MatLab unwrap  Itoh's algorithm
diffW=[0,diff(Wphase)];
Wdiff=atan2(sin(diffW),cos(diffW));
npts=length(Wphase);
phase=zeros(size(Wphase));
for n=2:npts
    phase(n)=phase(n-1)+Wdiff(n);
end
return