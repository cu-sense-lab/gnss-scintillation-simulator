function [f,malias]=alias(f,fNyquist)
%
malias=0;
outside=1;
while outside==1
    low=find(f<-fNyquist); L=0;
    if ~isempty(low)
        malias=malias+1;
        f(low)=f(low)+2*fNyquist;
        L=1;
    end
    high=find(f>fNyquist); H=0;
    if ~isempty(high)
        f(high)=f(high)-2*fNyquist;
        H=1;
        malias=malias-1;
    end
    outside=max(L,H); 
end
return