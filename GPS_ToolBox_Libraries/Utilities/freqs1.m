function f=freqs1(nsamp, fspan, iout)
%      Generate DFT frequency vector
%USAGE:  f=freqs1(nsamp, fspan, iout)
%INPUTS:
%   nsamp=  number of samples
%   fspan=  frequency span     (2*fNyq=1/dt)
%   iout=0 => natural order    (-fNyq:fNyq-df)
%   iout=1 => fft order        ([0:fNyq,-fNyq:-df])
%OUTPUTS:
%   f= frequency vector

%
%  Chuck Rino
%  Rino Consulting
%  July 2010

if 2*fix(nsamp/2)==nsamp
    %nsamp is even
    ff=([0:nsamp-1]/nsamp-1/2)*fspan;
    if iout==1
        f(1:nsamp/2)=ff(nsamp/2+1:nsamp);
        f(nsamp/2+1:nsamp)=ff(1:nsamp/2);
    else
        f=ff;
    end
else
    %nsamp is odd
    ff=([0:nsamp-1]-(nsamp-1)/2)*fspan/nsamp;
    if iout==1
        f(1:(nsamp-1)/2+1)=ff((nsamp-1)/2+1:nsamp);
        f((nsamp-1)/2+2:nsamp)=ff(1:(nsamp-1)/2);
    else
        f=ff;
    end
end
return
