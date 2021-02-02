function g=EulerC
%Euler's constant to as many decimal places as zeros
ndec=7;
n=10^ndec;
k=[1:n];
arg=1./k;
g=sum(arg)-log(n);
return