 function simturb=turbsim1(rootSDF1)
%        simturb=turbsim1(rootSDF1)
%   
% Generate realization of field with SDF rootSDF.^2
%
n2=length(rootSDF1);
xi=(randn(1,n2)+1i*randn(1,n2));  
simturb=sqrt(2)*real(fftshift((fft(fftshift(rootSDF1.*xi)))));
return