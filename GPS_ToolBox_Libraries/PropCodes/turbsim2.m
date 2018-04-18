 function simturb=turbsim2(rootSDF2)
%        simturb=turbsim2(rootSDF2)
%   
% Generate realization of field with SDF rootSDF.^2
%
[n1,n2]=size(rootSDF2);
xi=randn(n1,n2)+1i*randn(n1,n2);                       
simturb=2*real(fftshift(fft2(fftshift(rootSDF2.*xi)))); 
return