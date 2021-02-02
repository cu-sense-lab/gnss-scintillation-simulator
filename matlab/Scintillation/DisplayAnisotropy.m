function  DisplayAnisotropy(A,B,C,Y,Z,scale)
gcf;
D=A*C-B^2/4;
if D<0
    error('D<0')
end
nsamp=500;
phi =linspace(0,2*pi,nsamp);
Cphi=cos(phi); Sphi=sin(phi);
DEN=C*Cphi.^2+B*Cphi.*Sphi+A*Sphi.^2;
rho=scale*sqrt(D./DEN);
hold on
plot(Y+rho.*Cphi,Z+rho.*Sphi,'r')
return


