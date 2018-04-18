function  [X0]=freqExtrapolate(U,p1,mu0,p2,rhoOveff,f1,f2)
%Extrapolate parameters from f1 to f2  
if mu0>=1
    Cpp=U;
else
    Cpp=U/mu0^(p2-p1);
end
Cpp_f2=Cpp*(f1/f2)^(2+(p1-1)/2);
mu0_f2=mu0*sqrt(f1/f2);
rhoOveff_f2=rhoOveff*sqrt(f1/f2);
if mu0_f2>=1
    U_f2=Cpp_f2;
else
    U_f2=Cpp_f2*(mu0_f2)^(p2-p1);
end
X0=[U_f2, p1, mu0_f2, p2, rhoOveff_f2];
return