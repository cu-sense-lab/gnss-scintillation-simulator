function   [Pmu,nb]=generatePmu(Cpp,p1,p2,mu0,mu)
%USAGE:  [Pmu,nb]=generatePmu(Cpp,p1,p2,mu0,mu)
%
%Generate P(mu) for parameters Cpp, p1, p2, mu0

nmu=length(mu);
Pmu=zeros(1,nmu);
nb=[];
for n=1:nmu
     Pmu1=abs(mu(n))^(-p1);
     Pmu2=mu0^(p2-p1)*abs(mu(n))^(-p2);
     if Pmu1<=Pmu2
         Pmu(n)=Cpp*Pmu1;
         nb=n;
     else
         Pmu(n)=Cpp*Pmu2;
     end
end
return