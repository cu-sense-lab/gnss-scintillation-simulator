function      [SDF_P]=generateSDF_P_TEC(Cp,p1,p2,f0,fDop_P)
%USAGE:   [SDF_P,nb]=generateSDF_P_TEC(Cp,p1,p2,f0,fDop_P)
%Generate SDF(fDop) for parameters Cp, p1, p2, f0
%
 if isempty(p2)
     p2=p1;
 end
npts=length(fDop_P);
SDF_P=zeros(1,npts);

n1=fDop_P <=f0 ;
SDF_P(n1)=fDop_P(n1).^(-p1);
n2= fDop_P >f0 ;
SDF_P(n2)=f0^(p2-p1)*abs(fDop_P(n2)).^(-p2);
SDF_P=Cp*SDF_P;
return