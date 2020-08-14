function SDF=TwoComponentSDF(q,cs,p1,p2,q0)
%USAGE   SDF=TwoComponentSDF(q,cs,p1,p2,q0)
%
%PUURPOSE evaluate SDF(q)=cs*q^p1 for q<q0
%                         cs*q0^(p1/p2)*q^p2
%INPUTS:
%        q=spatial wavenumber         (vector)
%       cs=turbulent strength         (scalar)
%       p1=spectral index for q<=q0   (scalar)
%       p2=spectral index for q>q0    (scalar)
n=length(q); n1=find(q<=q0); 
n2=max(n1)+1:n;
if n1==1
    SDF=cs*abs(q).^(-p2);
else
    SDF=cs*abs(q(n1)).^(-p1);
    SDF=[SDF,cs*q0^(p2-p1)*abs(q(n2)).^(-p2)];
end
SDF(find(isnan(SDF)))=0;
return

    