function [ID,Cs1,ps1,Cs2,ps2,qB] = SpecType(P)
%SpecType
Cs1=P(2,1);  Cs2=P(2,2);
ps1=P(1,1)/10;  ps2=P(1,2)/10;
qB=log(Cs1/Cs2)/(ps2-ps1);
if ps1-ps2<=-0.2
    ID='Noise';
else
    if abs(Cs1-Cs2)<4 && abs(ps1-ps2)<0.2
        ID='One Component';
        Cs1=(Cs1+Cs2)/2; ps1=(ps1+ps2)/2;
        qB=0;
    else
        ID='Two Component';
    end
end
end

