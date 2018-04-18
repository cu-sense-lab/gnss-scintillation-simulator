function IDstring=generateIDstring(U0,p1,p2,mu0,rhoF,S4)
%Generate ID character string
IDU=num2str(fix(U0*100)/100);
IDp1=num2str(fix(p1*100)/100);
IDmu0=num2str(fix(mu0*100)/100);
IDp2=num2str(fix(p2*100)/100);
IDrho=num2str(fix(rhoF*100)/100);
IDS4=num2str(fix(S4*100)/100);
ID=['U=',IDU,' \eta_1=',IDp1,' \mu_0=',IDmu0,' \eta_2=',IDp2];
if ~isempty(rhoF)
    IDstring=[ID,'\rho_F=',IDrho];
else
    IDstring=ID;
    IDstring=strrep(IDstring,'mu0','qo');
end
if ~isempty(S4)
    IDstring=[IDstring,' S4=',IDS4];
end
return
