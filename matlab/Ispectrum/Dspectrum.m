function    [Dmu,mu,result]=Dspectrum(U,p1,p2,mu0,varargin)
%USAGE [Dmu,mu,result]=Dspectrum(U,p1,p2,mu0,varargin)
%Compute Dspectrum at Ispectrum mu values
if isempty(varargin)
    [~,mu]=Ispectrum(U,p1,p2,mu0);
else
    mu=varargin{1};
end
nsteps=length(mu);
DspecParams=generateIspecParams(U,p1,p2,mu0);
Dmu=zeros(1,nsteps);
muD=zeros(1,nsteps);
for nmu=1:nsteps
    DspecParamsM=[DspecParams,' ',num2str(mu(nmu))];
    [status,result]= system([[getGlobalpath2exe,'\Dspectrum.exe '],DspecParamsM]);
    if status~=0
        error('result')
    end
    data=importdata([pwd,'\dspectrum.dat']);
    [~,ndata]=size(data);
    if ndata~=3
        fclose('all');
        error('Dspectrum fault ')
    else
        muD(nmu)=data(:,1);  %Check should be mu
        Dmu(nmu)=data(:,2);
    end
    fclose('all');
end
Dmu=Dmu';
delete([pwd,'\dspectrum.dat'])
delete([pwd,'\dspectrum.log'])
return
