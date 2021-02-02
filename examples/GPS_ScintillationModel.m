function     [S4,p1,mu0,p2,Imu,mu]=GPS_ScintillationModel(U,varargin)
%Phase-screen theory predicts the intensity SDF in normalized units as a
%function of U, p1, mu0=q0*rhoF, and p2
%The scintillation model sets representative values of p1, mu0, and p2 for the specified value of U
%User input values apply to entire U range
%For single value of U, model outputs Imu,mu
Imu=[]; mu=[];
if isempty(varargin)
    n1=find(U<0.45);
    n2=setdiff(1:length(U),n1);
    S4=zeros(size(U));
    p1=zeros(size(U));
    mu0=zeros(size(U));
    p2=zeros(size(U));
    for n=n1
        p1(n)=2.5;
        mu0(n)=.5;
        p2(n)=2.5;
        [~,~,S4(n),~,~,~]=Ispectrum(U(n),p1(n),p2(n),mu0(n));
    end
    for n=n2
        p1(n)=2.5;
        mu0(n)=.5;
        p2(n)=3.5;
        [Imu,mu,S4(n),~,~,~]=Ispectrum(U(n),p1(n),p2(n),mu0(n));
    end
else
    p1=varargin{1};
    mu0=varargin{2};
    p2=varargin{3};
end
return

