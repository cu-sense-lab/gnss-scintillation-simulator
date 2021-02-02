function [phase_d,Doppler,DIFF2]=Phase_Detrend(t_sec,phase,varargin)
%Outputs smoothed Doppler & Detrended phase
%
if isempty(varargin)
    norder=10;
else
    norder=varargin{1};
end
[nfreqs,nsamps]=size(phase);
dt=diff(t_sec(1:2));
Doppler=zeros(size(phase));
phase_d=zeros(size(phase));
DIFF2=zeros(size(phase));
for nfreq=1:nfreqs
    Doppler(nfreq,2:nsamps)=diff(phase(nfreq,:));
    Doppler(nfreq,1)=Doppler(nfreq,2);
    DIFF2(nfreq,2:nsamps)=diff(Doppler(nfreq,:)/dt)/(2*pi)^2;
    DIFF2(nfreq,1)=DIFF2(nfreq,2);
    [P,~,MU]=polyfit(t_sec,Doppler(nfreq,:),norder);
    y=polyval(P,t_sec,[],MU);
    Doppler(nfreq,:)=Doppler(nfreq,:)-y;
    for n=2:nsamps
        phase_d(nfreq ,n)=phase_d(nfreq,n-1)+Doppler(nfreq,n);
    end
    Doppler(nfreq,:)=y/(2*pi*dt);
end