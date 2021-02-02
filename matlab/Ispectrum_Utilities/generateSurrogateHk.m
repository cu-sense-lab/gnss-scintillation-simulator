function       [psi,phase0,fracMom,fDop,mu,SDF0]=...
                      generateSurrogateHk(U,p1,p2,mu0,rhoOveff,Dt,nfft,SEED,fscale,varargin)
%USAGE: [psi,phase0,fracMom,fDop,mu,SDF0]=generateSurrogateHk(U,p1,p2,mu0,rhoOveff,Dt,nfft,SEED,nint)
%
%Generate time-domain surrogate realization of complex field
%INPUT:
%Phase screen parameters=U, p1, p2, mu0, rhoOveff
%Intensity time sampling      =Dt sec.
%DFT sampling                    =nfft
%Random number seed      =SEED
%
%OUTPUT
%psi          =complex field (hk)
%phase0  =equivalent phase screen realization
%fracMom=(1X5) fractional moments 1 through 5  
%Note:S4=sqrt(fracMom(2)-1)
%fDop      =Doppler frequency range
%mu        =normalizae frequency
%

if isempty(varargin)
    nint=1;
else
    nint=varargin{1};
end
%Doppler sampling
fDop=(-nfft/2:nfft/2-1)/(nfft*Dt); 

%Convert Doppler to mu
mu=2*pi*fDop*rhoOveff;
mu_p=mu(nfft/2+2:nfft);

if mu0>=1
    Cpp=U;
else
    Cpp=U/mu0^(p2-p1);
end
%Generate SDF(mu)
[SDF0,~]=generatePmu(Cpp,p1,p2,mu0,mu);
SDF0(nfft/2+1)=0;
dmu=mu_p(1);

%Generate phase realization
rng(SEED);  
phase0=turbsim1(sqrt(SDF0*dmu/2/pi));
%Remove linear trend
y=linex(1:nfft,1,nfft,phase0(1),phase0(nfft));
phase0=phase0-y;
phase0=phase0*fscale;
%Phase screen simulation
psi=exp(1i*phase0);
%Propagation factor
mu=fftshift(mu);
pfac=exp(-1i*(mu).^2/2);   %Parabolic approximation 
psi_hat=fft(psi);
psi_hat=psi_hat.*pfac;
psi=ifft(psi_hat);
if nint>1
    psi=fftInterp(psi,nint);
    phase0=real(fftInterp(phase0,nint));
end

fracMom=generate_fracMom(abs(psi).^2);
return