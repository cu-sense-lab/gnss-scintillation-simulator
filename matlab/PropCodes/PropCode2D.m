function  [psi,IStatSum,phase0]=PropCode2D(k,x,xProp,y,thetap,StructureParams,PhaseScreen,SF)
%USAGE:    [psi,IStatSum,phase0]=PropCode2D(k,x,xProp,y,thetap,StructureParams,PhaseScreen,SF)
%
%INPUT:        k         =wavenumber
%                    x         =propagation distance in layer
%                    xProp =Free space propagation steps 
%                    y         =distance array
%                    thetap=propagation angle (rad)
%                    StructureParams=Structure with fiellds Cp,p1,q0,p2
%                    PhaseScreen=0 for FPE integration 1 for Phase screen
%                    SF      =scale factor
%OUTPUT:
%                    psi      =complex field at max(xProp)
%                    IStatSum = cell array statistics summary versus [x,xProp]
%                    phase0    = integrated phase

%Written by Chuck Rino 
%                   Rino Consulting
%                   May 7, Version
%

%Phase screen propagation code
nPropSteps=length(xProp);
nx=length(x); dx=diff(x(1:2));
ny=length(y); dy=diff(y(1:2));
%Spectral domain sampling
dky=2*pi/ny/dy; ky=(-ny/2:ny/2-1)*dky; 

%Propagation direction 
c_thetap=cos(thetap); s_thetap=sin(thetap); t_thetap=tan(thetap);

%Propagation factor
gkx=sqrt(1-(ky/k+s_thetap).^2);
pfac = fftshift(exp(1i*(k*gkx-k*c_thetap+t_thetap*ky)));

%Split-step FPE Integration
fprintf('Split step FPE Integration %5i steps \n',nx)
if PhaseScreen==1
    fprintf('Phase Screen ON \n')
end

IStatSum=cell(1,nx+nPropSteps);
fracMom=ones(1,5);
IStatSum{1}=fracMom;
psi=ones(size(y));
fprintf('\nFPE Integration %5i Steps\n', nx)
for nstep=1:nx
   showprogress(nstep,1)
    if PhaseScreen~=1  
        if nstep>1
            %Propagation x(nstep-1) to x(nstep)
            psi_hat=fft(psi);
            psi_hat=psi_hat.*pfac.^dx;
            psi=ifft(psi_hat);
            fracMom=generate_fracMom(abs(psi).^2);
            IStatSum{nstep}=fracMom;
        end
    end
     phase=generatePhasePerturbation1D(nstep,ky,SF,thetap,StructureParams);
    if nstep==1
        phase0=phase;
    else
        phase0=phase0+phase;
    end
    psi=psi.*exp(1i*phase);
end
if PhaseScreen==1
    for nstep=1:nx
        IStatSum{nstep}=ones(1,5);
    end
end
fprintf('\n')

fprintf('Free-space propagation %5i Steps\n', nPropSteps-1)
for nProp=1:nPropSteps-1;
    showprogress(nProp,1)
    %Free space propagation
    psi_hat=fft(psi);
    psi_hat=psi_hat.*pfac.^(xProp(nProp+1)-xProp(nProp));
    psi=ifft(psi_hat);
    fracMom=generate_fracMom(abs(psi).^2);
    IStatSum{nx+nProp-1}=fracMom; 
end
fprintf('\n')
return