function         [psi,IStatSum,phase0]=PropCode3D(k,x,xProp,y,z,thetap,phip,StructureParams,PhaseScreen,SF)
%USAGE:    [psi,IStatSum,phase0]=PropCode2D(k,x,xProp,y,z,thetap,phip,StructureParams,PhaseScreen,SF)
%
%INPUT:        k         =wavenumber
%                    x         =propagation distance in layer
%                    xProp =Free space propagation steps 
%                    y         =distance array 
%                    z         =distance array
%                    thetap=propagation angle (rad)
%                    phip   =propagation azimuth (rad) from xy plane
%                    StructureParams=Structure with fiellds Cp,gnu1,q0, gnu2
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
nz=length(z); dz=diff(z(1:2));

%Spectral domain sampling
dky=2*pi/ny/dy; ky=(-ny/2:ny/2-1)*dky; 
dkz=2*pi/nz/dz; kz=(-nz/2:nz/2-1)*dkz; 
c_thetap=cos(thetap); s_thetap=sin(thetap);
c_phip=cos(phip); s_phip=sin(phip);
t_thetap=tan(thetap);

%Spectral domain propagator (nz,ny)
[Ky,Kz]=meshgrid(ky,kz);
ahatk=[c_thetap; s_thetap*c_phip; s_thetap*s_phip];
kx=k*(c_thetap-sqrt(ones(nz,ny)-((Ky+ahatk(2))/k).^2-((Kz+ahatk(3))/k).^2))...
    -Ky*(t_thetap*ahatk(2))-Kz*(t_thetap*ahatk(3));
Pfac=fftshift(exp(1i*kx));

%Split-step FPE Integration
fprintf('Split step FPE Integration %5i steps \n',nx)
if PhaseScreen==1
    fprintf('Phase Screen ON \n')
end

IStatSum=cell(1,nx+nPropSteps);
fracMom=ones(1,5);
IStatSum{1}=fracMom;

psi=ones(size(Ky));
fprintf('\nFPE Integration %5i Steps \n', nx)
for nstep=1:nx
    showprogress(nstep,1)
    if PhaseScreen~=1
        if nstep>1
            %Propagation x(nstep-1) to x(nstep)
            psi_hat=fft2(psi);
            psi_hat=psi_hat.*Pfac.^(x(nx)-x(nx-1));
            psi=ifft2(psi_hat);
            fracMom=generate_fracMom(abs(psi).^2);
            IStatSum{nstep}=fracMom;
        end
    end
    phase=generatePhasePerturbation2D(nstep,Ky,Kz,SF,thetap,StructureParams);
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
psi0=psi;
for nProp=1:nPropSteps-1;
    showprogress(nProp,1)
    %Free space propagation
    psi_hat=fft2(psi);
    psi_hat=psi_hat.*Pfac.^(xProp(nProp+1)-xProp(nProp));
    psi=ifft2(psi_hat);
    fracMom=generate_fracMom(abs(psi).^2);
    IStatSum{nx+nProp-1}=fracMom;
end
fprintf('\n')
return