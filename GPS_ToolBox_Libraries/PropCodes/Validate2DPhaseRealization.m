function Validate2DPhaseRealization(y,z,SF,thetap,StructureParams,PlotID)

Cp=StructureParams.Cp;
gnu1=StructureParams.gnu1;
gnu2=StructureParams.gnu2;
q0=StructureParams.q0;

a=StructureParams.a;
b=StructureParams.b;
A=StructureParams.A;
B=StructureParams.B;
C=StructureParams.C;


ny=length(y); dy=diff(y(1:2));
nz=length(z); dz=diff(z(1:2));
%Spectral domain sampling

dky=2*pi/ny/dy; ky=(-ny/2:ny/2-1)*dky;
dkz=2*pi/nz/dz; kz=(-nz/2:nz/2-1)*dkz;
[Ky,Kz]=meshgrid(ky,kz);
%dkyC=diff(Ky(1,1:2));dkzC=diff(Kz(1:2,1));

 [rootSDF]=root_phaseSDF2(Cp,gnu1,gnu2,q0,Ky,Kz,SF,A,B,C);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 phase=generatePhasePerturbation2D(1,Ky,Kz,SF,thetap,StructureParams);
 
Prms=sqrt(sum(rootSDF(:).^2)*dky/2/pi*dkz/2/pi);
PrmsID=sprintf('Prms=%5.2f, rms=%5.2f ',Prms,rms(phase(:)));

figure
imagesc(y/1000,z/1000,phase)
title(PlotID)
xlabel('y-km')
zlabel('z-km')
colorbar

phase_y=phase(nz/2+1,:);
PSDy=abs(fftshift(fft(phase_y))).^2/ny;
phase_z=phase(:,ny/2+1:ny);
PSDz=abs(fftshift(fft(phase_z))).^2/nz;

figure
subplot(2,1,1)
xx=log10(ky(ny/2+2:ny)/2/pi);
yy=dB10(PSDy(ny/2+2:ny));
plot(xx,yy,'r')
grid on
xlabel('log10(ky/(2\pi))')
title(PrmsID)
subplot(2,1,2)
xx=log10(kz(nz/2+2:nz)/2/pi);
yy=dB10(PSDz(nz/2+2:nz));
plot(xx,yy,'r')
grid on
xlabel('log10(kz/(2\pi))')

return