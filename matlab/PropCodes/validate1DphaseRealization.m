function     [log10q,SDFdB]=validate1DphaseRealization(y,SF,thetap,StructureParams,PlotID)

Cp=StructureParams.Cp;
p1=StructureParams.p1;
q0=StructureParams.q0;
p2=StructureParams.p2;

ny=length(y); dy=diff(y(1:2));
%Spectral domain sampling

dky=2*pi/ny/dy; ky=(-ny/2:ny/2-1)*dky;
rootSDF=root_phaseSDF1(Cp,p1,p2,q0,ky,SF*sec(thetap)^2);

[phase]=generatePhasePerturbation1D(1,ky,SF,thetap,StructureParams);

ny=length(ky);
PSD=abs(fftshift(fft(phase))).^2/ny;
Prms=sqrt(sum(rootSDF.^2)*dky/2/pi);
PrmsID=sprintf('Prms=%5.2f, rms=%5.2f ',Prms,rms(phase));
log10q=log10(ky(ny/2+2:ny)/2/pi);
SDFdB=dB10( rootSDF(ny/2+2:ny).^2);
figure
subplot(2,1,1)
plot(y/1000,phase,'r')
grid on
ylabel('phase-rad')
xlabel('y-km')
title(PlotID)
subplot(2,1,2)
plot(log10q,dB10(PSD(ny/2+2:ny)),'r')
grid on
hold on
plot(log10q,SDFdB,'b')
grid on
hold on
plot(log10(q0/2/pi),dB10(Cp*q0^(-p1)),'mp')
legend('PhasePSD','SDF')
xlabel('log10(q/(2*pi)')
ylabel('PSD-dB')
title(PrmsID)

return