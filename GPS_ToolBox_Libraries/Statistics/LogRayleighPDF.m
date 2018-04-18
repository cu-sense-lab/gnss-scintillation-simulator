y_counts=linspace(0,1,500);
y_dB=40*linspace(0,1,500)-60;
dy=diff(y_dB(1:2));
N=10^(-4);
dB2N=log(10)/10;
PDF=exp(-exp(dB2N*y_dB)/N).*exp(dB2N*y_dB)*dB2N/N;
Norm_Check=sum(PDF)*dy;
figure
plot(y_dB,log(max(PDF,exp(-20))),'r')
grid on
%axis([min(y) max(y) -20 0])
ymax=log(N)/dB2N;
hold on
plot(ymax,-10,'m^')

figure
plot(y_counts,log(PDF),'b')
grid on