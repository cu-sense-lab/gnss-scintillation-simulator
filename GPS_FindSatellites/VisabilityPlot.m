function VisabilityPlot(Elevation,Azimuth,nplt,plotID)
%
[nsets,nsamps]=size(Elevation);
ndec=ceil(nsamps/nplt);
figure
for nn=1:nsets
    rho=cos(Elevation(nn,1:ndec:nsamps));
    x=cos(Azimuth(nn,1:ndec:nsamps)).*rho;
    y= sin(Azimuth(nn,1:ndec:nsamps)).*rho;
    hold on
    plot(x,y,'r.')
    [~,np]=min(x);
    text(x(np)-0.075,y(np),num2str(nn))
end
grid on
axis equal
title(plotID)
bold_fig
xlabel('West-East')
ylabel('North-South')
