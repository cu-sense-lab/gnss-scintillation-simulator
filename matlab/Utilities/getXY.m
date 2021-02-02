function [X,Y]=getXY(nfig,npoints)
%Select npoints XY values from GUI
figure(nfig)
[X,Y]=ginput(npoints);
return