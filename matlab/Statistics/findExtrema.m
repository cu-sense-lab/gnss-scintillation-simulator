function         [Pts,nLocs]=findExtrema(y,Np)

%
[PtsU,nLocsU]=findpeaks(y,'SortStr','descend','MinPeakWidth',3,'Npeaks',Np/2);
[PtsL,nLocsL]=findpeaks(-y,'SortStr','descend','MinPeakWidth',3,'Npeaks',Np/2);
Pts=[PtsU,-PtsL];
nLocs=[nLocsU,nLocsL];
return