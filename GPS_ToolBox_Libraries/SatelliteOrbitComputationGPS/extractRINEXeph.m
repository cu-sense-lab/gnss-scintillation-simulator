function       eph=extractRINEXeph(year,day_of_year,PRN)
%USAGE:  eph=extractRINEXeph(year,day_of_year,PRN)
%Extract RINEX ephemeris file from website:
%ftp://cddis.gsfc.nasa.gov/gnss/data/daily/YYYY/DDD/YYn/brdcDDD0.YYn.Z


datadir='ftp://cddis.gsfc.nasa.gov/gnss/data/daily/';
YYYY=num2str(year);
DDD=num2str(day_of_year);
if length(DDD)==1
    DDD=['0','0',DDD];
elseif length(DDD)==2
    DDD=['0',DDD];
end
YY=YYYY(3:4);
RinexFile=[datadir,YYYY,'/',DDD,'/',YY,'n/brdc',DDD,'0.',YY,'n.Z'];
[~,ephfile]=fileparts(RinexFile);
if ~exist(ephfile,'file')
    if ~exist([ephfile,'.Z'],'file')
        gunzip(RinexFile,pwd)
    end
    fprintf('Unzip %s   then type dbcont \n',ephfile)
    keyboard
else
   fprintf('Using ephemeris file %s \n',ephfile)
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eph,~] =rinexe(ephfile); %Extract 21XN ephemeris files
if ~isempty(PRN)
    %eph(21,:) is toe for each 21-element ephemeris
    neph=find(eph(1,:)==PRN);
    if isempty(neph)
        error('RINEX error ')
    else
        eph=eph(:,neph);
    end
end
return