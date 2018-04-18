%GPS Visibility
clear
close all

dtr=pi/180;
nsec_day=24*60*60;  %seconds/day
ElevationMask=20;

origin=input('Input [latitude (rad), longitude(rad), height (m)] , CR for default ');
if isempty(origin)
    origin_llh=[37.4529598*dtr,-122.18172520000002*dtr, 22];  %Menlo Park
    LocID='Menlo Park';
else
    origin_llh=origin;
    LocID=input('Location ID string ');
end
datetime_start=datetime;
fprintf('Date %s \n',datetime_start)
[year,month,day,hh,mm,ss]=datevec(datetime_start);
%vector time=[year, monty, day, hour, minute, second+fraction]

day_of_year=datenum([year,month,day])-datenum([year,0,0]);
%Get ephemeris data for day_of_year
eph=extractRINEXeph(year,day_of_year,[]);   %Last input argument is prn []=> all prns
%eph(21,:) is toe for each 21-element ephemeris
%eph(21,1) is toe at start of day/month/year for ephemeris file

UT2eph=eph(21,1);
UT_sec=hh*60*60+mm*60+ss;   %Start time
%GPS time for datetime
[GPSTime_start,GPSweek]=UT2GPStime([year,month,day,hh,mm,ss]);
t_sec=(0:nsec_day-1); nsamps=length(t_sec);
GPSTime=GPSTime_start+t_sec;

sat_Elevation=NaN(32,nsamps);
sat_Azimuth=NaN(32,nsamps);

sat_ecf=zeros(3,nsamps);
for prn=1:32
    neph=find(eph(1,:)==prn);       %neph => eph for satellite prn
    [latency,neph_min]=min(abs(eph(21,neph)-UT2eph-UT_sec));
    fprintf('prn=%5i latency=%12.2f %5i \n',prn,latency,neph_min)
    for nsamp=1:nsamps              %Calculate ecf position & velocity
        [ sat_ecf(:,nsamp),~]=satposvel(GPSTime(nsamp),eph(:,neph(neph_min)));
    end
    sat_llh=ecf2llhT(sat_ecf);
    sat_tcs=llh2tcsT(sat_llh,origin_llh);
    ELEV=atan2(sat_tcs(3,:),sqrt(sat_tcs(1,:).^2+sat_tcs(2,:).^2));
    AZ=atan2(sat_tcs(1,:),sat_tcs(2,:)); 
    nvis=find(ELEV>ElevationMask*dtr);
    sat_Azimuth(prn,nvis)=AZ(nvis); 
    sat_Elevation(prn,nvis)=ELEV(nvis);
end

figure
imagesc(t_sec/60/60,[1:32],sat_Elevation/dtr)
axis xy
colorbar
ylabel('PRN')
xlabel('UT-hrs')
DATETIME=datetimeID(year,month,day,hh,mm)
title(['                          GPS Visibility   ',LocID,'     ',DATETIME,'            Elev-deg'])

VisabilityPlot(sat_Elevation,sat_Azimuth,1440,[LocID,'   ',DATETIME])