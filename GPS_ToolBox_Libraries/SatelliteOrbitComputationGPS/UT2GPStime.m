function               [GPStime_sec,GPSweek,GPSweek_z,leapsec]=UT2GPStime(UTdate_time_vector)
%USAGE           [GPStime_sec,GPSweek,GPSweek_z,leapsec]=UT2GPStime(UTdate_time_vector) 
%
%Modification of Matlab Central script utc2gps
%   Copyright 2008 Ian M. Howat, ihowat@gmail.com
%   $Version: 1.0 $  $Date: 23-Aug-2008 13:56:44 $
% by Charles Rino
%Rino Consulting

sec_day=24*60*60;
if length(UTdate_time_vector)~=6
    error('UTdate_time_vector error')
end
daynumber=datenum(UTdate_time_vector);
GPSepoch=datenum([1980,1,6,0,0,0]);
GPStime=(daynumber-GPSepoch);
GPSweek=floor(GPStime/7);
GPStime_sec=(GPStime-GPSweek*7)*sec_day;
GPSweek_z=mod(GPSweek,1024);

% ADD NEW LEAP DATES HERE:
%FROM 
stepdates = [...
    'Jan 6 1980'
    'Jul 1 1981'
    'Jul 1 1982'
    'Jul 1 1983'
    'Jul 1 1985'
    'Jan 1 1988'
    'Jan 1 1990'
    'Jan 1 1991'
    'Jul 1 1992'
    'Jul 1 1993'
    'Jul 1 1994'
    'Jan 1 1996'
    'Jul 1 1997'
    'Jan 1 1999'
    'Jan 1 2006'
    'Jan 1 2009'
    'Jan 1 2012'
    'Jan 1 2015'
    'Jan 1 2017'];

% Convert Steps to datenums and make step offsets
stepdates = datenum(stepdates)'; %step date coversion
steptime = (0:length(stepdates)-1)'./86400; %corresponding step time (sec)

if ~isempty(find(daynumber < stepdates(1), 1))%date0 must all be after GPS start date
    error('Input dates must be after 00:00:00 on Jan 6th 1980') 
end
leapsec=steptime(sum((daynumber - stepdates) >= 0,2));
return