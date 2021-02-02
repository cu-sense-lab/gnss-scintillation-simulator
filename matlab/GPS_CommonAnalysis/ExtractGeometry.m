function          [satGEOM]=ExtractGeometry(mfile,StationID,h_intercept,Drift)
%USEAGE:    [satGEOM]=ExtractGeometry(mfile,StationID,h_intercept,Drift)
%

c = 299792458;        %Speed of light (vacuum)

load(mfile)
[~,mfile_name]=fileparts(mfile);

%GPSTime=LOG.GPSTime;

%Time derived from mfile_name
[UTC_start,~,year,~,~]=extractUTC_time(mfile_name);
t_sec0=UTC_start(4)*60*60+UTC_start(5)*60+UTC_start(6);

UTC_start_day=datenum(UTC_start);
day_of_year=fix(UTC_start_day-datenum(year,0,0,0,0,0));

%Generate 1-s time samples for Geometry
nsampsG=length(GPSTime(1,:));
t_secG=t_sec0+(0:nsampsG-1);

%Extract veff and interpolate to t_sec50 for time to space conversion

%t_secG is time in seconds that includes t_sec (=> add one extra sec)
nsampsG=nsampsG+1;
t_secG(nsampsG)=t_secG(nsampsG-1)+1;

GSPTime(1,nsampsG)=GPSTime(1,nsampsG-1);
GPSTime(2,nsampsG)=GPSTime(2,nsampsG-1)+1;
satGEOM=GenerateGPSGeometry(t_secG,GPSTime,UTC_start,eph,origin_llh,StationID,h_intercept,Drift);

return
