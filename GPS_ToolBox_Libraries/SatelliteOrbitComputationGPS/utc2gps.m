function [GPS_week, GPS_sec, GPS_day] = utc2gps(UTC_time,leap_sec) 
 
% [GPS_week, GPS_sec, GPS_day] = utc2gps(UTC_time,leap_sec); 
%                          or 
% GPS_time = utc2gps(UTC_time,leap_sec);  
%   
% where GPS_time = [GPS_week GPS_sec] (nx2) 
% 
% Converts a UTC time matrix to the equivalent time in GPS weeks, 
% GPS seconds, and GPS days (GPS time started at 00:00:00 6 JAN 1980) 
% 
% Input:   
%   UTC_time - matrix of the form [year month day hour minute second] (nx6) 
%               with 4-digit (1980) or 2-digit (80) years, 
%               valid years are 1980 - 2079 (2-digit 80-79) 
%   leap_sec - leap seconds applied to UTC relative to GPS (optional) 
%               can be a 1x1 or an nx1, if not entered the function will 
%               use a look-up table to determine the number of leap seconds 
% Output:  
%   GPS_week - GPS week (nx1) (if 0 or 1 output parameters are used, 
%               this is filled with [GPS_week GPS_sec] (nx2).  See the  
%               alternative calling option from above. 
%               GPS week values are kept in linear time without the 1024 rollovers 
%               (e.g. UTC time [1999 10 1 0 0 0] = GPS time [1029 432013]) 
%               use the MATLAB "mod" function to obtain 1024 based weeks 
%                  new_GPS_week = mod(GPS_week,1024); 
%   GPS_sec  - seconds into the week measured from Sat/Sun midnight (nx1) 
%   GPS_day  - days since the beginning of GPS time (optional) (nx1) 
% 
% See also GPS2UTC, GPS2LEAP 
 
% Written by: Maria Evans/Jimmy LaMance 10/9/96 
% Copyright (c) 1998 by Constell, Inc. 
 
% functions called: ERR_CHK, UTC2LEAP 
 
%%%%% BEGIN VARIABLE CHECKING CODE %%%%% 
% declare the global debug variable 
global DEBUG_MODE 
 
% Initialize the output variables 
GPS_week=[]; GPS_sec=[]; GPS_day=[]; 
 
% Check the number of input arguments and issues a message if invalid 
msg = nargchk(1,2,nargin); 
if ~isempty(msg) 
  fprintf('%s  See help on UTC2GPS for details.\n',msg); 
  fprintf('Returning with empty outputs.\n\n'); 
  return 
end 
 
% Check for 2-digit years 
I_2 = find(UTC_time(:,1) < 100); 
if ~isempty(I_2)    
  I_1900 = find(UTC_time(I_2,1) >= 80); 
  I_2000 = find(UTC_time(I_2,1) < 80); 
  if ~isempty(I_1900) 
    UTC_time(I_2(I_1900),1) = UTC_time(I_2(I_1900),1) + 1900; 
  end  
  if ~isempty(I_2000) 
    UTC_time(I_2(I_2000),1) = UTC_time(I_2(I_2000),1) + 2000; 
  end                                                      
end 
   
if nargin < 2 & size(UTC_time,2) == 6 
  leap_sec = utc2leap(UTC_time);  
elseif nargin < 2 & size(UTC_time,2) ~= 6  
  leap_sec = ones(size(UTC_time,1),1) * leap_sec; 
end % if              
 
% Get the current Matlab version 
matlab_version = version; 
matlab_version = str2num(matlab_version(1)); 
 
% If the Matlab version is 5.x and the DEBUG_MODE flag is not set 
% then set up the error checking structure and call the error routine. 

 
%%%%% END VARIABLE CHECKING CODE %%%%% 
 
%%%%% BEGIN ALGORITHM CODE %%%%% 
% data matrix with the number of days per month              
% days in full months for leap year 
leapdays =   [0 31 60 91 121 152 182 213 244 274 305 335];   
% days in full months for standard year  
noleapdays = [0 31 59 90 120 151 181 212 243 273 304 334];                                                      
 
% Leap year flag   
% determine which input years are leap years 
leap_year = ~rem((UTC_time(:,1)-1980),4);      
I_leap = find(leap_year == 1);                % find leap years 
I_no_leap = find(leap_year == 0);             % find standard years 
 
% Round the first column of the almanac to get integer satellite numbers 
a_round = ceil(UTC_time(:,2)); 
d_diff = a_round - UTC_time(:,2); 
I_diff_non_zero = find(d_diff ~= 0); 
if ~isempty(I_diff_non_zero) 
  UTC_time(I_diff_non_zero,2) = a_round(I_diff_non_zero); 
  fprintf('Non-integer months numbers have been rounded up in UTC2GPS.\n') 
end % if ~isempty(I_diff_non_zero) 
 
% generate a matrix that has the days per completed month for both  
% leap and standard years 
if any(I_leap) 
  dayspermonth(I_leap) = leapdays(UTC_time(I_leap,2)); 
end % if any(I_leap) 
 
if any(I_no_leap) 
  dayspermonth(I_no_leap) = noleapdays(UTC_time(I_no_leap,2)); 
end % if any(I_no_leap) 
 
% compute the number of leap days encounted in past years,  
% need to add one to the fix computation to get the year correct                                                          
leapyrs = fix((UTC_time(:,1) - 1980) ./ 4) + eval('~leap_year'); 
 
% Compute the number of days in completed years since 1980 
gpsday = (UTC_time(:,1) - 1980) .* 365 + leapyrs; 
                                
% Add the number of days for each completed month 
gpsday = gpsday + dayspermonth'; 
 
% add the number of days for each completed day  
% (which is 1 less than the current day) 
gpsday = gpsday + (UTC_time(:,3) - 1); 
 
% add the fraction of days for each completed hour (hour/24) 
% seconds into the current day (avoids round-off errors) 
gpssec = UTC_time(:,4) .* 3600;  
 
% add the fraction of days for each completed minute (minute/1440) 
gpssec = gpssec + UTC_time(:,5) .* 60; 
 
% add the fraction of days each completed second (second/86400) 
gpssec = gpssec + UTC_time(:,6); 
 
% Subtract 5 days because the starting date is 00:00 6 January 1980 
gpsday = gpsday - 5;  
 
GPS_day = gpsday;         
 
% formal return parameter GPS week 
GPS_week = fix(GPS_day ./ 7);  
days_into_week = floor(rem(GPS_day,7)); 
 
% Add leap seconds to the gps time  
GPS_sec = days_into_week' * 86400 + gpssec' + leap_sec'; 
 
% check to make sure the leap seconds don't force a week rollover 
I_next = find(GPS_sec >= 86400 * 7); 
if any(I_next) 
  GPS_week(I_next) = GPS_week(I_next) + 1; 
  GPS_sec(I_next) = GPS_sec(I_next) - 86400 * 7; 
end % if 
 
% check the output arguments, if the out requested is length 0 ot 1,  
% return GPS week and GPS seconds in a single vector 
if nargout == 0 | nargout == 1 
  GPS_week = [GPS_week GPS_sec']; 
end % if nargout == 0 | nargout == 1 
 
% Return the GPS week, GPS seconds into the week, and the number of  
% completed days (GPS). formal return parameter GPS days since Jan 6 1980 
if nargout == 3 
  GPS_day = GPS_week * 7 + GPS_sec / 86400;    
end % if nargout == 3                      
 
%%%%% END ALGORITHM CODE %%%%% 
 
% end UTC2GPS                    

