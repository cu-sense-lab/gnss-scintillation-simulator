function [GAST]=siderealtime(timehack,timezone)

% This Sidereal Clock uses your computer clock's definition of Coordinated Universal Time (UTC)
% to compute the corresponding Mean Sidereal Time at the Greenwich meridian of 0°: Longitude.

% * Inputs
% * timehack [year month day hour min sec]

% * Output
% * Greenwich Apparent Sidereal Time


if nargin < 1
    timehack=clock;
    timezone = 6;  %UT = CST + 6
end;

if isscalar(timehack) == 1
    timezone = timehack;
    timehack = clock;
else
    timezone = 6;
end;

% User must set the default timezone


% Convert local time to UT
timehack(1,4)=timehack(1,4) + timezone;

% Calculate JD (Julian Date) at 0:0:0 UT on date of interest
JD=367*timehack(1,1)-floor(7*(timehack(1,1)+floor(timehack(1,2)+9)/12)/4)+floor(275*timehack(1,2)/9)...
    +timehack(1,3)+1721013.5;
D = JD - 2451545;

% Calculate UT in decimal hours
UT=timehack(1,4)+(timehack(1,5)+timehack(1,6)/60)/60;
% Calculate T (number of centuries since UT began)
T=(JD-2451545)/36525;
% Calculate T0
T0 = 6.697374558 + (2400.051336*T) + (0.000025862*T^2) + (UT*1.0027379093);

% Reduce T0 to a value between 0 and 24 by adding or subtracting multiples of 24
T0=mod(T0,24);
% Greenwich mean sidereal time in decimal hours
GMST=T0;

% Solve for Apparent Sidereal time by solving for the
% nutation in right ascension
Om=125.04-0.052954*D;                              %Longitude of the ascending node of the Moon
L=280.47+0.98565*D;                                %Mean Longitude of the Sun
eps=23.4393-0.0000004*D;                           %obliquity
delta_psi=-0.000319*sind(Om)-0.000024*sind(2*L);   %nutation in longitude
eqeq=delta_psi*cosd(eps);                          %equation of the equinoxes
GAST=GMST + eqeq;                                  %Greenwich Apparent Sidereal Time

% Convert from decimal hours to degs to rads to match input for conversion
% to/from Earth Centered Inertial/Earth Centered Rotational functions
GAST=GAST*360/24*pi/180;

% % Verification code:  LAST_deg:LAST_min:LAST_sec should compare to clocks
% % Convert from radians to hour degrees
% thetaGr=GAST*180/pi*24/360;
% % Subtract longitude of interest
% longitude =  0;   %west longitudes are negative
% LAST=thetaGr+(longitude/15);
% % Reduce LST to a value between 0 and 24 by adding or subtracting multiples of 24.
% LAST=mod(LAST,24);
% %Local Sidereal Time in deg,min,sec
% LAST_deg=floor(LAST);
% LAST_min=floor((LAST-floor(LAST))*60);
% LAST_sec=(((LAST-floor(LAST))*60)-floor((LAST-floor(LAST))*60))*60;
% test_out=strcat(num2str(LAST_deg),':',num2str(LAST_min),':',num2str(LAST_sec));
% clc;
% disp(test_out);

