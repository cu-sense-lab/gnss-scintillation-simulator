function [UTC_start,UTC_end,year,month,day,mfileID]=extractUTC_time(mfilename)
%Extract UTC_time from CSU format filename
%
%Jan 2016 Version
Ubars        =strfind(mfilename,'_');
mfileID      =mfilename(1:Ubars(1)-1);
date_str    =mfilename(Ubars(1)+1:(Ubars(2)-1));
year          =str2double(date_str(1:4));
month       =str2double(date_str(5:6));
day            =str2double(date_str(7:8));
hours0      =str2double(mfilename(Ubars(3)+1:Ubars(3)+2));
min0         =str2double(mfilename(Ubars(3)+3:Ubars(3)+4));
sec0         =str2double(mfilename(Ubars(3)+5:Ubars(3)+6));
hours1      =str2double(mfilename(Ubars(4)+1:Ubars(4)+2));
min1         =str2double(mfilename(Ubars(4)+3:Ubars(4)+4));
sec1         =str2double(mfilename(Ubars(4)+5:Ubars(4)+6));
UTC_start=[year,month,day,hours0,min0,sec0];
UTC_end =[year,month,day,hours1,min1,sec1];
return