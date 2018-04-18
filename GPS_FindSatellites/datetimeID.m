function     DATETIME=datetimeID(year,month,day,hh,mm)
%
if hh<10
    HH=['0',num2str(hh)];
else
    HH=num2str(hh);
end
HH=[HH,':'];
if mm<10
    MM=['0',num2str(mm)];
else
    MM=num2str(mm);
end
HHMM=[HH,MM];        
DATETIME=[num2str(month),'/',num2str(day),'/',num2str(year),' ',HHMM];
return