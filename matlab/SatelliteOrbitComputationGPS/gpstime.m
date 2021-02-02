% This function converts between UTC and GPS time. An input of [year month
% day hour minute second] in UTC will yeild an output of [wkno tow] GPS
% time. An input of either [wkno tow] or [gpsseconds] GPStime will yeild an
% output of [year month day hour minute second] UTC time.
% 
% @author Steve Taylor


function [out]=gpstime(in)

stepdates = [...
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
    'Jul 1 2012'
    'Jul 1 2015'
    'Jan 1 2017'];

stepdates=datenum(stepdates);

startdate=datenum(1980,01,06,0,0,0);

if(length(in)==6)
    indate=datenum(in);
    gpsdays=indate-startdate;
    gpssecs=gpsdays*(24*3600)+utcoffset(stepdates,indate);
    out(1)=floor(gpssecs/(168*3600));
    out(2)=gpssecs-out(1)*3600*168;
else
    out=0;
    if(length(in)==2)
        t=in(1)*168*3600+in(2);
        clear in;
        in=t;
    end
    if(length(in)==1)
        day=in/(3600*24);
        offset=utcoffset(stepdates,day+startdate);
        day=(in-offset)/(3600*24)+startdate;
        st=datestr(day,'dd-mm-yyyy HH:MM:SS');
        %[out(1) out(2)]=sscanf(st,'%d-%d-%d %d:%d:%d');
        out(3)=str2num(st(1:2));
        out(2)=str2num(st(4:5));
        out(1)=str2num(st(7:10));
        out(4)=str2num(st(12:13));
        out(5)=str2num(st(15:16));
        out(6)=str2num(st(18:19));
    end
end
end

function offset=utcoffset(stepdates,indate)
i=1;c=0;
while(i<length(stepdates))
    if(indate >= stepdates(i))
        c=c+1;
    end
    i=i+1;
end
offset=c;
end
