function   path2library=findpathLibraries(Library);
%
path2directories=fileparts(pwd);
nmfiles=strfind(path2directories,'mfiles');
if length(path2directories)>(nmfiles+5)
    %Active directory is subdirectory
    path2directories=path2directories(1:nmfiles+6);
else
    path2directories=[path2directories,'\'];
end
path2Libraries=[path2directories,'GPS_Libraries'];
path2library=[path2Libraries,'\',Library];
return
