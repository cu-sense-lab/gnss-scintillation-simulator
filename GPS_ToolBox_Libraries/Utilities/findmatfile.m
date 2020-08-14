function nfile=findmatfile(data_files,filename)
%
nfiles=length(data_files);
for nfile=1:nfiles
    if contains(data_files(nfile).name,filename)
        break
    end
end
    