function  [mfile,mfile_dir,m]=selectmfile(wcfilt)
%List all files in wcfile
%Select mfile m in dir list
%
[mfile_dir,~]=fileparts(wcfilt);
mfiles=dir(wcfilt);
if isempty(mfiles);
    mfile=[];
    m=0;
    return
else
    nfiles=length(mfiles);
    for m=1:nfiles
        fprintf('%3i %s \n',m,mfiles(m).name);
    end
    m=input('Select file ');
    mfile=mfiles(m).name;
end
return