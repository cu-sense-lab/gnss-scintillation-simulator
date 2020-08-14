%add Librarier to MatLab path
%Assumes \GPS_Libraries is a subdirectory in \mfiles
%
path2directories=fileparts(pwd);
nmfiles=strfind(path2directories,'mfiles');
if length(path2directories)>(nmfiles+5)
    path2directories=path2directories(1:nmfiles+6);
else
    path2directories=[path2directories,'\'];
end
mfile_dirs=dir(path2directories);
for ndir=1:length(mfile_dirs)
    if contains(mfile_dirs(ndir).name,'GPS_Libraries')
        fprintf('GPS Library Date%s \n',datestr(mfile_dirs(ndir).datenum))
    end
end
path2Libraries=[path2directories,'GPS_Libraries'];
addpath([path2Libraries,'\Ispectrum']);
setGlobalpath2exe([path2Libraries,'\Ispectrum']);
fprintf('%s \n',getGlobalpath2exe)

addpath([path2Libraries,'\IPE_Utilities'])
addpath([path2Libraries,'\Statistics'])
addpath([path2Libraries,'\Utilities']);
addpath([path2Libraries,'\PropCodes'])

addpath([path2Libraries,'\GPS_CommonAnalysis']);
addpath([path2Libraries,'\SatelliteOrbitComputationGPS']);
addpath([path2Libraries,'\GPS_Transforms\GPS_CoordinateXforms']);
addpath([path2Libraries,'\IGRF_Compston']);
addpath([path2Libraries,'\Scintillation'])
addpath([path2Libraries,'\GPS_Processing'])
addpath([path2Libraries,'\ConfigurationSpace'])
addpath([path2Libraries,'\UnwrapPhase'])

%add path to Wavelab850 subset
path2wavelets=fullfile(path2Libraries,'\WaveletScaleSpectra');
addpath([path2wavelets,'\WaveletScaleSpectraUtilities']);
addpath(genpath([path2wavelets,'\WaveletsSubset\']))
 