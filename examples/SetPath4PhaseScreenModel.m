%Set Path for GPS analysis and simulations
path2directories=fileparts(pwd);
path2Libraries=[path2directories,'\GPS_ToolBox_Libraries'];
fprintf('*****USING GPS_ToolBox_Libraries********* \n')

%
addpath([path2Libraries,'\Ispectrum']);
setGlobalpath2exe([path2Libraries,'\Ispectrum\']);  %Make Ispectrum Global
addpath([path2Libraries,'\Ispectrum_Utilities\'])

addpath([path2Libraries,'\GPS_CommonAnalysis']);
addpath([path2Libraries,'\Scintillation']);
addpath([path2Libraries,'\PropCodes']);
addpath([path2Libraries,'\Statistics']);
addpath([path2Libraries,'\Utilities']);

%add paths to orbit computation and GPS_coordinate transformations
addpath([path2Libraries,'\SatelliteOrbitComputationGPS']);
addpath([path2Libraries,'\GPS_Transforms\GPS_CoordinateXforms']);
addpath([path2Libraries,'\IGRF_Compston']);

%add path to Wavelab850 subset
path2wavelets=fullfile(path2Libraries,'\WaveletScaleSpectra');
addpath([path2wavelets,'\WaveletScaleSpectraUtilities']);
addpath(genpath([path2Libraries,'\waveletScaleSpectra\']))