%% Add the path to the libraries to current search path
path2directories=pwd;
path2Libraries=[path2directories,'/Libraries'];
addpath([path2Libraries,'/PropGeomCalc']);
addpath([path2Libraries,'/GenScintFieldRealization'])
addpath([path2Libraries,'/GPS_CoordinateXforms']);
addpath([path2Libraries,'/IGRF_Compston']);
addpath([path2Libraries,'/Utilities']);

%% User input parameters***************************************************
%Please specify date and time as [year month day hour minute second]
userInput.dateTime = [2014 01 02 10 00 00]; 

%Please choose data length for simulation\n (300s, 600s, or 900s)
userInput.length = 300;

% Please specify receiver position as [lat(rad), lon(rad), height(m)\n]'
userInput.RXPos = [0.3876 1.9942 59.6780]';

% Please input receiver Vel as [V1,V2,V3]'
%    V1 = east-west velocity on the earth arc (m/s, eastward +)
%    V2 = north-south velocity on the earch arc (m/s, northward +)
%    V3 = up-down velocity (m/s, up +)
userInput.RXVel = [100 0 0]';

% Please specify satellite PRN (0~32)
userInput.PRN = 1;

% Plotting figures of the simulated propagation geometry and scintillation intensity and phase? yes-1/no-0
userInput.plotSign = 1;

% Please specify how many GPS frequencies to simulate (1- GPS L1 only; 2 - GPS L1 and L2; 3 - GPS L1,L2, and L5)
userInput.frequencyNo = 3;

% Please specify the S4 index and tau0 for the ground observed scintillation 
userInput.S4 = 0.8; % S4 index (0~1)
userInput.tau0 = 0.7; % Signal intensity decorrelation time in sec.

%% Obtain the U and rhoVeff values based on the user input S4 and tau. 
[U_mapped,rhoFVeff_mapped] = ParaMapping(userInput);

%% Calculate the propagation geometry**************************************
if(sum(userInput.RXVel)~=0)
    satGEOM = RunPropGeomCalc(userInput,rhoFVeff_mapped);
end

%% Generate scintillation signal field realizations
[Scin_psi, Scin_amp, Scin_phi] = RunGenScintFieldRealization(userInput,satGEOM,U_mapped,rhoFVeff_mapped);