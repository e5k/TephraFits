% This is a demo run for TephraFits
% For help about the function, type:
% >> help tephraFits

%% Example 0: Plot the data
clear;
thickness	= [100, 50, 30, 20, 10, 5]; % Isopach thickness (cm)
areaT       = [7.0, 8.9, 12.3, 17.4, 21.34, 25.4]; % Square-root of area (km)
tephraFits(areaT, thickness);

%% Example 1: Isopach
clear;
thickness	= [100, 50, 30, 20, 10, 5]; % Isopach thickness (cm)
areaT       = [7.0, 8.9, 12.3, 17.4, 21.34, 25.4]; % Square-root of area (km)
isopach  	= tephraFits(areaT, thickness, {'exponential', 'powerlaw', 'weibull'}, 'BIS', 2, 'C', 150); % Poor weibull fit
isopach  	= tephraFits(areaT, thickness, {'exponential', 'powerlaw', 'weibull'}, 'BIS', 2, 'C', 150, 'lambdaRange', [.01 100], 'nRange', [.01 100] ); % Refined weibull fit

%% Example 2: Isomass
clear;
massAcc     = [100, 50, 30, 20, 10, 5] ./ 1e2 .* 1000; % Isomass accumulation (kg/m2)
areaM       = [7.0, 8.9, 12.3, 17.4, 21.34, 25.4]; % Square-root of area (km)
isomass     = tephraFits(areaM, massAcc, {'exponential', 'powerlaw', 'weibull'}, 'dataType', 'isomass', 'BIS', 2, 'C', 150, 'lambdaRange', [.1 100], 'nRange', [.1 100]); 

%% Example 3: Isopleth
clear;
diameter	= [0.8, 1.0, 1.6, 2.0, 3.0, 3.2, 4.0]; % Isopleth diameter (cm)
areaD		= [23.7, 22.9, 20.4, 18.3, 14.9, 14.1, 12.2]; % Square-root of area (km)
isopleth  	= tephraFits(areaD, diameter, {'exponential', 'weibull'}, 'dataType', 'isopleth', 'lambdaRange', [.1 100], 'nRange', [.1 100]);

%% Example 4: Thickness transects
clear;
% Downwind transect
thicknessDW = [35, 30, 32, 18, 23, 8, 6]; % Outcrop thickness (cm)
distanceDW  = [13.1, 13.4, 14.2, 16.0, 16.3, 23.9, 24.8]; % Outcrop distance (km)
transectDW	= tephraFits(distanceDW, thicknessDW, 'exponential', 'dataType', 'transect');
% Crosswind transect - White points in Fig. 1
thicknessXW = [115, 42, 33, 23, 13]; % Outcrop thickness (cm)
distanceXW  = [4.4, 5.0, 5.3, 8.3, 8.9]; % Outcrop distance (km)
transectXW	= tephraFits(distanceXW, thicknessXW, 'exponential', 'BIS', 3, 'dataType', 'transect');

%% Example 5: Uncertainty assessment
clear;
thickness	= [100, 50, 30, 20, 10, 5]; % Isopach thickness (cm)
area 		= [7.0, 8.9, 12.3, 17.4, 21.34, 25.4]; % Square-root of area (km)
CE 			= 20; % 20% error on the distal integration limit of the power law fit

% Example 5.1: uniform errors on xData and yData and uniform distribution of errors
thicknessE  = 10; % 10% error on all yData
areaE 		= 10; % 10% error on all xData
isopachP 	= tephraFits(area, thickness, {'exponential', 'powerlaw', 'weibull'}, 'BIS', 2, 'C', 200, 'lambdaRange', [.01 100], 'nRange', [.01 100], 'runMode', 'probabilistic', 'nbRuns', 100, 'errorType', 'uniform', 'xError', areaE, 'yError', thicknessE, 'CError',  CE);

% Example 5.2: specific errors for each value on xData and yData and Normal distribution of errors
thicknessE  = [10, 10, 20, 20, 30, 30]; % yError varying from 10% to 30%
areaE 		= [10, 10, 20, 20, 30, 30]; % xError varying from 10% to 30%
isopachP 	= tephraFits(area, thickness, {'exponential', 'powerlaw', 'weibull'}, 'BIS', 2, 'C', 200, 'lambdaRange', [.01 100], 'nRange', [.01 100], 'runMode', 'probabilistic', 'nbRuns', 100, 'errorType', 'normal', 'xError', areaE, 'yError', thicknessE, 'CError',  CE);

%% Example 6: Eruption classification
clear;
% Isopach data
thickness	= [100, 50, 30, 20, 10, 5]; % Isopach thickness (cm)
area 		= [7.0, 8.9, 12.3, 17.4, 21.34, 25.4]; % Square-root of area (km)
isopach  	= tephraFits(area, thickness, {'exponential', 'weibull'}, 'BIS', 3, 'lambdaRange', [.01 100], 'nRange', [.01 100], 'plotType', 'none'); % Isopach fit

% Isopleth data
diameter	= [0.8, 1.0, 1.6, 2.0, 3.0, 3.2, 4.0]; % Isopleth diameter (cm)
area 		= [23.7, 22.9, 20.4, 18.3, 14.9, 14.1, 12.2]; % Square-root of area (km)
isopleth  	= tephraFits(area, diameter, {'exponential', 'weibull'}, 'dataType', 'isopleth', 'lambdaRange', [.1 100], 'nRange', [.1 100], 'plotType', 'none'); % Isopleth fit

% Classification
tephraFits(isopach, isopleth); % Plot classification

