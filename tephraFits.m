function [V,C] = tephraFits(xData, yData, fitType, varargin)
%% TEPHRAVOLUME Calculate the volume/mass of tephra deposits
%   V = tephraVolume(xdata, ydata, method, deposit, ...) returns a structure
%   containing the volume [km3] (or mass [kg]), fit details and VEI.
%
%   Required input arguments:
%       xdata:   Square-root of isopach area (km) or distance from the vent (km)
%       ydata:   Thickness (cm) or mass load (g/m2)
%       fitType: 'exponential'   (Pyle 1989, Fierstein & Nathenson 1992)
%                'powerlaw'      (Bonadonna & Houghton 2005)
%                'weibull'       (Bonadonna & Costa 2012)
%       Multiple fits can be specified in a cell array (Ex. {'exponential','powerlaw'})
%
%   Fit-specific input arguments:
%       Exponential: 0 or 1 arguments
%           - 'BIS': 
%           - Fit thickness data with 1 exponential segment:
%               vExp = tephraVolume(thickness, area, 'exponential', 'volume') 
%
%           - Fit thickness data with 2 exponential segments, where the
%           break in slope is located between the 3rd and 4th points (in
%           decreasing thickness):
%               vExp = TEPHRAVOLUME(thickness, area, 'exponential', 'volume', 3) 
%
%           - Fit thickness data with 3 exponential segments, where the
%           break in slope are located after the 3rd and 6th points (in
%           decreasing thickness):
%               vExp = tephraVolume(thickness, area, 'exponential', 'volume', [3,6])
%
%       Power-law: 2 arguments
%           tephraVolume(thickness, area, 'powerlaw', 'volume', C, T0)
%           	C:  Distal integration limit (km)
%               T0: Intersection of the proximal exponential segment
%           Examples:
%           - Fit thickness data with a Power-Law up to 300 km and using 
%           the T0 from a 1-exponential segment approach:
%               vPL = tephraVolume(thickness, area, 'powerlaw', 'volume', 300, vExp.T0)
%
%           - Fit thickness data with a Power-Law up to 100 km and using 
%           the T0 from the proximal segment of a multiple exponential 
%           segment approach:
%               vPL = tephraVolume(thickness, area, 'powerlaw', 'volume', 100, vExp.T0(1))
%
%       Weibull: 1 or 2 arguments
%           Examples:
%           - Fit thickness data with a Weibull function using the volume
%           (km3) obtained with a power-law to constrain the optimization 
%            ranges of the lambda and n parameters:
%               vWBL = tephraVolume(thickness, area, 'weibull', 'volume', vPL.volume_km3)
%
%           - Fit thickness data with a Weibull function manually specifying
%           the optimization ranges of the lambda and n parameters:
%               vWBL = tephraVolume(thickness, area, 'weibull', 'volume', [0.1, 1000], [0.1, 1000])
%
%   Optional imput arguments passed as parameter pairs
%       'deposit':      Defines the type of deposit to fit
%                       - 'volume':     Calculates the tephra volume (km3) based on Ln(Thickness) vs. sqrt isopach area relationship (default)
%                                       xData: square-root of isopach area (km)
%                                       yData: isopach thickness (cm)
%                       - 'mass':       Calculates the tephra mass (kg) based on Ln(Mass accumulation) vs. sqrt isomass area relationship
%                                       xData: square-root of isomass area (km)
%                                       yData: isomass accumulation (kg/m2)
%                       - 'thinning':   Thinning profile based on Ln(Thickness) vs. distance relationship
%                                       xData: distance from source (km)
%                                       yData: thickness (cm)
%                       - 'fining':     Fining profile based on Ln(diameter of maximum clast) vs. sqrt isoleth area relationship
%                                       xData: square-root of isopleth area (km)
%                                       yData: clast diameter (cm)
%
%       'yscale':       Defines the scale of the yaxis   
%                       - 'ln':         Natural log (default)
%                       - 'log10':      Log base 10
%                       - 'linear':     Linear
%
%       'maxDistancee': Defines the maximum extent of extrapolation in distal part for plotting. 1 means 100%, i.e. the distance to the most distal point is doubled
%                       Default = 1
%
%       'fits2plot':    Defines which fits to plot. If passed, it should be passed as a 1xlength(fitTypes) boolean vector
%                       Example: if fitType = {'exponential', 'powerlaw} and 'maxDistance', [1,0] is specified, only the exponential fit will be plotted
%                       All fits are plotted by default
%
%       'plotType':     - 'subplot':    Plots are inside a subplot (default)
%                       - 'separate':   Individual plots
%
%       'runMode':      - 'single':     Single volume fit (default)
%                       - 'probabilistic': Monte Carlo simulations following the method of Biass et al. (2014). If  activated, the following arguments are required
%
%       'nbRuns':       Number of runs of the Monte Carlo simulation
%
%       'xError':       Error in % for each measurement in xData. Should be passed as a vector of the same size as xData.
%                       Default is 10% on all points
%
%       'yError':       Error in % for each measurement in yData. Should be passed as a vector of the same size as yData.
%                       Default is 10% on all points
%
%       'CError':       Error in % on the distal integration limit of the power-law. 
%
%       'errorType':    Error distribution around each point.
%                       - 'normal':     Gaussian distribution where the error is 3 sigma of the distribution (default)
%                       - 'uniform':    Uniform distribution
%
%       'errorBound':   Percentiles used for reporting the error. Default is [5,95] (i.e. 5th and 95th percentiles)
%
 

addpath('Dependencies/')
%% Define the main storage structures
C = struct;     % Configuration structure
V = struct;     % Volume structure

%% Check optional arguments
% Set default values and defines C as the configuration structure
C.deposit       = 'volume';
C.yscale        = 'ln';                             % Scale of the y axis
C.nbRuns        = 0;                                % Number of runs of the Monte-Carlo simulation
C.errorType     = 'normal';                         % Error envelop
C.errorBound    = [5,95];                           % Percentiles used for volume estimate
C.xError        = ones(size(xData)).*10;            % Error on xData
C.yError        = ones(size(yData)).*10;            % Error on yData
C.CError        = 10;                               % Error on the distal integration limit
C.runMode       = 'single';                         % Single/probabilistic mode        
C.maxDist       = 1;                                % Defines the maximum extent of extrapolation in distal part. 1 means 100%, i.e. the distance to the most distal point is doubled
C.fit2plot      = true(size(fitType));              % Defines which distributions to plot    
C.plotType      = 'subplot';                        % Defines if plots are part of a subplot or separate figures

% Check deposit
if ~isempty(findCell(varargin, 'deposit'))                                  
    if isempty(findCell({'volume', 'mass', 'thinning', 'fining'}, lower(varargin{findCell(varargin, 'deposit')+1})))
        error('deposit accepts ''volume'', ''mass'', ''thinning''or ''fining''');
    end
    C.deposit = varargin{findCell(varargin, 'deposit')+1};
end
% Check yscale
if ~isempty(findCell(varargin, 'yscale'))                               
    if isempty(findCell({'ln', 'log10', 'linear'}, lower(varargin{findCell(varargin, 'yscale')+1})))
        error('yscale accepts ''ln'', ''log10'' or ''linear''');
    end
    C.yscale = varargin{findCell(varargin, 'yscale')+1};
end
% Check runMode
if ~isempty(findCell(varargin, 'runMode'))                               
    if isempty(findCell({'single', 'probabilistic'}, lower(varargin{findCell(varargin, 'runMode')+1})))
        error('runMode accepts ''single'' or ''probabilistic''');
    end
    C.runMode = varargin{findCell(varargin, 'runMode')+1};
end
% Check maxDist
if ~isempty(findCell(varargin, 'maxDistance'))                               
    C.maxDist = varargin{findCell(varargin, 'maxDistance')+1};  
end
% Check Fits to plot
if ~isempty(findCell(varargin, 'fit2plot'))         
    C.fit2plot    = logical(varargin{findCell(varargin, 'fit2plot')+1});  
end
% Check plot type
if ~isempty(findCell(varargin, 'plottype'))         
    C.plotType    = varargin{findCell(varargin, 'plottype')+1};  
end

% In case the probabilistic mode is activated, extra checks
if strcmp(C.runMode, 'probabilistic')
    % Check errorType
    if ~isempty(findCell(varargin, 'errorType'))
        if isempty(findCell({'uniform', 'normal'}, lower(varargin{findCell(varargin, 'errorType')+1})))
            error('errorType accepts ''uniform'' or ''normal''');
        end
        C.runMode = varargin{findCell(varargin, 'errorType')+1};
    end
    
    % Check xError and yError
    if ~isempty(findCell(varargin, 'xError')) && ~isempty(findCell(varargin, 'yError'))
        if (size(xData,1) ~= size(varargin{findCell(varargin, 'xError')+1},1)) || (size(xData,2) ~= size(varargin{findCell(varargin, 'xError')+1},2)) || (size(yData,1) ~= size(varargin{findCell(varargin, 'yError')+1},1)) ||(size(yData,2) ~= size(varargin{findCell(varargin, 'yError')+1},2))
            error('The size and dimensions of xError and yError should be the same as xData and yData');
        end
        C.xError = varargin{findCell(varargin, 'xError')+1};
        C.yError = varargin{findCell(varargin, 'yError')+1};
    else; error('The probabilistic mode is activated and requires to define ''xError'' and ''yError''');
    end
    
    % Check number of runs
    if ~isempty(findCell(varargin, 'nbRuns'))
        C.nbRuns = varargin{findCell(varargin, 'nbRuns')+1};
    else; error('The probabilistic mode is activated and requires to define the number of runs (i.e. ''nbRuns'')');
    end
    
    % Check the error on C
    if ~isempty(findCell(varargin, 'CError'))
        C.CError = varargin{findCell(varargin, 'CError')+1};
    else; error('The probabilistic mode is activated and requires to define the error for the distal integration limit of the power-law fit (i.e. ''CError'', in %)');
    end
    
    % Check errorBound
    if ~isempty(findCell(varargin, 'errorBound'))
        if size(varargin{findCell(varargin, 'errorBound')+1},1) ~= 1 && size(varargin{findCell(varargin, 'errorBound')+1},2) ~= 2
            error('errorBound should be specified as a symetrical 1x2 vector')
        else;    C.errorBound = varargin{findCell(varargin, 'errorBound')+1};
        end
    end
end

%% Check required arguments
if ischar(fitType); fitType = {fitType}; end    % In case only one fitType
fitType                     = sort(fitType);    % Sort alphabetically

% Check fit types
if isempty(findCell(fitType, 'exponential')) && isempty(findCell(fitType, 'powerlaw')) && isempty(findCell(fitType, 'weibull'))
    error('fitType accepts ''exponential'', ''powerlaw'' or ''weibull''. Specify multiple fits using a cell array.\n Example: {''exponential'', ''weibull''}')
end

% Check fit-dependant arguments
if ~isempty(findCell(fitType, 'exponential')) 
    % Check if break-in-slope is specified
    if isempty(findCell(varargin, 'BIS'));  V.fitProps.EXP_BIS = [];
                                            warning('No break-in-slope index specified, using one exponential segment. Use ''BIS'' argument to specify the break-in-slope indices.')
    else                                   
        V.fitProps.EXP_BIS = varargin{findCell(varargin, 'BIS')+1}; 
        % Check if more than one instance of the exponential fit is called
        if length(findCell(fitType, 'exponential')) ~= size(V.fitProps.EXP_BIS)
            error('Specify as many break-in-slopes as instances of the exponential segment');
        end
    end
end    
if ~isempty(findCell(fitType, 'powerlaw'))
    % Check if distal integration limit is specified
    if isempty(findCell(varargin, 'C'));    error('The power-law method requires to provide the distal integration limit (i.e. argument C)')
    else;                                   V.fitProps.PL_C = varargin{findCell(varargin, 'C')+1}; end
    
    % Check if distal integration limit is specified -> double check if the T0 should be specified in ln !!!!!!!
    if isempty(findCell(fitType, 'exponential')) && isempty(findCell(varargin, 'T0'))
                                            error('The power-law method requires the y intersect from the exponential value. Either add ''exponential'' to the fit type or specify a ''T0'' argument (in the natural logarithm of the y axis value)')
    elseif ~isempty(findCell(varargin, 'T0'))
                                            V.fitProps.PL_T0 = varargin{findCell(varargin, 'T0')+1};
    end
end    
if strcmpi(fitType, 'weibull')
    % Check if ranges of lambda and n are specified
    
    if ~isempty(findCell(varargin, 'lambdaRange')) && ~isempty(findCell(varargin, 'nRange'))
                                            V.fitProps.WBL_lambdaRange  = varargin{findCell(varargin, 'lambdaRange')+1}; 
                                            V.fitProps.WBL_nRange       = varargin{findCell(varargin, 'nRange')+1}; 
    % If the ranges are not directly specified but other fits are requested AND the deposit type is mass or volume
    elseif isempty(findCell(varargin, 'lambdaRange')) && length(fitType) > 1 && (strcmp(C.deposit, 'volume') || strcmp(C.deposit, 'mass'))
                                            V.fitProps.WBL_lambdaRange  = varargin{findCell(varargin, 'lambdaRange')+1}; 
                                            V.fitProps.WBL_nRange       = varargin{findCell(varargin, 'nRange')+1}; 
                                            warning('No initial ranges of lambda/n specified for the Weibull optimization. Using the average volume of all other fits to estimate initial ranges.\n In case deposit was set to ''mass'', conversion to volume using a density of 1000 kg/m3');
    else;                                   error('Initial ranges of lambda and n must be specified for the Weibull optimization with the arguments ''lambdaRange'' and ''nRange'' as a 1x2 vector containing [min, max]');
    end
end

%% Prepare data
tmp         = [reshape(xData, length(xData), 1), reshape(yData, length(yData), 1)];
[tmp, idx]  = sortrows(tmp,1);
V.xData     = tmp(:,1);
V.yData     = tmp(:,2);

% If probabilstic mode
if strcmp(C.runMode, 'probabilistic')
    tmp      = [reshape(C.xError, length(C.xError), 1), reshape(C.yError, length(C.yError), 1)];
    C.xError = tmp(idx,1);
    C.yError = tmp(idx,2);
    V.xDataP = randomize(V.xData, C.xError, C.nbRuns, C.errorType);         % Add noise to xData
    V.yDataP = randomize(V.yData, C.yError, C.nbRuns, C.errorType);         % Add noise to yData
    if ~isempty(findCell(fitType, 'powerlaw'))                              % If PL, add noise to the distal integraion limit
        V.fitProps.PL_CP = randomize(V.fitProps.PL_C, C.CError, C.nbRuns, C.errorType);
    end
end

%% Run
for iF = 1:length(fitType)
    [Vtmp,V]        = fitMe(V,C,fitType{iF});
    V.(fitType{iF}) = prepareOutput(Vtmp, V,C,fitType{iF});
end

%% Plots
% Setup plot
% Colormap
%cmap = lines(length(fitType));
cmap = [0.3639    0.5755    0.7484
        0.9153    0.2816    0.2878
        0.4416    0.7490    0.4322
        1.0000    0.5984    0.2000
        0.6769    0.4447    0.7114];

% yscale
ydata = cell(length(fitType),2);
for iF = 1:length(fitType)
    if strcmp(C.yscale, 'ln')
        ydata{iF,1} = log(V.yData);
        ydata{iF,2} = log(V.(fitType{iF}).Y);
        yl          = 'Ln ';
    elseif strcmp(C.yscale, 'log10')
        ydata{iF,1} = log10(V.yData);
        ydata{iF,2} = log10(V.(fitType{iF}).Y);
        yl    = 'Log_1_0 ';
    else
        ydata{iF,1} = V.yData;
        ydata{iF,2} = V.(fitType{iF}).Y;
        yl    = '';
    end
end

% Deposit type

% Setup figure
f = figure;
if strcmpi(C.deposit, 'volume') || strcmpi(C.deposit, 'mass')
    ax = cell(2,3);
    for i = 1:6; ax{i} = subplot(2,3,i, 'Box', 'on'); hold on; end
else
    ax = cell(2,1);
    for i = 1:2; ax{i} = subplot(2,1,i, 'Box', 'on'); hold on; end
end

% Plot data
if strcmpi(C.deposit, 'volume') || strcmpi(C.deposit, 'mass')
    plot_fit(ax{1},V,C,fitType(C.fit2plot),ydata,cmap)
    plot_inVSout(ax{2}, V,C,fitType(C.fit2plot),cmap)
    plot_volume(ax{4}, V,C,fitType(C.fit2plot),cmap)
    plot_VEI(ax{5}, V,C,fitType(C.fit2plot),cmap)
else
    plot_fit(ax{1},V,C,fitType(C.fit2plot),ydata,cmap)
    plot_inVSout(ax{2}, V,C,fitType(C.fit2plot),cmap)
end

% If activated, extract subplot to separate figures
if strcmp(C.plotType, 'separate')
    if strcmpi(C.deposit, 'volume') || strcmpi(C.deposit, 'mass')
        iEnd = 6;
    else
        iEnd = 2;
    end          
    f2 = zeros(1,iEnd);
    for i = 1:iEnd
        f2(i)   = figure;
        ah      = copyobj([ax{i}, legend(ax{i})],f2(i));
        set(ah(1), 'Position',[.1,.1,.8,.8]);
    end
    delete(f)
end

function plot_VEI(ax,V,C,fitType,cmap)
axes(ax);
% if strcmpi(C.runMode, 'probabilistic')
%     for iF = 1:length(fitType)
%         bplot(V.(fitType{iF}).VEIP, iF, 'nomean', 'nolegend', 'outliers', 'whisker', C.errorBound(1), 'color', cmap(iF,:), 'linewidth',.5, 'width', .5);
%     end
% else
%     for iF = 1:length(fitType)
%         bar(ax, iF, V.(fitType{iF}).VEI, 'FaceColor', cmap(iF,:));
%     end 
%end

for iF = 1:length(fitType)
    bar(ax, iF, V.(fitType{iF}).VEI, 'FaceColor', cmap(iF,:));
end

% Set labels
xlim([0,length(fitType)+1]);
ax.XTick       = 0:length(fitType)+2;
lab            = cell(1, length(fitType)+2);
lab(2:end-1)   = fitType;
ax.XTickLabel  = lab;
ax.YTick       = 0:7;
ax.YGrid       = 'on';
ylim([0,7]);
[~,yl]         = getLabels(C,'VEI');
ylabel(ax, yl);

function plot_volume(ax,V,C,fitType,cmap)
axes(ax);
if strcmpi(C.runMode, 'probabilistic')
    for iF = 1:length(fitType)
        bplot(V.(fitType{iF}).volumeP_km3, iF, 'nomean', 'nolegend', 'outliers', 'whisker', C.errorBound(1), 'color', cmap(iF,:), 'linewidth',.5, 'width', .5);
    end
else
    for iF = 1:length(fitType)
        bar(ax, iF, V.(fitType{iF}).volume_km3, 'FaceColor', cmap(iF,:));
    end 
end

% Set labels
xlim([0,length(fitType)+1]);
ax.XTick       = 0:length(fitType)+2;
lab            = cell(1, length(fitType)+2);
lab(2:end-1)   = fitType;
ax.XTickLabel  = lab;

[~,yl]         = getLabels(C,'volume');
ylabel(ax, yl);

function plot_inVSout(ax,V,C,fitType,cmap)
axes(ax); % Set current axes
h       = zeros(length(fitType),1);
l       = cell(length(fitType),1);
h(1)    = plot(ax, sqrt([0,max(V.yData)]), sqrt([0,max(V.yData)]), '-k', 'linewidth', 1);
l{1}    = '1:1';

for iF = 1:length(fitType)   
    if strcmpi(C.runMode, 'probabilistic')
        tmpYmP = squeeze(V.(fitType{iF}).YmP);
        for i = 1:size(tmpYmP,2)
            plot(ax, sqrt(V.yData), sqrt(tmpYmP(:,i)), '.', 'MarkerSize', 8, 'Color', cmap(iF,:));
        end
    end
h(iF+1) = plot(ax, sqrt(V.yData), sqrt(V.(fitType{iF}).Ym), 'o', 'MarkerSize', 5,'MarkerFaceColor', cmap(iF,:), 'MarkerEdgeColor','k', 'LineWidth',.2);
l{iF+1} = fitType{iF};
end

for i = 2:length(h)
    uistack(h(i), 'top');
end
legend(h,l);

[xl,yl] = getLabels(C,'InOut');
xlabel(ax, xl);
ylabel(ax, yl, 'Interpreter', 'tex');

function plot_fit(ax,V,C,fitType,ydata,cmap)
axes(ax); % Set current axes
h       = zeros(length(fitType),1);
l       = cell(length(fitType),1);
h(1)    = plot(ax, V.xData, ydata{1}, '+k');
l{1}    = 'Observations';

for iF = 1:length(fitType)   
    if strcmpi(C.runMode, 'probabilistic')
        % Vectors to store probabilistic fits
        xP = zeros(size(V.(fitType{iF}).XP,2),2);
        yP = zeros(size(V.(fitType{iF}).XP,2),2);
        for i = 1:size(V.(fitType{iF}).XP,2)
            xP(i,1) = prctile(squeeze(V.(fitType{iF}).XP(1,i,:)), C.errorBound(1));
            xP(i,2) = prctile(squeeze(V.(fitType{iF}).XP(1,i,:)), C.errorBound(2));          
            yP(i,1) = prctile(squeeze(V.(fitType{iF}).YP(1,i,:)), C.errorBound(1));
            yP(i,2) = prctile(squeeze(V.(fitType{iF}).YP(1,i,:)), C.errorBound(2));
        end
        
        % Work on yScale
        if      strcmp(C.yscale, 'ln');         yP = log(yP);
        elseif  strcmp(C.yscale, 'log10');      yP = log10(yP);
        end
        
        % Do some cleaning of the Weibull data
        if      strcmpi(fitType{iF}, 'weibull');
            idx = sum(isinf(yP),2);  % Remove  infs
            yP  = yP(not(idx),:);
            xP  = xP(not(idx),:);
        end
%         x = [xP(:,1)', fliplr(xP(:,2)')];
%         y = [yP(:,1)', fliplr(yP(:,2)')];
%         fill(x,y,cmap(iF,:),'facealpha',.25, 'edgealpha', 0);                 %#plot filled area
            plot(ax, xP(:,1), yP(:,1), ':', 'LineWidth',.5, 'Color', cmap(iF,:));
            plot(ax, xP(:,2), yP(:,2), ':', 'LineWidth',.5, 'Color', cmap(iF,:));
    end
    
    h(iF+1) = plot(ax, V.(fitType{iF}).X, ydata{iF,2}, '-', 'LineWidth',1, 'Color', cmap(iF,:));  % Plot curve  
    l{iF+1} = fitType{iF};
end

uistack(h(1), 'top');
legend(h,l);

[xl,yl] = getLabels(C,'fits');
xlabel(ax, xl);
ylabel(ax, yl);

% Get x and y labels
function [xl, yl] = getLabels(C,plotType)
% Fits plots
if strcmp(plotType, 'fits')
    if      strcmp(C.yscale, 'ln');         yl = 'Ln ';
    elseif  strcmp(C.yscale, 'log10');      yl = 'Log_1_0 ';
    else;                                   yl = '';
    end
    
    if      strcmp(C.deposit,'volume');     yl = [yl, 'Thickness (cm)'];
                                            xl = 'Isopach area^{0.5} (km)';
    elseif  strcmp(C.deposit,'thinning');   yl = [yl, 'Thickness (cm)'];
                                            xl = 'Distance from source (km)';
    elseif  strcmp(C.deposit,'mass');       yl = [yl, 'Tephra accumulation (kg/m^2)'];
                                            xl = 'Isomass area^{0.5} (km)';
    elseif  strcmp(C.deposit,'fining');     yl = [yl, 'Diameter (cm)'];
                                            xl = 'Isopleth area^{0.5} (km)';
    end 
end

% In vs out plots
if strcmp(plotType, 'InOut')
    if  strcmp(C.deposit,'mass');           yl = 'Computed mass (kg)^{0.5}';
                                            xl = 'Observed mass (kg)^{0.5}';
    elseif  strcmp(C.deposit,'fining');     yl = 'Computed diameter (cm)^{0.5}';
                                            xl = 'Observed diameter (cm)^{0.5}';
    else;                                   yl = 'Computed thickness (cm)^{0.5}';
                                            xl = 'Observed thickness (cm)^{0.5}';
    end 
end

% Volume plots
if strcmp(plotType, 'volume')
    if  strcmp(C.deposit,'mass');           yl = 'Mass (kg)'; xl = [];
    else;                                   yl = 'Volume (km^3)'; xl = [];
    end 
end
    
% VEI plots
if strcmp(plotType, 'VEI')
    if  strcmp(C.deposit,'mass');           yl = 'Magnitude'; xl = [];
    else;                                   yl = 'VEI'; xl = [];
    end 
end

function out = prepareOutput(Vtmp,V,C,fitType)
%% Prepare output
out = struct;

out.name            = fitType;
out.X               = Vtmp.X;
out.Y               = Vtmp.Y;
out.Ym              = Vtmp.Ym;
% Deterministic 
if strcmpi(C.deposit, 'volume')
    out.volume_km3  = Vtmp.volume;
    out.VEI         = floor(log10(Vtmp.volume)+5);
elseif strcmpi(C.deposit, 'mass')
   out.mass_kg      = Vtmp.volume*1e11;
   out.magnitude    = log10(out.mass_kg)-7;
end
if strcmpi(fitType, 'exponential')
    out.T0          = exp(Vtmp.F(:,1));
    out.k           = Vtmp.F(:,2);
    out.I           = Vtmp.I;
    if strcmpi(C.deposit, 'fining')
        out.bc      = .391066419./-1.*out.k;
    else
        out.bt      = .391066419./-1.*out.k;
    end
elseif strcmpi(fitType, 'powerlaw')
    out.m           = -1*Vtmp.F(2);
    out.TPl         = 10^Vtmp.F(1);
elseif strcmpi(fitType, 'weibull')
    out.theta       = Vtmp.F(1);
    out.lambda      = Vtmp.F(2);
    out.n           = Vtmp.F(3);
    % Display a warning if parameters are at the edge of bounds
    if nnz(round(V.fitProps.WBL_lambdaRange,3) == round(out.lambda,3)) > 0 || nnz(round(V.fitProps.WBL_nRange,3) == round(out.n,3)) > 0
        warning('The Weibull parameters obtained by the optimization algorithm converge to initial bounds. You might want to expand them.')
    end
end
out.r2              = Vtmp.r2;

% Probabilistic
if strcmpi(C.runMode, 'probabilistic')
    out.XP              = Vtmp.XP;
    out.YP              = Vtmp.YP;
    out.YmP             = Vtmp.YmP;
    if strcmpi(C.deposit, 'volume')
        out.volume_range= [prctile(Vtmp.volumeP, C.errorBound(1)), prctile(Vtmp.volumeP, C.errorBound(2))];
        out.volumeP_km3 = Vtmp.volumeP;
        out.VEIP        = log10(Vtmp.volumeP)+5;
    elseif strcmpi(deposit, 'mass')
        out.mass_range  = [prctile(Vtmp.volume*1e11, C.errorBound(1)), prctile(Vtmp.volume*1e11, C.errorBound(2))];
        out.massP_kg    = Vtmp.volumeP.*1e11;
        out.magnitudeP  = log10(out.massP_kg)-7;
    end
    if strcmpi(fitType, 'exponential')
        % Here fits are stored as a mxn matrix, where m is the number of
        % segments and n is the number of runs
        out.T0P         = exp(reshape(Vtmp.FP(:,1,:), size(Vtmp.FP,1),size(Vtmp.FP,3)));
        out.kP          = reshape(Vtmp.FP(:,2,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
        out.IP          = reshape(Vtmp.IP, size(Vtmp.FP,1),size(Vtmp.FP,3));
        if strcmpi(C.deposit, 'fining')
            out.bcP     = .391066419./-1.*out.k;
        else
            out.btP     = .391066419./-1.*out.k;
        end
    elseif strcmpi(fitType, 'powerlaw')
        out.mP          = -1.*reshape(Vtmp.FP(:,2,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
        out.TPlP        = 10.^reshape(Vtmp.FP(:,1,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
    elseif strcmpi(fitType, 'weibull')
        out.thetaP      = reshape(Vtmp.FP(:,1,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
        out.lambdaP     = reshape(Vtmp.FP(:,2,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
        out.nP          = reshape(Vtmp.FP(:,3,:), size(Vtmp.FP,1),size(Vtmp.FP,3));
    end
    out.r2P             = reshape(Vtmp.r2P, size(Vtmp.FP,1),size(Vtmp.FP,3));
end

function [R,V] = fitMe(V,C,fitType)
%% Fits & volume

R = struct;

% Deterministic approach
if strcmpi(fitType, 'exponential')
    [R.F,R.X,R.Y,R.r2,R.I,R.Ym] = fitEXP(V.xData, V.yData, V.fitProps.EXP_BIS, C);
    R.(C.deposit)               = volEXP(exp(R.F(:,1)), R.F(:,2));
    
elseif strcmpi(fitType, 'powerlaw')
    % Check T0
    if isfield(V.fitProps, 'PL_T0');    T0 = V.fitProps.PL_T0;
    else;                               T0 = V.exponential.T0(1);          % Here I take the T0 from the proximal segment
    end
    [R.F,R.X,R.Y,R.r2,R.Ym]     = fitPL(V.xData, V.yData, T0, C);    
    R.(C.deposit)               = volPL(T0, -1*R.F(2), 10^R.F(1), V.fitProps.PL_C);
    
elseif strcmpi(fitType, 'weibull')
    % If optimization ranges are not defined, use volume
    if ~isfield(V.fitProps, 'WBL_lambdaRange')
        tmp = [];
        if isfield(V.('exponential'), 'volume_km3')
            tmp(length(tmp)+1) = V.exponential.volume_km3;
        elseif isfield(V.('exponential'), 'mass_kg')
            tmp(length(tmp)+1) = V.exponential.mass_kg/1e3;                    % Conversion from mass to volume using density of 1000 kg/m3
        end
        if isfield(V.('powerlaw'), 'volume_km3')
            tmp(length(tmp)+1) = V.powerlaw.volume_km3;
        elseif isfield(V.('powerlaw'), 'mass_kg')
            tmp(length(tmp)+1) = V.powerlaw.mass_kg/1e3;                    % Conversion from mass to volume using density of 1000 kg/m3
        end
        [V.fitProps.WBL_lambdaRange, V.fitProps.WBL_nRange] = get_WBL_ranges(mean(tmp));
    end       
    [R.F,R.X,R.Y,R.r2,R.Ym]     = fitWBL(V.xData, V.yData, V.fitProps.WBL_lambdaRange, V.fitProps.WBL_nRange, C);
    R.(C.deposit)               = volWBL(R.F(1), R.F(2), R.F(3));
end

% Probabilistic approach
fprintf(1,'\n');
disp(['Fitting ', fitType]);
waittext(0,'init');
if strcmp(C.runMode, 'probabilistic')
    for iR = 1:C.nbRuns
        waittext(iR,'percent');
        if strcmpi(fitType, 'exponential')
            [R.FP(:,:,iR), R.XP(:,:,iR), R.YP(:,:,iR), R.r2P(:,:,iR), R.IP(:,:,iR), R.YmP(:,:,iR)] =...
                fitEXP(V.xDataP(:,:,iR), V.yDataP(:,:,iR), V.fitProps.EXP_BIS, C);
            R.([C.deposit,'P'])(iR) = volEXP(exp(R.FP(:,1,iR)), R.FP(:,2,iR));
            
        elseif strcmpi(fitType, 'powerlaw')
            if isfield(V.fitProps, 'PL_T0');    T0 = ones(1,1,C.nbRuns).*V.fitProps.PL_T0;
            else;                               T0 = V.exponential.T0P(1,:);
            end     
            [R.FP(:,:,iR), R.XP(:,:,iR), R.YP(:,:,iR), R.r2P(:,:,iR), R.YmP(:,:,iR)]  =...
                fitPL(V.xDataP(:,:,iR), V.yDataP(:,:,iR), T0(iR), C);
            R.([C.deposit,'P'])(iR) = volPL(T0(iR), -1*R.FP(:,2,iR), 10^R.FP(:,1,iR), V.fitProps.PL_CP(iR));
            
        elseif strcmpi(fitType, 'weibull')
            [R.FP(:,:,iR), R.XP(:,:,iR), R.YP(:,:,iR), R.r2P(:,:,iR), R.YmP(:,:,iR)]  =...
                fitWBL(V.xDataP(:,:,iR), V.yDataP(:,:,iR), V.fitProps.WBL_lambdaRange, V.fitProps.WBL_nRange, C);
            R.([C.deposit,'P'])(iR) = volWBL(R.FP(:,1,iR), R.FP(:,2,iR), R.FP(:,3,iR));
        end
    end
end

%% MISC FUNCTIONS
% Finds the cell index of a string in a cell array
function idx = findCell(cell2find, str2find)
cell2find   = cellfun(@num2str, cell2find, 'UniformOutput', false);
%[~,idx]     = find(not(cellfun('isempty', strfind(cell2find, str2find)))==1);
idx     = find(not(cellfun('isempty', strfind(cell2find, str2find))),1);

%% STOCHASTIC FUNCTIONS
function dataP = randomize(data, xError, nbRuns, errorType)
dataP = repmat(data, 1, 1, nbRuns);
dataE = repmat(xError, 1, 1, nbRuns);

if strcmpi(errorType, 'normal')
    dataP = dataP + dataP .* ((dataE./3).*randn(size(dataE))./100);
elseif strcmpi(errorType, 'uniform')
    dataP = dataP + dataP .* ((-dataE+(dataE-(-dataE)) .* rand(size(dataE)))./100);
end

% Returns the percentile p of the vector X
function yi = prctile(X,p)
x=X(:);
if length(x)~=length(X)
    error('please pass a vector only');
end
n = length(x);
x = sort(x);
Y = 100*(.5 :1:n-.5)/n;
x = [min(x); x; max(x)];
Y = [0 Y 100];
yi = interp1(Y,x,p);

%% EXPONENTIAL FUNCTIONS
function V = volEXP(T0, k)
% Calculate the volume with the method of Fierstein and Nathenson (1992)
% T0    Extrapolated thickness at A = 0 (cm)
% k     Slope of the segment

T0 = T0/10^5;

%% 1 segment
if size(T0,1) == 1
    % Equation (12) of Fierstein and Nathenson (1992) 
    V = 2*(T0)/k^2;       
end    

%% 2 segments
if size(T0,1) == 2
    A = (log(T0(2))-log(T0(1)))/(k(1)-k(2));
    k = -k;
    % Equation (18) of Fierstein and Nathenson (1992) 
    V = 2*T0(1)/k(1)^2 + 2*T0(1)*((k(2)*A+1)/k(2)^2 - (k(1)*A+1)/k(1)^2)*exp(-k(1)*A);                % Equation 18 of Fierstein and Nathenson (1992)
end

%% 3 segments   
if size(T0,1) == 3
    A1 = (log(T0(2))-log(T0(1)))/(k(1)-k(2));
    A2 = (log(T0(3))-log(T0(2)))/(k(2)-k(3));
    k = -k;
    % Equation (3) of Bonadonna and Houghton (2005)
    V = 2*T0(1)/k(1)^2 + 2*T0(1)*((k(2)*A1+1)/k(2)^2 - (k(1)*A1+1)/k(1)^2)*exp(-k(1)*A1) + 2*T0(2)*((k(3)*A2+1)/k(3)^2 - (k(2)*A2+1)/k(2)^2)*exp(-k(2)*A2);
end

function [F, X, Y, r2, I, Ym] = fitEXP(xdata, ydata, Aip, C)
% EXPONENTIAL fit for volume calulation with the method of Fierstein and Nathenson (1992)
% xdata: Square root of the area (km)
% ydata: Thickness (cm)
% Aip:   Vector containing the break-in-slopes as indices

ydata = log(ydata);

% Define segment indices
if isempty(Aip) %length(Aip) == 1 && Aip == 0         % One segment
    idx     = [1, length(xdata)];
else                                    % Multiple segments
    idx     = zeros(length(Aip)+2,1);
    idx(1)  = 1;
    idx(end)= length(xdata);
    for i = 1:length(Aip)
        idx(i+1) = Aip(i);
    end
end

% Fit
F	= zeros(length(idx)-1,2); % Fits
I	= zeros(length(idx)-1,1); % Intersection

for i = 1:length(idx)-1
    % Define range of data to perform the fit
    if i == 1  
        idxS = idx(i);          % Here, the intersection of the two segments is not fixed at a given point
    else
        idxS = idx(i)+1;        % In case there are multiple segments, this ensures that the same point is not used twice to define the fit.
    end
    idxE    = idx(i+1);
    
    % Fit data
    N       = size(xdata(idxS:idxE),1);
    F(i,2)  = (sum(xdata(idxS:idxE).*ydata(idxS:idxE)) - sum(xdata(idxS:idxE))*sum(ydata(idxS:idxE))/N ) / ( sum(xdata(idxS:idxE).^2) - sum(xdata(idxS:idxE))^2/N ); % Slope
    F(i,1)  = mean(ydata(idxS:idxE))-F(i,2)*mean(xdata(idxS:idxE));
end

% Calculate segments intersections
X 	= cell(length(idx)-1,1);  % X data for extrapolation
Y   = cell(length(idx)-1,1);  % Y(X)
Ym  = cell(length(idx)-1,1);  % Y(X) for observed points
r2  = zeros(length(idx)-1,1); % R-square
Xt  = 0:1:1000;

% 1 segment
if isempty(Aip)
    X{1} = [0, xdata(end) + xdata(end)*C.maxDist];
    Y{1} = exp(F(1,1)).*exp(F(1,2).*X{1});
    Ym{1}= exp(F(1,1)).*exp(F(1,2).*xdata);
    r2(1)= rsquare(exp(ydata), exp(F(i,1)).*exp(F(i,2).*xdata)); 
% Multiple segments
else
    % Find intersection of 2 segments
    for i = 1:length(idx)-1
        % First segment
        if i == 1                       
            I(i)    = polyxpoly(Xt, exp(F(i,1)).*exp(F(i,2).*Xt), Xt, exp(F(i+1,1)).*exp(F(i+1,2).*Xt));
            X{i}    = [0,I(i)];
            idxS = idx(i);
        end       
        % Last segment
        if i == length(idx)-1           
            X{i}    = [X{i-1}(end), xdata(end) + xdata(end)*C.maxDist];
            idxS = idx(i)+1;
        end        
        % Middle segments
        if i > 1 && i < length(idx)-1   
            I(i)    = polyxpoly(Xt, exp(F(i,1)).*exp(F(i,2).*Xt), Xt, exp(F(i+1,1)).*exp(F(i+1,2).*Xt));
            X{i}    = [X{i-1}(end), I(i)];
            idxS = idx(i)+1;
        end
        idxE = idx(i+1);
        % Calculate Y(X)
        Y{i}    = exp(F(i,1)).*exp(F(i,2).*X{i});
        Ym{i}   = exp(F(i,1)).*exp(F(i,2).*xdata(idxS:idxE));
        % R-Square
        r2(i)   = rsquare(exp(ydata(idxS:idxE)), exp(F(i,1)).*exp(F(i,2).*xdata(idxS:idxE))); 
    end
end
% Prepare output data
X   = reshape(cell2mat(X'), size(cell2mat(X),1)*size(cell2mat(X),2), 1)';
Y   = reshape(cell2mat(Y'), size(cell2mat(Y),1)*size(cell2mat(Y),2), 1)';
Ym  = cell2mat(Ym)';

%% POWER-LAW FUNCTIONS
function V = volPL(T0, m, Tpl, C)
% Calculate the volume with the method of Bonadonna and Houghton (2005)
% Theta, lambda, n: Fits obtained with the bc2012 function
% T0:    Maximum thickness
% m:     Power law exponent
% TPl:   Coefficient (= Tpl)
% C:     Distal integration limit (km)

T0    = T0/10^5;
Tpl   = Tpl/10^5;

% Equation 6 of Bonadonna and Houghton (2005)
B = ((T0/Tpl)^(-1/m));
V = (2*Tpl/(2-m)) * (C^(2-m) - B^(2-m));       

function [F, X, Y, r2, Ym] = fitPL(xdata, ydata, T0, C)
% POWER-LAW fits for volume calulation with the method of Bonadonna and Houghton (2005)
% Note that we use here the commonly used approximation of an exponential
% fit of log10(x), log10(y) to approximate a power-law
% xdata: Square root of the area (km)
% ydata: Thickness (cm)
% T0:    Intercept obtained with the method of Fierstein and Nathenson (1992)

N       = size(xdata,1);
xdata   = log10(xdata);
ydata   = log10(ydata);

F(2)    = ( sum(xdata.*ydata) - sum(xdata)*sum(ydata)/N ) / ( sum(xdata.^2) - sum(xdata)^2/N ); % Slope
F(1)    = mean(ydata)-F(2)*mean(xdata);

B       = (T0/10^F(1))^(-1/(-1*F(2)));                                      % Calculate proximal integration limit
X       = linspace(B, 10^xdata(end) + 10^xdata(end)*C.maxDist, 100);        % X vector
Y       = 10^F(1).*X.^(F(2));                                               % Y(X)
Ym      = 10^F(1).*10.^xdata.^(F(2));
r2      = rsquare(10.^ydata, 10^F(1).*(10.^xdata).^(F(2)));                 % R-square

%% WEIBULL FUNCTIONS
function V = volWBL(theta, lambda, n)
% Calculate the volume with the method of Bonadonna and Costa (2012)
% Theta, lambda, n: Fits obtained with the bc2012 function
% Theta:    Thickness scale (cm)
% Lambda:   Decay length scale of deposit thinning (km)
% n:        Shape parameter
V = (2*theta/10^5*lambda^2)/n;       % Equation 3 of Bonadonna and Costa (2012)

function [F, X, Y, r2, Ym] = fitWBL(xdata, ydata, lam_r, n_r, C)
wbl     = @(x02)weibe(x02, xdata, ydata); % Weibull function
F       = fminsearchbnd(wbl, [.5, .5, .5], [.1, lam_r(1) ,n_r(1)], [5000, lam_r(2), n_r(2)]); % Optimization algorithm

% Theta calculated from eq 2 of Daggit et al.
F(1)    = F(2)^(F(3)-2) * sum( (xdata.^(F(3)-2) ./ ydata).*exp(-(xdata/F(2)).^F(3)) ) .* ...
        (sum( ((xdata.^(F(3)-2) ./ ydata).*exp(-(xdata./F(2)).^F(3)) ).^2)).^-1; 

X       = linspace(0, xdata(end) + xdata(end)*C.maxDist, 100);              % X vector
Y       = F(1).*(X/F(2)).^(F(3)-2) .* exp(-((X/F(2)).^F(3)));               % Y(X)
Ym      = F(1).*(xdata/F(2)).^(F(3)-2) .* exp(-((xdata/F(2)).^F(3)));
r2      = rsquare(ydata, F(1).*(xdata/F(2)).^(F(3)-2) .* exp(-((xdata/F(2)).^F(3))));  % R-square

function F = weibe(x, xdata, ydata)
% WEIBULL fits for volume calulation with the method of Bonadonna and Costa (2012)
% xdata: Square root of the area (km)
% ydata: Thickness (cm)
% lam_r: Range of lambda values used in optimisation algorithm
% n_r:   Range of n values used in optimisation algorithm

x(1) = x(2)^(x(3)-2) * sum( (xdata.^(x(3)-2) ./ ydata).*exp(-(xdata/x(2)).^x(3)) ) .* ...
    (sum( ((xdata.^(x(3)-2) ./ ydata).*exp(-(xdata./x(2)).^x(3)) ).^2)).^-1; % Theta calculated from eq 2 of Daggit et al.
y  = x(1)*(xdata/x(2)).^(x(3)-2).*exp(-1*(xdata/x(2)).^x(3));
F  = sum(((ydata - y)./ydata).^2) + log(sum(((ydata - y)./ydata).^2));

function [lam_r, n_r] = get_WBL_ranges(vol)
% In case the user did not specify it, this function returns typical ranges 
% of lambda and n values considering a VEI based on the mean of the volume 
% values obtained with the methods of Fierstein and Nathenson (1992) 
% and Bonadonna and Houghton (2005) to estimate the VEI.
% vol:  Volume (km3)

if vol >= 10
    lam_r   = [50, 1000];
    n_r     = [10, 100];
elseif vol < 10 && vol >= 1
    lam_r   = [20, 400];
    n_r     = [10, 100];
elseif vol < 1 && vol >= 0.1
    lam_r   = [5, 100];
    n_r     = [5, 50];
elseif vol < 0.1 && vol >= 0.01
    lam_r   = [1, 20];
    n_r     = [1, 10];
elseif vol < 0.01 && vol >= 0.001
    lam_r   = [0.5, 5];
    n_r     = [1, 10];
elseif vol < 0.001
    lam_r   = [0.2, 20];
    n_r     = [1, 10];
end
    
    
