function V = MER(height, model, varargin)

% Example
% MER(12, {'WW87','M09','DB12'}, 'runMode', 'probabilistic', 'heightError', 20, 'nbRuns', 10, 'cWW87Error', 20, 'cM09Error', 20, 'windSpeed', 10, 'windSpeedError', 20, 'ventElevation', 1, 'tropHeight', 12, 'tropHeightError', 20)

if ischar(model)
    model = {model}; 
end    

%% Check input arguments
C.nbRuns        = 0;                                % Number of runs of the Monte-Carlo simulation
C.errorType     = 'normal';                         % Error envelop
C.errorBound    = [5,95];                           % Percentiles used for volume estimate
C.runMode       = 'single';                         % Single/probabilistic mode  

% Check arguments for Wilson and Walker (1987)
if ~isempty(findCell(model, 'WW87'))
    if ~isempty(findCell(varargin, 'cWW87'))
        % Check constant
        C.cWW87 = varargin{findCell(varargin, 'cWW87')+1};
    else
        C.cWW87 = .236;
        warning('No constant for Wilson and Walker (1987) was defined. A value of 0.236 was assumed.');
    end
end
% Check arguments for Mastin et al. (2009)
if ~isempty(findCell(model, 'M09'))
    % Check constant
    if ~isempty(findCell(varargin, 'cM09'))
        C.cM09 = varargin{findCell(varargin, 'cM09')+1};
    else
        C.cM09 = 2;
        warning('No constant for Mastin et al. (2009) was defined. A value of 2 was assumed.');
    end
end
% Check arguments for Degruyter and Bonadonna (2012)
if ~isempty(findCell(model, 'DB12'))
    % Check wind speed
    if ~isempty(findCell(varargin, 'windSpeed'))
        C.windSpeed = varargin{findCell(varargin, 'windSpeed')+1};
    else
        error('The model of Degruyter and Bonadonna requires the wind speed at the tropopause');
    end
    % Check vent elevation
    if ~isempty(findCell(varargin, 'ventElevation'))
        C.ventElevation = varargin{findCell(varargin, 'ventElevation')+1};
    else
        error('The model of Degruyter and Bonadonna requires the vent elevation');
    end
    % Check tropopause height
    if ~isempty(findCell(varargin, 'tropHeight'))
        C.tropHeight = varargin{findCell(varargin, 'tropHeight')+1};
    else
        C.tropHeight = 12;
        warning('The model of Degruyter and Bonadonna requires the tropopause height. A value of 12 km was assumed');
    end
end
% Check runMode
if ~isempty(findCell(varargin, 'runMode'))                               
    if isempty(findCell({'single', 'probabilistic'}, lower(varargin{findCell(varargin, 'runMode')+1})))
        error('runMode accepts ''single'' or ''probabilistic''');
    end
    C.runMode = varargin{findCell(varargin, 'runMode')+1};
else
    C.runMode = 'single';                         % Single/probabilistic mode        
end



% In case the probabilistic mode is activated, extra checks
if strcmp(C.runMode, 'probabilistic')
    % Check errorType
    if ~isempty(findCell(varargin, 'errorType'))
        if isempty(findCell({'uniform', 'normal'}, lower(varargin{findCell(varargin, 'errorType')+1})))
            error('errorType accepts ''uniform'' or ''normal''');
        end
        C.errorType = varargin{findCell(varargin, 'errorType')+1};
    else
        C.errorType     = 'normal';                         % Error envelop
        warning('No error envelop was specified, using a normal distribution');
    end
    
    % Check number of runs
    if ~isempty(findCell(varargin, 'nbRuns'))
        C.nbRuns = varargin{findCell(varargin, 'nbRuns')+1};
    else; error('The probabilistic mode is activated and requires to define the number of runs (i.e. ''nbRuns'')');
    end
    
    % Check errorBound
    if ~isempty(findCell(varargin, 'errorBound'))
        C.errorBound = varargin{findCell(varargin, 'errorBound')+1};
    else
        C.errorBound = [5,95];
        warning('No error bounds were specified, using the 5th and 95th percentiles');     
    end
    
    % Check error on plume height
    if ~isempty(findCell(varargin, 'heightError'))
        C.heightError = varargin{findCell(varargin, 'heightError')+1};
    else
        error('The probabilistic mode is activated and requires to define the error on the plume height (i.e. ''heightError'')');     
    end
    
    % Check arguments for Wilson and Walker (1987)
    if ~isempty(findCell(model, 'WW87'))  
       if ~isempty(findCell(varargin, 'cWW87Error'))
           % Check constant
           C.cWW87Error = varargin{findCell(varargin, 'cWW87Error')+1};
       else
           error('The probabilistic mode is activated and requires to define the error on the Wilson and Walker (1987) constant (i.e. ''cWW87Error'')');
       end
    end
    
    % Check arguments for Mastin et al. (2009)
    if ~isempty(findCell(model, 'M09'))  
       if ~isempty(findCell(varargin, 'cM09Error'))
           % Check constant
           C.cM09Error = varargin{findCell(varargin, 'cM09Error')+1};
       else
           error('The probabilistic mode is activated and requires to define the error on the Mastin et al. (2009) constant (i.e. ''cM09Error'')');
       end
    end
    
    % Check arguments for Degruyter and Bonadonna (2012)
    if ~isempty(findCell(model, 'DB12'))
        % Check wind speed
        if ~isempty(findCell(varargin, 'windSpeedError'))
            C.windSpeedError = varargin{findCell(varargin, 'windSpeedError')+1};
        else
            error('The probabilistic mode is activated and requires to define the error on the wind speed at the tropopause for the model of Degruyter and Bonadonna (2012)');
        end

        % Check tropopause height
        if ~isempty(findCell(varargin, 'tropHeightError'))
            C.tropHeightError = varargin{findCell(varargin, 'tropHeightError')+1};
        else
            C.tropHeight = 12;
            error('The probabilistic mode is activated and requires to define the error on the tropopause height for the model of Degruyter and Bonadonna (2012)');
        end
    end
end


%% Run calculations
% Single run
for iM = 1:length(model)
    V.height = height;
    if strcmp(model{iM}, 'WW87')
        V.WW87.MER          = (height./C.cWW87).^4;
        V.WW87.constant     = C.cWW87;
    elseif strcmp(model{iM}, 'M09')
        V.M09.MER           = ((height./C.cM09).^(1/.241)).*2500;
        V.M09.constant      = C.cM09;
    elseif strcmp(model{iM}, 'DB12')
        V.DB12.MER          = get_MER_DB12(height-C.ventElevation, C.windSpeed, C.tropHeight);
        V.DB12.ventElevation= C.ventElevation;
        V.DB12.windSpeed    = C.windSpeed;
        V.DB12.tropHeight   = C.tropHeight;
    end
end

% Probabilistic run
for iM = 1:length(model)
    V.heightP = randomize(V.height,C.heightError, C.nbRuns, C.errorType);
    if strcmp(model{iM}, 'WW87')
        V.WW87.constantP    = randomize(V.WW87.constant,C.cWW87Error, C.nbRuns, C.errorType);
        V.WW87.MERP         = (V.heightP./V.WW87.constantP).^4;
    elseif strcmp(model{iM}, 'M09')
        V.M09.constantP     = randomize(V.M09.constant,C.cM09Error, C.nbRuns, C.errorType);
        V.M09.MERP          = ((V.heightP./C.cM09Error).^(1/.241)).*2500;
    elseif strcmp(model{iM}, 'DB12')
        V.DB12.windSpeedP   = randomize(V.DB12.windSpeed,C.windSpeedError, C.nbRuns, C.errorType);
        V.DB12.tropHeightP  = randomize(V.DB12.tropHeight,C.tropHeightError, C.nbRuns, C.errorType);
        V.DB12.MERP         = zeros(C.nbRuns,1);
        for i = 1:C.nbRuns
            V.DB12.MERP(i)  = get_MER_DB12(V.heightP(i) - V.DB12.ventElevation, V.DB12.windSpeedP(i), V.DB12.tropHeightP(i));
        end
    end
end

%% STOCHASTIC FUNCTIONS
function dataP = randomize(data, dataE, nbRuns, errorType)
%dataP = repmat(data, 1, 1, nbRuns);
%dataE = repmat(xError, 1, 1, nbRuns);

if strcmpi(errorType, 'normal')
    dataP = data + data .* ((dataE./3).*randn(nbRuns,1)./100);
elseif strcmpi(errorType, 'uniform')
    dataP = data + data .* ((-dataE+(dataE-(-dataE)) .* rand(nbRuns,1))./100);
end