clear

%% Isopach
isopach_thick   = [90, 40, 20, 10, 8, 5, 4, 3];
isopach_area    = [0.4829, 1.4988, 2.6885, 6.2340, 9.7683, 16.2802, 26.1925, 36.3183];
%tephraFits(sqrt(isopach_area),isopach_thick,{'exponential', 'powerlaw', 'weibull'}, 'C', 20, 'CError', 20, 'deposit', 'volume', 'BIS', 3, 'runMode', 'probabilistic', 'nbRuns', 10, 'xError', ones(size(isopach_thick)).*10, 'yError', ones(size(isopach_thick)).*10)
tephraFits(sqrt(isopach_area),isopach_thick,{'exponential', 'powerlaw', 'weibull'}, 'C', 20, 'CError', 20, 'deposit', 'isopach', 'BIS', 3)

% Thinning
lobe_thickness  = [12.5000   96.0000   40.0000   16.0000   14.0000   16.0000    6.0000   13.0000    7.0000   12.0000   10.0000   26.0000];
lobe_distance   = [0.5278    0.5905    0.9964    1.4822    1.6171    1.6600    1.8847    2.0655    2.0872    2.4066    2.5733    1.2893];
%tephraVolume3(sqrt(lobe_distance),lobe_thickness,{'exponential', 'powerlaw'},'deposit', 'thinning', 'C', 20)

% Isopleth
isopleth_diam   = [9,7,6,5];
isopleth_area   = [1.5500, 5.7500, 10.0000, 26.5000];

