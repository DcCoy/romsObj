% Set your romsObj options here
addpath scripts

% PATHS
rootPath      = '/data/project1/demccoy/romsObj/';                  % where romsObj, romsDiag, romsComp, romsOpt live
simPath       = '/data/project2/model_output/';                     % where ROMS output is stored (or linked)
scriptsPath   = '/data/project2/demccoy/romsObj/scripts/';          % where additional matlab scripts live 
valiPath      = '/data/project1/demccoy/ROMS/validation/products/'; % where validation data lives
tmpfigsPath   = '/data/project1/demccoy/tmpfigs/';                  % where temporary figures will be stored 

% DEFAULT FIGURE OPTIONS
figsFormat    = 'png';            % format of output figures (jpg, png, pdf)
figVisible    = 'off';            % use 'on' to make figure windows visible 
background    = rgb('LightGray'); % color to be used in backgrounds of plots 
coastcolor    = rgb('DimGray');   % color to be used for NaN points (coasts, bathymetery, etc) 
fontsize      = 8;                % fontsize in figures
figtype       = 'mfig';           % 'sfig' = 9cm wide, 'mfig' = 14cm, 'lfig' = 19cm   
coast         = 'coast';          % choose from: coast, crude, low, high, intermediate, or full (coast map resolution)  
ticks         = 1;                % to show lon/lat ticks on x/y axes
