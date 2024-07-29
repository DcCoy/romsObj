% Script to extract ROMS fields and compare them against various validation productions
% (see romsDiag.m for a list of available diagnostics)
% If changes are desired (i.e., different depth limits, colormaps, etc) apply them via diag_overrides.m
%
% The routine will only process either:
%   (1) monthly-averaged output, in the form of 12 output files
%   (2) annually-averaged output, in the form of 1 output file

% Grab paths
run('../romsOpt.m');

% Run options
if (0)
    % Pacmed
    simName    = 'pacmed_0p25';
    runName    = 'pdamien_loop3';
    file       = 1; 
else
    % Peru
    simName    = 'peru_chile_0p1';
    runName    = 'microbes_eth_obligate_tune0';
    file       = 109:120;
end

% Plot choices (see romsDiag)
plotchoice = [1:7];

% Initialize object
obj = initROMS(romsObj,simName,runName);

% Call romsDiag
obj = romsDiag(obj,file,[plotchoice]);
