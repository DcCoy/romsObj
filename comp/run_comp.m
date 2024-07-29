% Script to extract ROMS fields from two similar simulations and compare them 
% (see romsComp.m for a list of available diagnostics)
% If changes are desired (i.e., different depth limits, colormaps, etc) apply them via comp_overrides.m
%
% The routine will only process either:
%   (1) monthly-averaged output, in the form of 12 output files
%   (2) annually-averaged output, in the form of 1 output file

% Grab paths
run('../romsOpt.m');

% Run options
simName    = 'peru_chile_0p1';
runNames   = {'microbes_eth_facultative_tune0','microbes_eth_obligate_tune0'};
file       = {109:120,49:60}; %Y059 M1-12 vs Y054 M1-12
plotchoice = [1 2 3 4 5 6 7];

% Initialize objects
for i = 1:length(runNames)
	obj(i) = initROMS(romsObj,simName,runNames{i});
end

% Call romsComp
obj = romsComp(obj,file,[plotchoice]);
