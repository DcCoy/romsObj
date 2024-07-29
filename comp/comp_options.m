%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% romsComp OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If possible, don't edit these settings.
% Instead, apply simulation-specific overrides
% via 'comp_overrides' under simName folders
%
%
% General description of options
% .opt   = switch to plot fields (1 == yes)
% .vars  = which ROMS vars to plot (same length as .opt, .lims, .dlims, .cmaps)
% .zdeps = choice of depths
% .lims  = colorbar limits
% .dlims = colorbar limits for difference plot
% .xlims = Lon or Lat degree limits for slicePlot(s)
% .zlims = Depth limits for slicePlot(s)

% (1) zslices
plots(1).opt   = [   1     1     1   ... % temp, salt, rho
					 1     1     1   ... % O2, NOX, PO4
					 1     1     1   ... % NO2, N2O, NH4
					 1     1     1   ... % nstar, Fe, SiO3   
					 1     1];           % DIC, Alk   
plots(1).vars  = {'temp' ,'salt','rho',  ...
				  'O2'   ,'NOX' ,'PO4',  ...
				  'NO2'  ,'N2O' ,'NH4',  ...
				  'nstar','Fe'  ,'SiO3', ...
				  'DIC'  ,'Alk'};
plots(1).zdeps = [   0            50          150         300        450    ];
plots(1).lims  = {[0      30],[0      30],[0      25],[0      18],[0      18]; ... % temp
				  [30     37],[30     37],[33.4 36.5],[33.4 35.5],[33.4 35.5]; ... % salt
				  [19     27],[19     27],[24     28],[26.5 28.6],[26.5 28.6]; ... % rho
				  [0     350],[0     350],[0     350],[0     350],[0     350]; ... % O2
				  [0      40],[0      40],[0      40],[0      50],[0      50]; ... % NOX
				  [0       3],[0       3],[0       3],[0       5],[0       5]; ... % PO4
				  [0       1],[0       1],[0       3],[0       3],[0       3]; ... % NO2
				  [0    0.02],[0    0.05],[0    0.08],[0    0.08],[0    0.08]; ... % N2O
				  [0       3],[0       3],[0       1],[0       1],[0       1]; ... % NH4
				  [-25    25],[-25    25],[-25    25],[-25    25],[-25    25]; ... % nstar
				  [0    1e-3],[0    1e-3],[0    1e-3],[0    2e-3],[0    2e-3]; ... % Fe
				  [0      50],[0      50],[0      50],[0     100],[0     100]; ... % SiO3
				  [1900 2250],[1900 2250],[2000 2400],[2000 2400],[2000 2400]; ... % DIC
				  [2200 2450],[2200 2450],[2300 2450],[2300 2450],[2300 2450]};    % Alk
plots(1).dlims = {[-3       3]; ... % temp
				  [-0.6   0.6]; ... % salt
				  [-0.6   0.6]; ... % rho
				  [-80     80]; ... % O2
				  [-10     10]; ... % NOX
				  [-1       1]; ... % PO4
				  [-0.6   0.6]; ... % NO2
				  [-0.05 0.05]; ... % N2O 
				  [-0.6   0.6]; ... % NH4
				  [-10     10]; ... % nstar
				  [-1e-3 1e-3]; ... % Fe
				  [-50     50]; ... % SiO3
				  [-300   300]; ... % DIC
				  [-300   300]};    % Alk
plots(1).cmaps =  {'thermal' ; ... % temp  
				   'haline'  ; ... % salt
				   'dense'   ; ... % rho
				   '-ice'    ; ... % O2
				   'tempo'   ; ... % NOX
				   'tempo'   ; ... % PO4
			  	   'tempo'   ; ... % NO2
				   'tempo'   ; ... % N2O 
				   'tempo'   ; ... % NH4
				   'balance' ; ... % nstar
				   'tempo'   ; ... % Fe
				   'tempo'   ; ... % SiO3
				   'deep'    ; ... % DIC
				   'deep'    };    % Alk

% (2) Zonal transects 
plots(2).opt   = [   1     1     1   ... % temp, salt, rho
					 1     1     1   ... % O2, NOX, PO4
					 1     1     1   ... % NO2, N2O, NH4
					 1     1     1   ... % nstar, Fe, SiO3
					 1     1];           % DIC, Alk
plots(2).vars  = {'temp' ,'salt','rho',  ...
				  'O2'   ,'NOX' ,'PO4',  ...
				  'NO2'  ,'N2O' ,'NH4',  ...
				  'nstar','Fe'  ,'SiO3', ...
				  'DIC'  ,'Alk'};
plots(2).choice = 'lat';
plots(2).lats   = 0;
plots(2).xlims  = [140 260];
plots(2).zlims  = [0 1000];
plots(2).lims   = {[0       30]; ... % temp  
				   [33.5  35.5]; ... % salt
				   [22      32]; ... % rho
				   [0      300]; ... % O2
				   [0       50]; ... % NOX
				   [0        5]; ... % PO4
				   [0        1]; ... % NO2
				   [0     0.08]; ... % N2O 
				   [0        3]; ... % NH4
				   [-25     25]; ... % nstar
				   [0     1e-3]; ... % Fe
				   [0      100]; ... % SiO3
				   [1900  2400]; ... % DIC
				   [2200  2450]};    % Alk
plots(2).dlims =  {[-5       5]; ... % temp
				   [-1.0   1.0]; ... % salt
				   [-0.6   0.6]; ... % rho
				   [-100   100]; ... % O2
				   [-10     10]; ... % NOX
				   [-1       1]; ... % PO4
				   [-0.6   0.6]; ... % NO2
				   [-0.05 0.05]; ... % N2O
				   [-0.6   0.6]; ... % NH4
				   [-10     10]; ... % nstar
				   [-1e-3 1e-3]; ... % Fe
				   [-50     50]; ... % SiO3
				   [-300   300]; ... % DIC
				   [-300   300]};    % Alk
plots(2).cmaps =  {'thermal' ; ... % temp
				   'haline'  ; ... % salt
				   'dense'   ; ... % rho
				   '-ice'    ; ... % O2
				   'tempo'   ; ... % NOX
				   'tempo'   ; ... % PO4
				   'tempo'   ; ... % NO2
				   'tempo'   ; ... % N2O
				   'tempo'   ; ... % NH4
				   'balance' ; ... % nstar
				   'tempo'   ; ... % Fe
				   'tempo'   ; ... % SiO3
				   'deep'    ; ... % DIC
				   'deep'    };    % Alk

% (3) Meridional transects
plots(3).opt   = [   1     1     1   ... % temp, salt, rho
					 1     1     1   ... % O2, NOX, PO4
					 1     1     1   ... % NO2, N2O, NH4
					 1     1     1   ... % nstar, Fe, SiO3
					 1     1];           % DIC, Alk
plots(3).vars  = {'temp' ,'salt','rho',...
				  'O2'   ,'NOX' ,'PO4',...
				  'NO2'  ,'N2O' ,'NH4',...
				  'nstar','Fe'  ,'SiO3',...
				  'DIC'  ,'Alk'};
plots(3).choice = 'lon';
plots(3).zlims  = [0 1000];
plots(3).lons   = [  210    255   272  ];
plots(3).xlims  = [-30 50;-40 20;-40 10];
plots(3).lims   = {[0       30]; ... % temp  
				   [33.5  35.5]; ... % salt
				   [22      32]; ... % rho
				   [0      300]; ... % O2
				   [0       50]; ... % NOX
				   [0        5]; ... % PO4
				   [0        1]; ... % NO2
				   [0     0.08]; ... % N2O 
				   [0        3]; ... % NH4
				   [-25     25]; ... % nstar
				   [0     1e-3]; ... % Fe
				   [0      100]; ... % SiO3
				   [1900  2400]; ... % DIC
				   [2200  2450]};    % Alk
plots(3).dlims  = {[-5       5]; ... % temp
				   [-1.0   1.0]; ... % salt
				   [-0.6   0.6]; ... % rho
				   [-100   100]; ... % O2
				   [-10     10]; ... % NOX
				   [-1       1]; ... % PO4
				   [-0.6   0.6]; ... % NO2
				   [-0.05 0.05]; ... % N2O
				   [-0.6   0.6]; ... % NH4
				   [-10     10]; ... % nstar
				   [-1e-3 1e-3]; ... % Fe
				   [-50     50]; ... % SiO3
				   [-300   300]; ... % DIC
				   [-300   300]};    % Alk
plots(3).cmaps =  {'thermal' ; ... % temp
				   'haline'  ; ... % salt
				   'dense'   ; ... % rho
				   '-ice'    ; ... % O2
				   'tempo'   ; ... % NOX
				   'tempo'   ; ... % PO4
				   'tempo'   ; ... % NO2
				   'tempo'   ; ... % N2O
				   'tempo'   ; ... % NH4
				   'balance' ; ... % nstar
				   'tempo'   ; ... % Fe
				   'tempo'   ; ... % SiO3
				   'deep'    ; ... % DIC
				   'deep'    };    % Alk

% (4)  Surface fields
plots(4).opt   = [   1         1          1            1            1            1           1             1     ];
plots(4).vars  = { 'MLD'  ,  'SSH'   ,  'NPP'   ,   'FG_N2O' ,   'sustr'  ,   'svstr'  ,    'ws'    ,    'wsc'   };
plots(4).cmaps = { 'deep' ;  'deep'  ; 'algae'  ;  'balance' ;  'balance' ;  'balance' ;   'amp'    ; 'balance'  };
plots(4).lims  = {[0  150];[0    1.5];[0   1000];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[0     1e-4];[-2e-7 2e-7]};
plots(4).dlims = {[-60 60];[-0.3 0.3];[-700 700];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[-1e-4 1e-4];[-2e-7 2e-7]};

% (5) surface chlA
plots(5).opt       = [1];
plots(5).vars      = {'SFC_CHL'};
plots(5).cmaps     = {'algae'};
plots(5).absLevs   = sqrt([0:0.01:1.0]);
plots(5).absLbls   = plots(5).absLevs.^2;
plots(5).diffLevs  = linspace(-0.4,0.4,42);

% (6) OMZ thickness
plots(6).opt       = [1];
plots(6).cmap      = {'amp'};
plots(6).omzthresh = [     0          5            10          20            50     ]; 
plots(6).lims      = {[0    500],[0     1000],[0     1000],[0     2000],[0     2000]};
plots(6).dlims     = {[-500 500],[-1000 1000],[-1000 1000],[-2000 2000],[-2000 2000]};

% (7) Transect maps
plots(7).lons = [210 255 272];
plots(7).lats = [0];

