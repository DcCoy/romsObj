function obj = romsDiag(obj,file,plotchoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to automatically generate diagnostic plots for individual ROMS simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% - obj = romsDiag(obj,file,plotchoice)
%
% Inputs:
% - obj  = roms object that is already initialized
% - file = file(s) to diagnose 
%
%    ------------------------------
%    ---- STANDARD DIAGNOSTICS ----
%    ------------------------------
%
%    1  == depth sections
%    2  == zonal transects
%    3  == meridional transects
%    4  == surface comparisons
%
%    -----------------------
%    -- OTHER DIAGNOSTICS --
%    -----------------------
%
%    5  == surface chlA map (nonlinear colorbar scale)
%    6  == OMZ thickness     (0, 5, 10, 20, 50uM)
%    7  == transect locations
%    8  == POC flux validation
%    9  == Gridded N cycle tracers obs
%    10  == deep depth sections (currently, Fe only)
%    11  == Fe vs Geotraces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle
addpath /data/project1/demccoy/ROMS/validation/n2o
addpath /data/project1/demccoy/ROMS/scripts

% Clear objects
obj = clearROMS(obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get plotchoice
tmpchoice = zeros(1,100);
tmpchoice(plotchoice) = 1;
plotchoice = tmpchoice;

% PLOT options
% (1) zslices
pltcnt = 1;
plots(pltcnt).on = plotchoice(pltcnt);
   plots(pltcnt).opt   = [   1     1     1   ... % temp, salt, rho
                             1     1     1   ... % O2, NOX, PO4
                             1     1     1   ... % NO2, N2O, NH4
                             1     1     1   ... % nstar, Fe, SiO3   
                             1     1];           % DIC, Alk   
   plots(pltcnt).vars  = {'temp' ,'salt','rho',  ...
                          'O2'   ,'NOX' ,'PO4',  ...
                          'NO2'  ,'N2O' ,'NH4',  ...
                          'nstar','Fe'  ,'SiO3', ...
                          'DIC'  ,'Alk'};
   plots(pltcnt).zdeps = [   0            50          150         300        450    ];
   plots(pltcnt).lims  = {[0      30],[0      30],[0      25],[0      18],[0      18]; ... % temp
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
   plots(pltcnt).dlims = {[-3       3]; ... % temp
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
   plots(pltcnt).cmaps =  {'thermal' ; ... % temp  
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
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
   plots(pltcnt).opt   = [   1     1     1   ... % temp, salt, rho
                             1     1     1   ... % O2, NOX, PO4
                             1     1     1   ... % NO2, N2O, NH4
                             1     1     1   ... % nstar, Fe, SiO3
                             1     1];           % DIC, Alk
   plots(pltcnt).vars  = {'temp' ,'salt','rho',  ...
                          'O2'   ,'NOX' ,'PO4',  ...
                          'NO2'  ,'N2O' ,'NH4',  ...
                          'nstar','Fe'  ,'SiO3', ...
                          'DIC'  ,'Alk'};
   plots(pltcnt).choice = 'lat';
   plots(pltcnt).lats   = 0;
   plots(pltcnt).xlims  = [140 260];
   plots(pltcnt).zlims  = [0 1000];
   plots(pltcnt).lims   = {[0       30]; ... % temp  
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
   plots(pltcnt).dlims =  {[-5       5]; ... % temp
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
   plots(pltcnt).cmaps =  {'thermal' ; ... % temp
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
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
   plots(pltcnt).opt   = [   1     1     1   ... % temp, salt, rho
                             1     1     1   ... % O2, NOX, PO4
                             1     1     1   ... % NO2, N2O, NH4
                             1     1     1   ... % nstar, Fe, SiO3
                             1     1];           % DIC, Alk
   plots(pltcnt).vars  = {'temp' ,'salt','rho',...
                          'O2'   ,'NOX' ,'PO4',...
                          'NO2'  ,'N2O' ,'NH4',...
                          'nstar','Fe'  ,'SiO3',...
                          'DIC'  ,'Alk'};
   plots(pltcnt).choice = 'lon';
   plots(pltcnt).zlims  = [0 1000];
   plots(pltcnt).lons   = [  210    255   272  ];
   plots(pltcnt).xlims  = [-30 50;-40 20;-40 10];
   plots(pltcnt).lims   = {[0       30]; ... % temp  
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
    plots(pltcnt).dlims = {[-5       5]; ... % temp
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
   plots(pltcnt).cmaps =  {'thermal' ; ... % temp
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
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
   plots(pltcnt).opt   = [   1         1          1            1            1            1           1             1     ];
   plots(pltcnt).vars  = { 'MLD'  ,  'SSH'   ,  'NPP'   ,   'FG_N2O' ,   'sustr'  ,   'svstr'  ,    'ws'    ,    'wsc'   };
   plots(pltcnt).cmaps = { 'deep' ;  'deep'  ; 'algae'  ;  'balance' ;  'balance' ;  'balance' ;   'amp'    ; 'balance'  };
   plots(pltcnt).lims  = {[0  150];[0    1.5];[0   1000];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[0     1e-4];[-2e-7 2e-7]};
   plots(pltcnt).dlims = {[-60 60];[-0.3 0.3];[-700 700];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[-1e-4 1e-4];[-2e-7 2e-7]};
    
% (5) surface chlA
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).opt       = [1];
    plots(pltcnt).vars      = {'SFC_CHL'};
    plots(pltcnt).cmaps     = {'algae'};
    plots(pltcnt).absLevs   = sqrt([0:0.01:1.0]);
    plots(pltcnt).absLbls   = plots(pltcnt).absLevs.^2;
	plots(pltcnt).diffLevs  = linspace(-0.4,0.4,42);

% (6) OMZ thickness
% NOTE: Validation products are hard coded to omzthresh
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).opt       = [1];
    plots(pltcnt).cmap      = cbrewer('seq','YlOrRd',40);
    plots(pltcnt).cmap(1,:) = [1 1 1];
    plots(pltcnt).omzthresh = [     0          5            10          20            50     ]; 
    plots(pltcnt).lims      = {[0    500],[0     1000],[0     1000],[0     2000],[0     2000]};
    plots(pltcnt).dlims     = {[-500 500],[-1000 1000],[-1000 1000],[-2000 2000],[-2000 2000]};

% (7) Transect maps
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).lons = [210 255 272];
    plots(pltcnt).lats = [0];

% (8) POC_FLUX_IN comparisons
% NOTE: Validation product depths are hard coded
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).opt   = [1];
    plots(pltcnt).cmaps = {'deep'};
    plots(pltcnt).levs  = {linspace(0,5e-4,40)};
    plots(pltcnt).dlevs = {linspace(-5e-4,5e-4,41)};

% (9) Gridded obs vs ROMS (including ratios)
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).xlims = [100 300];
    plots(pltcnt).ylims = [-40  60];

% (10) Deep Fe depth slices vs Tagliabue
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).opt   = [1];
    plots(pltcnt).vars  = {'Fe'};
    plots(pltcnt).cmaps = {'tempo'};
    plots(pltcnt).zdeps = [    600           800        1000         1500          2000    ];
    plots(pltcnt).lims  = {[0     1e-3],[0     1e-3],[0     1e-3],[0     1e-3],[0     1e-3]};
    plots(pltcnt).dlims = {[-1e-3 1e-3],[-1e-3 1e-3],[-1e-3 1e-3],[-1e-3 1e-3],[-1e-3 1e-3]};

% (11) Fe lon/lat slices vs GEOTRACES
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
    plots(pltcnt).sections = {'GP02','GP13','GP16','GP18','GP19'};
    plots(pltcnt).coord    = { 'lat', 'lat', 'lat', 'lon', 'lon'};
    plots(pltcnt).cmaps    = {'tempo'};
    plots(pltcnt).lims     = {[0     1e-3]};
    plots(pltcnt).dlims    = {[-1e-3 1e-3]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change to diagDir and load overrides options
diagDir = ['/data/project1/demccoy/ROMS/',obj.info.simName,'/analysis/diag/'];
mkdir(diagDir);
cd(diagDir);
run(['diag_overrides.m']);

% Reset pltcnt
pltcnt = 0;

% Clear workspace and begin
clearvars -except obj plots pltcnt file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PHYSICAL DIAGNOSTICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = 0;

% Get variables list
% phy only
if isempty(obj.info.bgc_avg) & isempty(obj.info.dia_avg) 
    all_vars = {obj.info.phy_avg(file(1)).Variables.Name};
% dia only
elseif isempty(obj.info.phy_avg) & isempty(obj.info.bgc_avg) 
    all_vars = {obj.info.dia_avg(file(1)).Variables.Name};
% bgc only
elseif isempty(obj.info.phy_avg) & isempty(obj.info.dia_avg) 
    all_vars = {obj.info.bgc_avg(file(1)).Variables.Name};
% bgc + dia only
elseif isempty(obj.info.phy_avg) & ~isempty(obj.info.bgc_avg) & ~isempty(obj.info.dia_avg) 
    all_vars = {obj.info.bgc_avg(file(1)).Variables.Name ...
                obj.info.dia_avg(file(1)).Variables.Name};
% bgc + phy only
elseif isempty(obj.info.dia_avg) & ~isempty(obj.info.bgc_avg) & ~isempty(obj.info.phy_avg) 
    all_vars = {obj.info.phy_avg(file(1)).Variables.Name ...
                obj.info.bgc_avg(file(1)).Variables.Name};
% all available
else
    all_vars = {obj.info.phy_avg(file(1)).Variables.Name ...
                obj.info.bgc_avg(file(1)).Variables.Name ...
                obj.info.dia_avg(file(1)).Variables.Name};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D diagnostics (zslices)
% P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
   % Variable loop
   for v = 1:length(vars)
      if ~opt(v)
         disp(['...skipping ',vars{v},'...']);
         continue
      end
      % Get data
      obj = clearROMS(obj);
      try
          obj = zslice(obj,vars(v),zdeps,file);
      catch
          disp(['...skipping ',vars{v},'...']);
          continue
      end
      obj = loadDiag(obj,vars(v),zdeps);
      % Depth loop
      for z = 1:length(zdeps)
         for d = 1:length(obj.diag.(vars{v}))
            % Get levels
            levs  = linspace(lims{v,z}(1),lims{v,z}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42); 
            % Extract data
            romsdat    = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,z,:)),3);
            diagdat    = nanmean(squeeze(obj.diag.(vars{v})(d).slice(:,:,z,:)),3);
            close all
            % Make figures
            [figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
            % ROMS figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_roms'],'-m5');
            close(figs(1));
            % Diag figure
            set(0,'CurrentFigure',figs(2));
            title([obj.diag.(vars{v})(d).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diag',num2str(d)],'-m5');
            close(figs(2));   
            % Diff figure
            set(0,'CurrentFigure',figs(3));
            title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diff',num2str(d)],'-m5');   
         end
      end
   end
   % Clear data
   obj = clearROMS(obj);
   clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zonal transects
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
   % Variable loop
   for v = 1:length(vars);
      if ~opt(v) 
         disp(['...skipping ',vars{v},'...']);
         continue
      end
      % Transect loop
      for l = 1:length(lats)
         % Get data
         obj = clearROMS(obj);
         try
             obj = sliceROMS(obj,vars(v),choice,lats(l),file,...
                'zdep',obj.grid.z_avg_dep(obj.grid.z_avg_dep<=zlims(end)));
         catch
             disp(['...skipping ',vars{v},'...']);
             continue
         end
         obj = sliceDiag(obj,vars(v),choice,lats(l),'zlim',zlims(end));
         for d = 1:length(obj.diag.(vars{v}))
             % Extract data
             romsdat = obj.data.avg.(vars{v}).slice;
             diagdat = obj.diag.(vars{v})(d).slice;
             % Get levels
             levs  = linspace(lims{v}(1),lims{v}(2),40);
             dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
             % Make figures
             [figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
                'figdim',0.5,'levels',levs,'difflevels',dlevs);
             % ROMS figure
             set(0,'CurrentFigure',figs(1));
             title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
             ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
             export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_roms'],'-m5');
             close(figs(1))
             % Diag figure
             set(0,'CurrentFigure',figs(2));
             title([obj.diag.(vars{v})(d).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
             ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
             export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diag',num2str(d)],'-m5');
             close(figs(2))
             % Difference figure
             set(0,'CurrentFigure',figs(3));
             title(['Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
             ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
             export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diff',num2str(d)],'-m5');
             close(figs(3))
         end
      end
   end
   % Clear data
   obj = clearROMS(obj);
   clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% longitude transects 
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
   % Variable loop
   for v = 1:length(vars);
      if ~opt(v) %| ~ismember(vars{v},all_vars);
         disp(['...skipping ',vars{v},'...']);
         continue
      end
      % Transect loop
      for l = 1:length(lons)
         % Get data
         obj = clearROMS(obj);
%		 try
			 obj = sliceROMS(obj,vars(v),choice,lons(l),file,...
				'zdep',obj.grid.z_avg_dep(obj.grid.z_avg_dep<=zlims(end)));
			 obj = sliceDiag(obj,vars(v),choice,lons(l),'zlim',zlims(end));
%		 catch
%             disp(['...skipping ',vars{v},'...']);
%             continue
%         end
         for d = 1:length(obj.diag.(vars{v}))
            % Extract data
            romsdat = obj.data.avg.(vars{v}).slice;
            diagdat = obj.diag.(vars{v})(d).slice;
            % Get levels
            levs  = linspace(lims{v}(1),lims{v}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
            % Make figures
            [figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
               'figdim',0.5,'levels',levs,'difflevels',dlevs);
            % ROMS figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
            ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_roms'],'-m5');
            close(figs(1))
            % Diag figure
            set(0,'CurrentFigure',figs(2));
            title([obj.diag.(vars{v})(d).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
            ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diag',num2str(d)],'-m5');
            close(figs(2))
            % Difference figure
            set(0,'CurrentFigure',figs(3));
            title(['Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
            ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diff',num2str(d)],'-m5');
            close(figs(3))
         end
      end
   end
   % Clear data
   obj = clearROMS(obj);
   clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface diagnostics
% P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
    % Variable loop
    for v = 1:length(vars);
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
        % Get data
        obj = clearROMS(obj);
%		try
			obj = loadData(obj,vars(v),file);
			if strcmp(vars{v},'FG_N2O')
				obj.data.avg.FG_N2O.data = -obj.data.avg.FG_N2O.data; % loss = positive
			end
			obj = loadDiag(obj,vars(v),0);
			% Diagnostic loop (sometimes more than 1 available)
			for d = 1:length(obj.diag.(vars{v}))
				close all
				% Extract data
				romsdat    = nanmean(obj.data.avg.(vars{v}).data,3);
				if strcmp(vars{v},'sustr');
					romsdat = romsObj.u2rho(romsdat);
				elseif strcmp(vars{v},'svstr');
					romsdat = romsObj.v2rho(romsdat);
				end
				diagdat    = nanmean(obj.diag.(vars{v})(d).slice,3);
				% Get levels
				levs  = linspace(lims{v}(1),lims{v}(2),40);
				dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
				% Make figures
				[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
				% ROMS figure
				set(0,'CurrentFigure',figs(1));
				title(['ROMS ',obj.data.avg.(vars{v}).name],'Interpreter','Latex');
				ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
				export_fig('-png',[obj.paths.plots.diag,vars{v},'_roms'],'-m5');
				close(figs(1));
				% Diag figure
				set(0,'CurrentFigure',figs(2));
				title([obj.diag.(vars{v})(d).name],'Interpreter','Latex');
				ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
				export_fig('-png',[obj.paths.plots.diag,vars{v},'_diag',num2str(d)],'-m5');
				close(figs(2));   
				% Diff figure
				set(0,'CurrentFigure',figs(3));
				title(['Difference'],'Interpreter','Latex');
				ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
				export_fig('-png',[obj.paths.plots.diag,vars{v},'_diff',num2str(d)],'-m5');
				close(figs(3));
			end
%		catch
 %            disp(['...skipping ',vars{v},'...']);
  %           continue
   %     end
    end
    % Clear data
    obj = clearROMS(obj);
    clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface chla
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
	try
		obj = loadData(obj,vars,file);
		obj = loadDiag(obj,vars,0);
        for d = 1:length(obj.diag.(vars{1}))
            % Reduce data
            romsdat = nanmean(obj.data.avg.(vars{1}).data,3);
            diagdat = nanmean(obj.diag.(vars{1})(d).slice,3);
			romsdat(romsdat<0) = 0;
            diffdat = romsdat - diagdat;
			romsdat = sqrt(romsdat);
			diagdat = sqrt(diagdat);
            % Plot
            [figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
            % ROMS figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS ',obj.data.avg.(vars{1}).name,': sfc'],'Interpreter','Latex');
            ylabel(cbs(1),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
			cbs(1).XTick      = absLevs([1:10:101]);
            cbs(1).XTickLabel = absLbls([1:10:101]); 
            export_fig('-png',[obj.paths.plots.diag,vars{1},'_roms'],'-m5');
            close(figs(1));
            % Diag figure
            set(0,'CurrentFigure',figs(2));
            title([obj.diag.(vars{1}).name,': sfc'],'Interpreter','Latex');
            ylabel(cbs(2),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
			cbs(2).XTick      = absLevs([1:10:101]);
            cbs(2).XTickLabel = absLbls([1:10:101]); 
            export_fig('-png',[obj.paths.plots.diag,vars{1},'_diag',num2str(d)],'-m5');
            close(figs(2));
            % Diff figure
            set(0,'CurrentFigure',figs(3));
            title(['Difference'],'Interpreter','Latex');
            ylabel(cbs(3),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{1},'_diff',num2str(d)],'-m5');
            close(figs(3));
		end
		% Clear data
		obj = clearROMS(obj);
		clearvars -except obj plots pltcnt file all_vars
	catch
		 disp(['...skipping ',vars{1},'...']);
	end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMZ thickness
% P6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
   % Get OMZ thickness
    if ~opt(1) | ~ismember('O2',all_vars);
        disp(['...skipping ',vars{1},'... (no variable)']);
    else
		try
			obj = computeVar(obj,{'OMZ'},file,'thresh',omzthresh);
		catch
             disp(['...skipping OMZ thickness...']);
        end	
        obj = loadDiag(obj,{'OMZ'},omzthresh);
        % Make comparison plots
        for d = 1:length(obj.diag.OMZ);
            for th = 1:length(omzthresh)
                levs  = linspace(lims{th}(1),lims{th}(2),length(cmap)+1);
                dlevs = linspace(dlims{th}(1),dlims{th}(2),length(cmap)+1);
                if ndims(obj.data.avg.OMZ.int)==4
                    romsdat = nanmean(squeeze(obj.data.avg.OMZ.int(:,:,th,:)),3) .* obj.grid.mask_rho;
                else
                    romsdat = squeeze(obj.data.avg.OMZ.int(:,:,th)) .* obj.grid.mask_rho;
                end
                diagdat = squeeze(obj.diag.OMZ(d).slice(:,:,th));
                [figs,cbs] = mapCmp(obj,romsdat,diagdat,'levels',levs,'difflevels',dlevs);
                % ROMS figure
                set(0,'CurrentFigure',figs(1));
                title(['ROMS: ',obj.data.avg.OMZ.name,'(O$_2$ $<$ ',num2str(omzthresh(th)),' mmol $m^{-3}$)'],'Interpreter','Latex');
                ylabel(cbs(1),obj.data.avg.OMZ.units,'Interpreter','Latex')
                set(gcf,'ColorMap',cmap);
                export_fig('-png',[obj.paths.plots.diag,'OMZ_roms_th',num2str(omzthresh(th))],'-m5');
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                title([obj.diag.OMZ(d).name],'Interpreter','Latex');
                ylabel(cbs(2),obj.data.avg.OMZ.units,'Interpreter','Latex')
                set(gcf,'ColorMap',cmap);
                export_fig('-png',[obj.paths.plots.diag,'OMZ_diag_',num2str(d),'_th',num2str(omzthresh(th))],'-m5');
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['Difference'],'Interpreter','Latex');
                ylabel(cbs(3),obj.data.avg.OMZ.units,'Interpreter','Latex')
                export_fig('-png',[obj.paths.plots.diag,'OMZ_diff_',num2str(d),'_th',num2str(omzthresh(th))],'-m5');
                close all
            end
        end
    end
    % Clear data
    obj = clearROMS(obj);
    clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transects on map
% P7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Get lon lines
    for i = 1:length(lons)
        lony{i} = [-90:0.1:90];
        lonx{i} = [lons(i)*ones(size(lony{i}))];
    end
    % Get lat lines
    for i = 1:length(lats)
        latx{i} = [0:0.1:360];
        laty{i} = [lats(i)*ones(size(latx{i}))];
    end
    % Plot
    [fig] = quickMap(obj,'ticks',1,'fontsize',8);
    hold on
    m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'k','linewidth',2);
    for i = 1:length(lons)
        [in,~] = inpolygon(lonx{i},lony{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
        m_plot(lonx{i}(in==1),lony{i}(in==1),'--k');
    end
    for i = 1:length(lats)
        [in,~] = inpolygon(latx{i},laty{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
        m_plot(latx{i}(in==1),laty{i}(in==1),'--k');
    end
    title(['Location of Transects'],'Interpreter','Latex');
    export_fig('-png',[obj.paths.plots.diag,'trans_locations'],'-m5');
    % Clear data
    obj = clearROMS(obj);
    clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps of POC_FLUX_IN
% P8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));

   % Get POC FLUX IN 
    if ~opt(1) | ~ismember('POC_FLUX_IN',all_vars);
        disp(['...skipping POC_FLUX_IN... (no variable)']);
    else
        % Load POC FLUX IN
        obj = clearROMS(obj);
        obj = zslice(obj,{'POC_FLUX_IN'},75,file);
        obj = loadDiag(obj,{'POC_FLUX_IN'},75);

        % Depth loop
        for d = 1:length(obj.diag.POC_FLUX_IN); % skip 100m estimate from Clements
            close all
            romsdat    = nanmean(squeeze(obj.data.avg.POC_FLUX_IN.slice),3);
            diagdat    = nanmean(obj.diag.POC_FLUX_IN(d).slice,3);
            [figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{1},'levels',levs{1},'difflevels',dlevs{1});
            % ROMS figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS ',obj.data.avg.POC_FLUX_IN.name,': 75m'],'Interpreter','Latex');
            ylabel(cbs(1),obj.data.avg.POC_FLUX_IN.units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_roms'],'-m5');
            close(figs(1));
            % Diag figure
            set(0,'CurrentFigure',figs(2));
            title([obj.diag.POC_FLUX_IN(d).name,': Euphotic'],'Interpreter','Latex');
            ylabel(cbs(2),obj.data.avg.POC_FLUX_IN.units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_diag',num2str(d)],'-m5');
            close(figs(2));
            % Diff figure
            set(0,'CurrentFigure',figs(3));
            title(['Difference'],'Interpreter','Latex');
            ylabel(cbs(3),obj.data.avg.POC_FLUX_IN.units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_diff',num2str(d)],'-m5');
        end
    end
    % Clear data
    obj = clearROMS(obj);
    clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gridded obs vs ROMS (including ratios)
% P9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Call extract_3d_obs
    extract_3d_obs(obj,xlims,ylims);
    % Call extract_3d_roms
    extract_3d_roms(obj,xlims,ylims,file);
    % Clear data
    obj = clearROMS(obj);
    clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D deep diagnostics (zslices)
% P10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
   % Variable loop
   for v = 1:length(vars)
      if ~opt(v) | ~ismember(vars{v},all_vars);
         disp(['...skipping ',vars{v},'...']);
         continue
      end
      % Get data
      obj = clearROMS(obj);
	  try	
		  obj = zslice(obj,vars(v),zdeps,file);
      catch
          disp(['...skipping OMZ thickness...']);
          continue
      end		  
      obj = loadDiag(obj,vars(v),zdeps);
      % Depth loop
      for z = 1:length(zdeps)
         for d = 1:length(obj.diag.(vars{v}))
            % Get levels
            levs  = linspace(lims{v,z}(1),lims{v,z}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42); 
            % Extract data
            romsdat    = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,z,:)),3);
            diagdat    = nanmean(squeeze(obj.diag.(vars{v})(d).slice(:,:,z,:)),3);
            close all
            % Make figures
            [figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
            % ROMS figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_roms'],'-m5');
            close(figs(1));
            % Diag figure
            set(0,'CurrentFigure',figs(2));
            title([obj.diag.(vars{v})(d).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diag',num2str(d)],'-m5');
            close(figs(2));   
            % Diff figure
            set(0,'CurrentFigure',figs(3));
            title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diff',num2str(d)],'-m5');   
         end
      end
   end
   % Clear data
   obj = clearROMS(obj);
   clearvars -except obj plots pltcnt file all_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe comparisons against GEOTRACES 
% P11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));

    % Go through each transect
    fdir = ['/data/project1/demccoy/ROMS/validation/Fe/'];
    for i = 1:length(sections)
        % Load section
        this_section = load([fdir,sections{i},'.mat']);

        % Sort profiles by increasing lon or lat
        if strcmp(coord{i},'lon');
            [B,I] = sort(this_section.lat0);
        elseif strcmp(coord{i},'lat');
            [B,I] = sort(this_section.log0);
        end

        % Remake section
        obs.lat0 = this_section.lat0(I);
        obs.log0 = this_section.log0(I);
        obs.bot0 = this_section.bot0(I);
        obs.dep0 = this_section.dep0(I,:);
        obs.Fe0  = this_section.Fe0(I,:);
        
        this_section = obs;
        
        % Call getProfile at each station
        obj = getProfile(obj,{'Fe'},this_section.log0,this_section.lat0,file);
        % Cycle through each station and interpolate ROMS data
        tmp = nan([size(this_section.Fe0) 12]);
        for j = 1:length(this_section.log0)
            % Get good points
            this_depth = this_section.dep0(j,:);
            idx = ~isnan(this_depth);
            good_depths = this_depth(idx);
            for t = 1:size(obj.profile.depth,3)
                % Get data for interpolation
                tmpdep       = flipud(squeeze(obj.profile.depth(j,:,t))');
                tmpdat       = flipud(squeeze(obj.data.avg.Fe.profile(j,:,t))');
                % Override surface point to be 0
                tmpdep(1)  = 0;
                tmpout       = interp1(tmpdep,tmpdat,-good_depths)';
                tmp(j,find(idx==1),t) = tmpout; 
            end    
        end
        % Get annual average
        out{i}.roms = tmp; 
        out{i}.diag = repmat(this_section.Fe0 ./ 1000, [1 1 12]); %nmol/L to mmol/m3
        out{i}.dep  = repmat(this_section.dep0, [1 1 12]);
        out{i}.lon  = repmat(obj.profile.lon',[1 size(out{i}.dep,2) 12]);
        out{i}.lat  = repmat(obj.profile.lat',[1 size(out{i}.dep,2) 12]);;
        if i < length(sections)
            obj = clearROMS(obj);
        end
    end

    % Make scatter figures
    for i = 1:length(sections)    
        if strcmp(coord{i},'lon');
            deg = nanmean(out{i}.lat,3);
        elseif strcmp(coord{i},'lat');
            deg = nanmean(out{i}.lon,3);
        end    
        dep     = nanmean(out{i}.dep,3);
        romsdat = nanmean(out{i}.roms,3);
        diagdat = nanmean(out{i}.diag,3);
        levs  = linspace(lims{1}(1),lims{1}(2),40);
        dlevs = linspace(dlims{1}(1),dlims{1}(2),41);

        % ROMS
        fig = romsObj.piofigs('mfig',0.5);
        C = contourf(deg,dep,romsdat,levs,'linestyle','none');    
        hold on
        S = scatter(deg(:),dep(:),10,diagdat(:),'filled');
        caxis([levs(1) levs(end)]);
        set(gca,'YDir','Reverse');
        yl = max(dep(~isnan(romsdat)));
        ylim([0 yl]);
        ylabel('Depth ($m$)','Interpreter','Latex');
        cb = colorbar;
        ylabel(cb,obj.data.avg.Fe.units,'Interpreter','Latex');
        xlbl = get(gca,'XTickLabel');
        if strcmp(coord{i},'lat');
            for j = 1:length(xlbl)
                if str2num(xlbl{j}) > 180
                    newlbl{j} = [num2str(-(str2num(xlbl{j}) - 360)),char(176),'W'];
                else
                    newlbl{j} = [num2str((str2num(xlbl{j}))),char(176),'E'];
                end
            end
        else
            for j = 1:length(xlbl)
                if str2num(xlbl{j})>=0
                    newlbl{j} = [num2str(str2num(xlbl{j})),char(176),'N'];
                else
                    newlbl{j} = [num2str(abs(str2num(xlbl{j}))),char(176),'S'];
                end
            end
        end
        set(gca,'XTickLabel',newlbl);
        set(gca,'Color',rgb('DimGray'));
        set(gcf,'inverthardcopy','off');
        title(['ROMS ',obj.data.avg.Fe.name,': ',sections{i}],'Interpreter','Latex');
        export_fig('-png',[obj.paths.plots.diag,'Fe_',sections{i},'_roms_v_obs'],'-m5');
        close all
    end
   % Clear data
   obj = clearROMS(obj);
   clearvars -except obj plots pltcnt file all_vars
end
