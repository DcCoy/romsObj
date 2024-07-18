function obj = romsComp(obj,file,plotchoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to automatically generate comparison plots between ROMS simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% - obj = romsComp(obj,file1,file2,plotchoice)
%
% Inputs:
% - obj   = roms objects (2) that are already initialized
%			... comparisons will be obj(1) - obj(2)
% - file  = cell array of files to load  
% -         ... file{1} = obj(1) files etc.
% - plotchoice (see below) 
%
%	------------------------------
%	---- STANDARD DIAGNOSTICS ----
%	------------------------------
%
%   1  == depth sections
%   2  == zonal transects
%   3  == meridional transects
%   4  == surface comparisons
%
%    -----------------------
%    -- OTHER DIAGNOSTICS --
%    -----------------------
%
%    5  == surface chlA map (log scale)
%    6  == OMZ thickness     (0, 5, 10, 20, 50uM)
%
% Example:
% - obj = romsComp(obj,{[1],[1]},[2:11]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle
addpath /data/project1/demccoy/ROMS/validation/n2o

% Clear objects
for i = 1:length(obj)
	obj(i) = clearROMS(obj(i));
end

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
   plots(pltcnt).opt   = [   1         1          1            1            1            1           1             1       ];
   plots(pltcnt).vars  = { 'MLD'  ,  'SSH'   ,  'NPP'   ,   'FG_N2O' ,   'sustr'  ,   'svstr'  ,    'ws'    ,    'wsc'     };
   plots(pltcnt).cmaps = { 'deep' ;  'deep'  ; 'algae'  ;  'balance' ;  'balance' ;  'balance' ;   'amp'    ; 'balance'    };
   plots(pltcnt).lims  = {[0  150];[0    1.5];[0   1000];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[0     1e-4];[-2e-10 2e-10]};
   plots(pltcnt).dlims = {[-60 60];[-0.3 0.3];[-700 700];[-5e-7 5e-7];[-1e-4 1e-4];[-1e-4 1e-4];[-1e-4 1e-4];[-2e-10 2e-10]};

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
% NOTE: Diagnostics are hard coded to omzthresh
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change to diagDir and load eoverrides options
compDir = ['/data/project1/demccoy/ROMS/',obj(1).info.simName,'/analysis/comp/'];
mkdir(compDir)
cd(compDir);
try; run(['comp_overrides.m']);
catch; disp('No overrides'); 
end
comppath = [compDir,'plots/',obj(1).info.runName,'_VS_',obj(2).info.runName,'/'];
mkdir(comppath)

% Clear workspace and begin
clearvars -except obj plots comppath file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STANDARD DIAGNOSTICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = 0;

% Get variables list
roms1_vars = {obj(1).info.phy_avg(file{1}(1)).Variables.Name ...
			  obj(1).info.bgc_avg(file{1}(1)).Variables.Name ...
			  obj(1).info.dia_avg(file{1}(1)).Variables.Name};
roms2_vars = {obj(2).info.phy_avg(file{2}(1)).Variables.Name ...
			  obj(2).info.bgc_avg(file{2}(1)).Variables.Name ...
			  obj(2).info.dia_avg(file{2}(1)).Variables.Name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D diagnostics (zslices)
% P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	% Inter-ROMS comparisons
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		try
			for i = 1:length(obj)
				obj(i) = clearROMS(obj(i));
				obj(i) = zslice(obj(i),vars(v),zdeps,file{i});
				tmp{i} = obj(i).data.avg.(vars{v}).slice;
			end
		catch
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for z = 1:length(zdeps)
			% Get levels
            levs  = linspace(lims{v,z}(1),lims{v,z}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
			close all
			roms1dat    = nanmean(squeeze(tmp{1}(:,:,z,:)),3);
			roms2dat    = nanmean(squeeze(tmp{2}(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runName,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms1'],'-m5');
			close(figs(1));
			% ROMS2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runName,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms2'],'-m5');
			close(figs(2));
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_z',num2str(zdeps(z)),'_roms_diff'],'-m5');
			close(figs(3));
		end
	end
	% Clear data
	for i = 1:length(obj)
		obj(i) = clearROMS(obj(i));
	end
	clearvars -except obj plots comppath file pltcnt roms1_vars roms2_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zonal transects
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v) 
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for l = 1:length(lats)
			try
				for i = 1:length(obj)
					obj(i) = clearROMS(obj(i));
					obj(i) = sliceROMS(obj(i),vars(v),choice,lats(l),file{i},...
						'zdep',obj(i).grid.z_avg_dep(obj(i).grid.z_avg_dep<=zlims(end)));
					tmp{i} = obj(i).data.avg.(vars{v}).slice;
				end
			catch
				disp(['...skipping ',vars{v},'...']);
				continue
			end
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
			% Get levels
            levs  = linspace(lims{v}(1),lims{v}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
            % Make figures
            [figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
                'figdim',0.5,'levels',levs,'difflevels',dlevs);
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runName,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms1'],'-m5');
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runName,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms2'],'-m5');
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lat',num2str(lats(l)),'_roms_diff'],'-m5');
			close(figs(3))
		end
	end
	% Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
	clearvars -except obj plots comppath file pltcnt roms1_vars roms2_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% longitude transects
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
    % Inter-ROMS comparisons
    for v = 1:length(vars)
        if ~opt(v)
            disp(['...skipping ',vars{v},'...']);
            continue
        end
		for l = 1:length(lons)
			% Get data
			try
				for i = 1:length(obj)
					obj(i) = clearROMS(obj(i));
					obj(i) = sliceROMS(obj(i),vars(v),choice,lons(l),file{i},...
						'zdep',obj(i).grid.z_avg_dep(obj(i).grid.z_avg_dep<=zlims(end)));
					tmp{i} = obj(i).data.avg.(vars{v}).slice;
				end
			catch
				disp(['...skipping ',vars{v},'...']);
				continue
			end
			% Get levels
			levs  = linspace(lims{v}(1),lims{v}(2),40);
			dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
			roms1dat = tmp{1}; 
			roms2dat = tmp{2}; 
            % Make figures
            [figs,cbs] = sliceCmp(obj(1),roms1dat,roms2dat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
               'figdim',0.5,'levels',levs,'difflevels',dlevs);
			% ROMS1 figure
			set(0,'CurrentFigure',figs(1));
			title([obj(1).info.runName,' ',obj(1).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms1'],'-m5');
			close(figs(1));
			% ROM2 figure
			set(0,'CurrentFigure',figs(2));
			title([obj(2).info.runName,' ',obj(2).data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms2'],'-m5');
			close(figs(2));
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['ROMS Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[comppath,vars{v},'_lon',num2str(lons(l)),'_roms_diff'],'-m5');
			close(figs(3))
		end
	end
	% Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
	clearvars -except obj plots comppath file pltcnt roms1_vars roms2_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface diagnostics
% P4
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
		try
			for i = 1:length(obj)
				obj(i) = clearROMS(obj(i));
				obj(i) = loadData(obj(i),vars(v),file{i});
				if strcmp(vars{v},'FG_N2O')
					obj(i).data.avg.FG_N2O.data = -obj(i).data.avg.FG_N2O.data; % loss = positive
				end		
			end
		catch
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		% Get levels
		levs  = linspace(lims{v}(1),lims{v}(2),40);
		dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
		% Extract data
		roms1dat = nanmean(obj(1).data.avg.(vars{v}).data,3);
		roms2dat = nanmean(obj(2).data.avg.(vars{v}).data,3);
		% Make figures
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
	    % ROMS figure
		set(0,'CurrentFigure',figs(1));
		title([obj(1).info.runName,' ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms1'],'-m5');
		close(figs(1));
		set(0,'CurrentFigure',figs(2));
		title([obj(2).info.runName,' ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms2'],'-m5');
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms_diff'],'-m5');
		close(figs(3));
	end
    % Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
	clearvars -except obj plots comppath file pltcnt roms1_vars roms2_vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface chla
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	try
		for i = 1:length(obj)
			obj(i) = clearROMS(obj(i));
			obj(i) = loadData(obj(i),vars(1),file{i});
		end
		close all
		roms1dat = nanmean(obj(1).data.avg.(vars{v}).data,3);
		roms2dat = nanmean(obj(2).data.avg.(vars{v}).data,3);
		roms1dat(roms1dat<0) = 0;
		roms2dat(roms2dat<0) = 0;
		roms1dat = sqrt(roms1dat);
		roms2dat = sqrt(roms2dat);
		% Plot
		[figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title([obj(1).info.runName,' ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		cbs(1).XTick      = absLevs([1:10:101]);
		cbs(1).XTickLabel = absLbls([1:10:101]); 
		export_fig('-png',[comppath,vars{v},'_roms1'],'-m5');
		close(figs(1));
		set(0,'CurrentFigure',figs(2));
		title([obj(2).info.runName,' ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
		ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
		cbs(2).XTick      = absLevs([1:10:101]);
		cbs(2).XTickLabel = absLbls([1:10:101]); 
		export_fig('-png',[comppath,vars{v},'_roms2'],'-m5');
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['ROMS Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
		export_fig('-png',[comppath,vars{v},'_roms_diff'],'-m5');
		close(figs(3));
		% Clear data
		for i = 1:length(obj)
			obj(i) = clearROMS(obj(i));
		end
		clearvars -except obj plots comppath file pltcnt roms1_vars roms2_vars
	catch
		 disp(['...skipping ',vars{1},'...']);
	end	
end
