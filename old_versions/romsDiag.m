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
% - file = file(s) to diagnose, see obj.info.avg.files  
% - plotchoice (see below) 
%
%	------------------------------');
%	---- PHYSICAL DIAGNOSTICS ----');
%	------------------------------');
%	1  == 3d physical diagnostics');
%	2  == equator physical slices');
%	3  == surface physical fields');
%	4  == longitude physical slices');
%	5  == u velocity at equator');
%	-------------------------');
%	---- BGC DIAGNOSTICS ----');
%	-------------------------');
%	6  == 3d bgc diagnostics');
%	7  == equator bgc slices');
%	8  == surface bgc fields');
%	9  == longitude bgc slices');
%	10 == surface chlA');
%	11 == OMZ thickness');
%	12 == N-cycle profile comparisons');
%	13 == ROMS vs OCIM (N2O)
%	14 == Gridded obs vs ROMS (including ratios) 
%	15 == Fe depth slices vs Tagliabue
%	16 == Fe transects vs GEOTRACES
%	17 == POC flux comparisons at 75m
%	------------------------------');
%	------ OTHER DIAGNOSTICS -----');
%	------------------------------');
%	18 == slice degree maps');
%
% Example for year 10:
% - obj = romsDiag(obj,10,[1:11]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all
addpath /data/project1/demccoy/ROMS/
addpath /data/project1/demccoy/ROMS/validation/ncycle
addpath /data/project1/demccoy/ROMS/validation/n2o

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

% PHYSICAL PLOT options
% (1) PHYS: 3D zslices
pltcnt = 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1];
	plots(pltcnt).vars  = {'temp','salt','sigma'};
	plots(pltcnt).cmaps = {'thermal','haline','dense'};
	plots(pltcnt).zdeps = [0 50 150 300];
	plots(pltcnt).levs  = {linspace(3,30,40),linspace(3,30,40),linspace(0,25,40),linspace(2,18,40);...
						   linspace(30,37,40),linspace(30,37,40),linspace(33.4,36.5,40),linspace(33.5,35.5,40);...
						   linspace(19,27,40),linspace(19,27,40),linspace(24,28,40),linspace(26.5,28.6,40)};
	plots(pltcnt).dlevs = {linspace(-3,3,40),linspace(-0.6,0.6,40),linspace(-0.6,0.6,40)};

% (2) PHYS: Equator slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1];
	plots(pltcnt).vars   = {'temp','salt','sigma'};
	plots(pltcnt).cmaps  = {'thermal','haline','dense'};
	plots(pltcnt).choice = 'lat';
	plots(pltcnt).lats   = 0;
	plots(pltcnt).xlims  = [140 260];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).levs   = {linspace(3,30,40),linspace(34.5,35.5,40),linspace(22,32,40)};
	plots(pltcnt).dlevs  = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (3) PHYS: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1];
	plots(pltcnt).vars  = {'MLD','SSH'};
	plots(pltcnt).cmaps = {'deep','deep'};
	plots(pltcnt).bal   = [2 2];
	plots(pltcnt).levs  = {linspace(0,150,40),linspace(0,1.5,40)};
	plots(pltcnt).dlevs = {linspace(-60,60,40),linspace(-0.3,0.3,40)};

% (4) PHYS: Longitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1];
	plots(pltcnt).vars   = {'temp','salt','sigma'};
	plots(pltcnt).cmaps  = {'thermal','haline','dense'};
	plots(pltcnt).choice = 'lon';
	plots(pltcnt).lons   = [210 255 272];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).xlims  = [-30 50;-40 20;-40 10];
	plots(pltcnt).levs   = {linspace(3,30,40),linspace(3,30,40),linspace(3,30,40);...
						    linspace(33,36.5,40),linspace(33.5,36.5,40),linspace(33.5,36.5,40);...
						    linspace(22,32,40),linspace(22,32,40),linspace(22,32,40)};
	plots(pltcnt).dlevs  = {linspace(-3,3,40),linspace(-0.7,0.7,40),linspace(-0.7,0.7,40)};

% (5) PHYS: Zonal equator velocity
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).zlims = [50 500];
	plots(pltcnt).xlims = [-20 20];
	plots(pltcnt).vars  = {'u'};
	plots(pltcnt).cmaps = {'balance'};
	plots(pltcnt).levs  = {linspace(-0.3,0.3,40)};
	plots(pltcnt).dlevs = {linspace(-0.3,0.3,40)};

% BIOGEOCHEMICAL PLOT options
% (6) BGC: 3D zslices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1 1 1 1 1 1 1];
	plots(pltcnt).vars  = {'O2','NOX','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).zdeps = [0 50 150 300 450];
	plots(pltcnt).bal   = [2 2 2 2 2 2];
    plots(pltcnt).levs  = {linspace(0,350,40),linspace(0,350,40),linspace(0,350,40),linspace(0,350,40),linspace(0,350,40);...
                           linspace(0,40,40),linspace(0,40,40),linspace(0,40,40),linspace(0,40,40),linspace(0,40,40);...
                           linspace(0,3,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
                           linspace(0,1,40),linspace(0,1,40),linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
                           linspace(0,0.02,40),linspace(0,0.05,40),linspace(0,0.08,40),linspace(0,0.08,40),linspace(0,0.08,40);...
						   linspace(0,3,40),linspace(0,3,40),linspace(0,1,40),linspace(0,1,40),linspace(0,1,40);...
						   linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40),linspace(-25,25,40)};
    plots(pltcnt).dlevs = {linspace(-80,80,41),linspace(-10,10,41),linspace(-1,1,41),...
                           linspace(-0.6,0.6,41),linspace(-0.05,0.05,41),linspace(-0.6,0.6,41),linspace(-10,10,41)};

% (7) BGC: Equator slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1 1 1 1 1];
	plots(pltcnt).vars   = {'O2','NOX','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps  = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).choice = 'lat';
	plots(pltcnt).lats   = 0;
	plots(pltcnt).xlims  = [140 260];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).bal    = [2 2 2 2 2 2];
	plots(pltcnt).levs   = {linspace(0,300,40);...
						    linspace(0,45,40);...
						    linspace(0,3,40);...
						    linspace(0,1,40);...
						    linspace(0,0.08,40);...
						    linspace(0,3,40);...
						    linspace(-25,25,41)};
	plots(pltcnt).dlevs  = {linspace(-80,80,41),linspace(-15,15,41),linspace(-0.5,0.5,41),...
						    linspace(-1,1,41),linspace(-0.05,0.05,41),linspace(-1,1,41),linspace(-10,10,41)};

% (8) BGC: Surface fields
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1];
	plots(pltcnt).bal   = [2];
	plots(pltcnt).vars  = {'NPP'};
	plots(pltcnt).cmaps = {'algae'};
	plots(pltcnt).levs  = {linspace(0,1000,40)};
	plots(pltcnt).dlevs = {linspace(-700,700,40)};

% (9) BGC: Longitude slices
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt    = [1 1 1 1 1 1 1];
	plots(pltcnt).vars   = {'O2','NOX','PO4','NO2','N2O','NH4','nstar'};
	plots(pltcnt).cmaps  = {'-ice','tempo','tempo','tempo','tempo','tempo','-balance'};
	plots(pltcnt).choice = 'lon';
	plots(pltcnt).lons   = [210 255 272];
	plots(pltcnt).zlims  = [0 1000];
	plots(pltcnt).xlims  = [-30 50;-40 20;-40 10];
	plots(pltcnt).bal    = [2 2 2 2 2 2];
	plots(pltcnt).levs   = {linspace(0,300,40),linspace(0,300,40),linspace(0,300,40);...
						    linspace(0,45,40),linspace(0,45,40),linspace(0,45,40);...
						    linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						    linspace(0,0.4,40),linspace(0,0.4,40),linspace(0,0.4,40);...
						    linspace(0,0.08,40),linspace(0,0.08,40),linspace(0,0.08,40);...
						    linspace(0,3,40),linspace(0,3,40),linspace(0,3,40);...
						    linspace(-25,25,41),linspace(-25,25,41),linspace(-25,25,41)};
	plots(pltcnt).dlevs  = {linspace(-80,80,41),linspace(-15,15,41),linspace(-0.5,0.5,41),...
						    linspace(-1,1,41),linspace(-0.05,0.05,41),linspace(-1,1,41),linspace(-10,10,41)};

% (10) BGC: surface chlA
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).vars      = {'SFC_CHL'};
	plots(pltcnt).cmaps     = {'algae'};
	plots(pltcnt).absLevs   = log10([0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0 6.0 10.0]);
	plots(pltcnt).absLbls   = [0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0 6.0 10.0];
	plots(pltcnt).absCaxis  = real(log10([0.01 10]));
	plots(pltcnt).diffLevs  = [0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0];
	plots(pltcnt).diffLbls  = [0.005 0.01 0.02 0.05 0.1 0.2 0.4 0.8 1.5 3.0];
	plots(pltcnt).diffLevs  = [-fliplr(plots(pltcnt).diffLevs) plots(pltcnt).diffLevs];
	plots(pltcnt).diffLbls  = [-fliplr(plots(pltcnt).diffLbls) plots(pltcnt).diffLbls];
	plots(pltcnt).diffLevs  = romsMaster.dfloglevs(plots(pltcnt).diffLevs,0.001);
	plots(pltcnt).diffCaxis = [plots(pltcnt).diffLevs(1) plots(pltcnt).diffLevs(end)];

% (11) BGC: OMZ thickness
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).cmap      = cbrewer('seq','YlOrRd',40);
	plots(pltcnt).cmap(1,:) = [1 1 1];
	plots(pltcnt).omzthresh = [20 50];
	plots(pltcnt).levs      = {linspace(0,1000,40),linspace(0,2000,40)};
	plots(pltcnt).dlevs     = {linspace(-1000,1000,41),linspace(-2000,2000,41)};;


% (12) BGC: Ncycle profiles
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).vars         = {'NOX','NO2','N2O','NH4'};
	plots(pltcnt).region_name  = {     'Eq',            'ETSP',              'ETNP',         'CCS',           'NPSG'   }; 
	plots(pltcnt).obs_regions  = {[217 222 -2 3],[278.0 281.0 -15 -12],  [250 255 14 19],[237 242 30 35],[206 211 28 33]};
	plots(pltcnt).roms_regions = {[217 222 -2 3],[275.2 278.2 -6.7 -3.7],[255 260 13 18],[237 242 30 35],[206 211 28 33]};
	plots(pltcnt).tracer_lims  = {    [0 50],            [0 50],              [0 50],		[0 50],          [0 50];...    % NO3
									  [0 1.5],		     [0 7],				  [0 6],		[0 0.5],         [0 0.5];...   % NO2
									  [0 50],		     [0 160],			  [0 70],		[0 60],          [0 60];...    % N2O*1000
									  [0 1.5],		     [0 3.5],			  [0 1.5],		[0 2.5],         [0 0.5];...   % NH4
									  [0 250],		     [0 250],			  [0 220],		[0 300],         [0 300]};     % O2
	plots(pltcnt).ratio_lims   = {	  [0 1.5],			 [0 7],				  [0 6],        [0 0.5],         [0 0.5];...   % NO2
									  [0 1.5],			 [0 3.5],			  [0 1.5],      [0 2.5],         [0 0.5];...   % NH4
									  [1E-3 1E3] ,		 [1E-3 1E3],		  [1E-3 1E3],   [1E-3 1E3],      [1E-3 1E3];...% NO2vNH4
									  [1E-3 1E3] ,		 [1E-3 1E3],		  [1E-3 1E3],   [1E-3 1E3],      [1E-3 1E3];...% NH4vNO2
									  [0 210],			 [0 250],			  [0 220],      [0 300],         [0 300]};	   % O2

% (13) BGC: ROMS vs OCIM
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).regions = {[225 300 -30 30]};

% (14) BGC: Gridded obs vs ROMS (including ratios) 
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).xlims = [100 300];
	plots(pltcnt).ylims = [-40  60];

% (15) BGC: Fe depth slices vs Tagliabue 
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).opt   = [1];
	plots(pltcnt).vars  = {'Fe'};
	plots(pltcnt).cmaps = {'thermal'};
	plots(pltcnt).zdeps = [0 50 150 300 450 600 800 1000 1500 2000];
	plots(pltcnt).bal   = [2];
	plots(pltcnt).levs  = {linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),...
						   linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),...
						   linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40),linspace(0,1e-3,40)};
    plots(pltcnt).dlevs = {linspace(-1e-3,1e-3,41)};

% (16) BGC: Fe lon/lat slices vs GEOTRACES 
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	%plots(pltcnt).sections = {'GP02','GP13','GP16','GP18','GP19'};
	plots(pltcnt).sections = {'GP16'};
	plots(pltcnt).coord    = {'lat'};
	plots(pltcnt).cmaps    = {'thermal'};
	plots(pltcnt).levs     = {linspace(0,1e-3,40)};
	plots(pltcnt).dlevs    = {linspace(-1e-3,1e-3,41)};

% (17) POC_FLUX_IN comparisons 
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).cmaps = {'deep'};
	plots(pltcnt).levs  = {linspace(0,5e-4,40)};
	plots(pltcnt).dlevs = {linspace(-5e-4,5e-4,41)};

% (18) OTHER: Slice maps
pltcnt = pltcnt + 1;
plots(pltcnt).on = plotchoice(pltcnt);
	plots(pltcnt).lons = [210 255 272];
	plots(pltcnt).lats = [0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change to diagDir and load overrides options
diagDir = ['/data/project1/demccoy/ROMS/',obj.info.simName,'/analysis/diag/'];
mkdir(diagDir);
cd(diagDir);
try; run(['diag_overrides.m']);
catch; disp('No overrides');
end

% Reset pltcnt
pltcnt = 0;

% Clear workspace and begin
clearvars -except obj plots pltcnt file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PHYSICAL DIAGNOSTICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D physical diagnostics (zslices)
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
		obj = clearROMS(obj);
		obj = zslice(obj,vars(v),zdeps,file);
		obj = loadDiag(obj,vars(v),zdeps);
		% Depth loop
		for z = 1:length(zdeps)
			close all
			romsdat    = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,z,:)),3);
			diagdat    = nanmean(squeeze(obj.diag.(vars{v}).slice(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_roms'],'-m5');
			close(figs(1));
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diag'],'-m5');
			close(figs(2));	
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diff'],'-m5');	
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical equatorial section plots
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	for v = 1:length(vars);
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for l = 1:length(lats)
			obj = clearROMS(obj);
			obj = sliceROMS(obj,vars(v),choice,lats(l),file,...
				'zdep',obj.grid.z_avg_dep(obj.grid.z_avg_dep<zlims(end)));
			obj = sliceDiag(obj,vars(v),choice,lats(l),'zlim',zlims(end));
			romsdat = obj.data.avg.(vars{v}).slice;
			diagdat = obj.diag.(vars{v}).slice;
			[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_roms'],'-m5');
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diag'],'-m5');
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diff'],'-m5');
			close(figs(3))
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface physical diagnostics
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	obj = clearROMS(obj);
	obj = loadData(obj,vars(opt==1),file);
	obj = loadDiag(obj,vars(opt==1),0);
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for d = 1:length(obj.diag.(vars{v}))
			close all
			romsdat    = nanmean(obj.data.avg.(vars{v}).data,3);
			diagdat    = nanmean(obj.diag.(vars{v})(d).slice,3);
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'bal',bal(v),'levels',levs{v},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_roms'],'-m5');
			close(figs(1));
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v})(d).name],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v})(d).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_diag_',num2str(d)],'-m5');
			close(figs(2));	
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_diff_',num2str(d)],'-m5');
			close(figs(3));
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical longitude section plots
% P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	for v = 1:length(vars);
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for l = 1:length(lons)
			obj = clearROMS(obj);
			obj = sliceROMS(obj,vars(v),choice,lons(l),file,...
				'zdep',obj.grid.z_avg_dep(obj.grid.z_avg_dep<zlims(end)));
			obj = sliceDiag(obj,vars(v),choice,lons(l),'zlim',zlims(end));
			romsdat = obj.data.avg.(vars{v}).slice;
			diagdat = obj.diag.(vars{v}).slice;
			[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_roms'],'-m5');
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diag'],'-m5');
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diff'],'-m5');
			close(figs(3))
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equator u-velocity slices 
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	obj = equatorUcmp(obj,'179E_160W');
	diagdat = obj.diag.u.slice;
	romsdat = nanmean(obj.data.avg.u.slice,3);
	[figs,cbs] = sliceCmp(obj,romsdat,diagdat,file,'zlims',zlims,'xlims',xlims,...
		'figdim',0.5,'cmap','balance','levels',levs{1},'difflevels',dlevs{1});
	% ROMS figure
	set(0,'CurrentFigure',figs(1));
	title(['ROMS ',obj.data.avg.(vars{1}).name,': ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(1),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[obj.paths.plots.diag,'equ_roms'],'-m5');
	close(figs(1))
	% Diag figure
	set(0,'CurrentFigure',figs(2));
	title([obj.diag.(vars{1}).name,': ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(2),obj.diag.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[obj.paths.plots.diag,'equ_diag'],'-m5');
	close(figs(2))
	% Difference figure
	set(0,'CurrentFigure',figs(3));
	title(['Difference: ',num2str(obj.slice.sect-360),'$^oW$'],'Interpreter','Latex');
	ylabel(cbs(3),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
	export_fig('-png',[obj.paths.plots.diag,'equ_diff'],'-m5');
	close(figs(3))
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% BGC DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3d bgc diagnostics
% P6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	% Go through each compare with diagnostics
	% Variable loop
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		obj = clearROMS(obj);
		if length(file) == 1
			obj = zslice(obj,vars(v),zdeps,file);
		else
			obj = getAvgData(obj,vars(v),file,'zlvl',zdeps);
		end
		obj = loadDiag(obj,vars(v),zdeps);
		% Depth loop
		for z = 1:length(zdeps)
			close all
			if length(file) == 1
				romsdat = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,z,:)),3);
			else
				romsdat = obj.data.avg.(vars{v}).slice;
			end
			diagdat    = nanmean(squeeze(obj.diag.(vars{v}).slice(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_roms'],'-m5');
			close(figs(1));
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diag'],'-m5');
			close(figs(2));	
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diff'],'-m5');	
		end
		obj = clearROMS(obj);
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc equatorial section plots
% P7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	for v = 1:length(vars);
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		for l = 1:length(lats)
			obj = clearROMS(obj);
			obj = sliceROMS(obj,vars(v),choice,lats(l),file,'type','z_avg','zlim',zlims(end));
			obj = sliceDiag(obj,vars(v),choice,lats(l),'zlim',zlims(end));
			romsdat = obj.data.avg.(vars{v}).slice;
			diagdat = obj.diag.(vars{v}).slice;
			if size(romsdat,3) ~= size(diagdat,3)
				romsdat = squeeze(romsdat(:,:,1,:));
				diagdat = squeeze(nanmean(diagdat,4));
				romsdat = repmat(romsdat,[1 1 1 12]);
				diagdat = repmat(diagdat,[1 1 1 12]);
			end
			[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_roms'],'-m5');
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diag'],'-m5');
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lat',num2str(lats(l)),'_diff'],'-m5');
			close(figs(3))
		end
		obj.data.avg = [];
		obj.diag = [];
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface field plots
% P8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	obj = clearROMS(obj);
	obj = loadData(obj,vars,file);
	obj = loadDiag(obj,vars,0);
	for v = 1:length(vars)
		for d = 1:length(obj.diag.(vars{v}))
			close all
			romsdat    = nanmean(obj.data.avg.(vars{v}).data,3);
			diagdat    = nanmean(obj.diag.(vars{v})(d).slice,3);
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs{v},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_roms'],'-m5');
			close(figs(1));
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v})(d).name],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v})(d).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_diag_',num2str(d)],'-m5');
			close(figs(2));	
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_diff_',num2str(d)],'-m5');
			close(figs(3));
		end
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc longitude section plots
% P9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	for v = 1:length(vars);
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		obj = clearROMS(obj);
		obj = sliceROMS(obj,vars(v),choice,lons,file,'type','z_avg','zlim',zlims(end));
		obj = sliceDiag(obj,vars(v),choice,lons,'zlim',zlims(end));
		romsdat = obj.data.avg.(vars{v}).slice;
		diagdat = obj.diag.(vars{v}).slice;
		if size(romsdat,3) ~= size(diagdat,3)
			romsdat = squeeze(romsdat(:,:,1,:));
			diagdat = squeeze(nanmean(diagdat,4));
			romsdat = repmat(romsdat,[1 1 1 12]);
			diagdat = repmat(diagdat,[1 1 1 12]);
		end
		for l = 1:length(lons)
			[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims,...
				'figdim',0.5,'slice',l,'levels',levs{v,l},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_roms'],'-m5');
			close(figs(1))
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diag'],'-m5');
			close(figs(2))
			% Difference figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(lons(l)-360),'$^oW$'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_lon',num2str(lons(l)),'_diff'],'-m5');
			close(figs(3))
		end
		obj.data.avg = [];
		obj.diag = [];
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc surface chla
% P10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	obj = loadData(obj,vars,file);
	obj = loadDiag(obj,vars,0);
	% Reduce data
	romsdat = nanmean(obj.data.avg.(vars{1})(1).data,3);
	diagdat = nanmean(obj.diag.(vars{1})(1).slice,3);
	diffdat = romsdat - diagdat;
	romsdat = real(log10(romsdat));
	diagdat = real(log10(diagdat));
	diffdat = romsMaster.dfloglevs(diffdat,0.01);
	% Plot
	[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{1},'levels',absLevs,'difflevels',diffLevs);
	% ROMS figure
	set(0,'CurrentFigure',figs(1));
	title(['ROMS ',obj.data.avg.(vars{1}).name,': sfc'],'Interpreter','Latex');
	ylabel(cbs(1),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
	cbs(1).XTickLabel = absLbls; 
	export_fig('-png',[obj.paths.plots.diag,vars{1},'_roms'],'-m5');
	close(figs(1));
	% Diag figure
	set(0,'CurrentFigure',figs(2));
	title([obj.diag.(vars{1}).name,': sfc'],'Interpreter','Latex');
	ylabel(cbs(2),obj.diag.(vars{1}).units,'Interpreter','Latex');
	cbs(2).XTickLabel = absLbls; 
	export_fig('-png',[obj.paths.plots.diag,vars{1},'_diag'],'-m5');
	close(figs(2));
	% Diff figure
	set(0,'CurrentFigure',figs(3));
	title(['Difference'],'Interpreter','Latex');
	ylabel(cbs(3),obj.data.avg.(vars{1}).units,'Interpreter','Latex');
	cbs(3).XTick = diffLevs;
	cbs(3).XTickLabel = diffLbls; 
	cbs(3).Limits = diffCaxis;
	export_fig('-png',[obj.paths.plots.diag,vars{1},'_diff'],'-m5');
	close(figs(3));
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMZ thickness
% P11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
    % Get OMZ thickness
	obj = computeVar(obj,{'OMZ'},file,'thresh',omzthresh);
	obj = OMZthick(obj,omzthresh);
    % Make comparison plots
	for d = 1:length(obj.diag.OMZ);
		for th = 1:length(omzthresh)
			if ndims(obj.data.avg.OMZ.int)==4
				romsdat = nanmean(squeeze(obj.data.avg.OMZ.int(:,:,th,:)),3);
			else
				romsdat = squeeze(obj.data.avg.OMZ.int(:,:,th));
			end
			diagdat = squeeze(obj.diag.OMZ(d).int(:,:,th));
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'levels',levs{th},'difflevels',dlevs{th});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS: ',obj.data.avg.OMZ.name,'($O_2$ $<$ ',num2str(omzthresh(th)),' $mmol$ $m^{-3}$)'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.OMZ.units,'Interpreter','Latex')
			set(gcf,'ColorMap',cmap);
			export_fig('-png',[obj.paths.plots.diag,'OMZ_roms_th',num2str(omzthresh(th))],'-m5');
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.OMZ(d).name],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.OMZ(d).units,'Interpreter','Latex')
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
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ncycle subsurface comparisons
% P12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	% Call data_locations
	if exist([obj.paths.plots.diag,'prof_data_locations.png'])==0
		disp('Plotting profile locations');
		fig = data_locations(obj,obs_regions,roms_regions);
		export_fig('-png',[obj.paths.plots.diag,'prof_data_locations'],'-m5');
		close all
	end
	% Call extract_obs	
    for i = 1:length(obs_regions)
		if exist(['data/obs_tracers_region_',region_name{i},'.mat'])==0
			mkdir data;
			disp('Extracting obs data');
			extract_obs_ncycle(obs_regions{i},region_name{i});
		end
    end
	% Call extract roms
	for i = 1:length(roms_regions)
		if (1)
			mkdir(['data/',obj.info.runName]);
			disp('Extracting ROMS data');
			extract_roms_ncycle(obj,vars,roms_regions{i},region_name{i},file);
		end
	end
	
	% Call obs_vs_roms
	obs_vs_roms
	% Call obs_vs_roms_ratios
	obs_vs_roms_ratios
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROMS vs OCIM N2O 
% P13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	% Call ocim_locations
	if exist([obj.paths.plots.diag,'ocim_locations.png'])==0
		fig = ocim_locations(obj,regions);
		export_fig('-png',[obj.paths.plots.diag,'ocim_locations'],'-m5');
		close all
	end
	% Extract OCIM	
	for i = 1:length(regions)
		OCIM_n2o_yield(obj,regions{i},i);
		ROMS_n2o_yield(obj,regions{i},i,file);
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gridded obs vs ROMS (including ratios)
% P14
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
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe depth slices 
% P15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));
	% Go through each compare with diagnostics
	% Variable loop
	for v = 1:length(vars)
		if ~opt(v)
			disp(['...skipping ',vars{v},'...']);
			continue
		end
		obj = clearROMS(obj);
		obj = zslice(obj,vars(v),zdeps,file);
		obj = loadDiag(obj,vars(v),zdeps);
		% Depth loop
		for z = 1:length(zdeps)
			close all
			romsdat    = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,z,:)),3);
			diagdat    = nanmean(squeeze(obj.diag.(vars{v}).slice(:,:,z,:)),3);
			[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs{v,z},'difflevels',dlevs{v});
			% ROMS figure
			set(0,'CurrentFigure',figs(1));
			title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_roms'],'-m5');
			close(figs(1));
			% Diag figure
			set(0,'CurrentFigure',figs(2));
			title([obj.diag.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(2),obj.diag.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diag'],'-m5');
			close(figs(2));	
			% Diff figure
			set(0,'CurrentFigure',figs(3));
			title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
			ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
			export_fig('-png',[obj.paths.plots.diag,vars{v},'_z',num2str(zdeps(z)),'_diff'],'-m5');	
		end
		obj.data.avg = [];
		obj.diag = [];
	end
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fe comparisons against GEOTRACES 
% P16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));

	% Go through each transect
	fdir = ['/data/project1/demccoy/ROMS/validation/Fe/'];
	for i = 1:length(sections)
		this_section = load([fdir,sections{i},'.mat']);
		% Call getProfile at each station
		obj = getProfile(obj,{'Fe'},this_section.log0,this_section.lat0,file);
		% Cycle through each station and interpolate ROMS data
		tmp = nan([size(this_section.Fe0) 12]);
		for j = 1:length(this_section.log0)
			this_depth = this_section.dep0(j,:);
			idx = ~isnan(this_depth);
			good_depths = this_depth(idx);
			for t = 1:size(obj.profile.depth,3)
				% Get data for interpolation
				tmpdep       = squeeze(obj.profile.depth(j,:,t));
				tmpdat       = squeeze(obj.data.avg.Fe.profile(j,:,t));
				% Override surface point to be 0
				tmpdep(end)  = 0;
				tmpout       = interp1(tmpdep,tmpdat,-good_depths);
				tmp(j,find(idx==1),t) = tmpout; 
			end	
		end
		% Get annual average
		out{i}.roms = tmp; 
		out{i}.diag = repmat(this_section.Fe0 ./ 1000, [1 1 12]); %nmol/L to mmol/m3
		out{i}.dep  = repmat(this_section.dep0, [1 1 12]);
		out{i}.lon  = obj.profile.lon;
		out{i}.lat  = obj.profile.lat;
	end

	% Make slice figures
	for i = 1:length(sections)
		obj.slice.depth = out{i}.dep;
		if strcmp(coord{i},'lon');
			obj.slice.coord = 'longitude';
			obj.slice.deg   = out{i}.lat';
		elseif strcmp(coord,'lat');
			obj.slice.coord = 'latitude';
			obj.slice.deg   = out{i}.lon';
		end
		obj.slice.deg   = repmat(obj.slice.deg,[1 size(obj.slice.depth,2)]);
		romsdat = out{i}.roms;	
		diagdat = out{i}.diag;
		[figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{1},'xlims',[min(obj.slice.deg(:)) max(obj.slice.deg(:))],'zlims',[0 5000],...
			'figdim',0.5,'slice',1,'levels',levs{1},'difflevels',dlevs{1});
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj.data.avg.Fe.name,': ',sections{i}],'Interpreter','Latex');
		ylabel(cbs(1),obj.data.avg.Fe.units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'Fe_',sections{i},'_roms'],'-m5');
		close(figs(1))
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title(['GEOTRACES: ',sections{i}],'Interpreter','Latex');
		ylabel(cbs(2),obj.data.avg.Fe.units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'Fe_',sections{i},'_diag'],'-m5');
		close(figs(2))
		% Difference figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference: ',sections{i}],'Interpreter','Latex');
		ylabel(cbs(3),obj.data.avg.Fe.units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'Fe_',sections{i},'_diff'],'-m5');
		close(figs(3))
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps of POC_FLUX_IN
% P17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
	unpackStruct(plots(pltcnt));

	% Load POC FLUX IN
	obj = clearROMS(obj);
	obj = zslice(obj,{'POC_FLUX_IN'},75,file);
	obj = loadDiag(obj,{'POC_FLUX_IN'},75);

	% Depth loop
	for d = 2:length(obj.diag.POC_FLUX_IN); % skip 100m estimate from Clements
		close all
		romsdat    = nanmean(obj.data.avg.POC_FLUX_IN.slice,3);
		diagdat    = nanmean(obj.diag.POC_FLUX_IN(d).slice,3);
		[figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{1},'levels',levs{1},'difflevels',dlevs{1});
		% ROMS figure
		set(0,'CurrentFigure',figs(1));
		title(['ROMS ',obj.data.avg.POC_FLUX_IN.name,': 75m'],'Interpreter','Latex');
		ylabel(cbs(1),obj.diag.POC_FLUX_IN(1).units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_roms'],'-m5');
		close(figs(1));
		% Diag figure
		set(0,'CurrentFigure',figs(2));
		title([obj.diag.POC_FLUX_IN(d).name,': Euphotic'],'Interpreter','Latex');
		ylabel(cbs(2),obj.diag.POC_FLUX_IN(1).units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_diag',num2str(d)],'-m5');
		close(figs(2));
		% Diff figure
		set(0,'CurrentFigure',figs(3));
		title(['Difference'],'Interpreter','Latex');
		ylabel(cbs(3),obj.diag.POC_FLUX_IN(1).units,'Interpreter','Latex');
		export_fig('-png',[obj.paths.plots.diag,'POC_FLUX_IN_diff',num2str(d)],'-m5');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OTHER DIAGNOSTICS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps
% P18
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
	title(['Location of depth slices'],'Interpreter','Latex');
	export_fig('-png',[obj.paths.plots.diag,'trans_locations'],'-m5');
	% Clear data
	obj = clearROMS(obj);
	clearvars -except obj plots pltcnt file
end

