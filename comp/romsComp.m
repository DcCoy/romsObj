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
%            ... comparisons will be obj(1) - obj(2)
% - file  = cell array of files to load  
% -         ... file{1} = obj(1) files etc.
% - plotchoice (see below) 
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
%    5  == surface chlA map (log scale)
%    6  == OMZ thickness     (0, 5, 10, 20, 50uM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all

% Load options
romsOpt;
addpath /data/project1/demccoy/ROMS/validation/ncycle

% Clear objects
for i = 1:length(obj)
    obj(i) = clearROMS(obj(i));
end

% Process plotchoice(s)
for i = 1:6
    plots(i).on = plotchoice(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load diagnostic options
warning off;
run([compPath,'/comp_options.m']);

% Change to simulation diagDir and load overrides options
compDir  = [compPath,obj(1).info.simName];
run([compDir,'/comp_overrides.m']);

% Make plots directory
mkdir([compDir,'/plots']);
plotsDir = [compDir,'/plots/',obj(1).info.runName,'_vs_',obj(2).info.runName,'/'];
mkdir(plotsDir);

% Start plot counter
pltcnt = 0;

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
		% Depth loop
        for z = 1:length(zdeps)
            % Get levels
            levs  = linspace(lims{v,z}(1),lims{v,z}(2),40);
            dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
			% Extract data
            roms1dat   = nanmean(squeeze(tmp{1}(:,:,z,:)),3);
            roms2dat   = nanmean(squeeze(tmp{2}(:,:,z,:)),3);
            [figs,cbs] = mapCmp(obj(1),roms1dat,roms2dat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
            % ROMS1 figure
            set(0,'CurrentFigure',figs(1));
            title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_z',num2str(z),'_roms1'],figsQuality);
            close(figs(1));
            % ROMS2 figure
            set(0,'CurrentFigure',figs(2));
            title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_z',num2str(z),'_roms2'],figsQuality);
            close(figs(2));
            % Diff figure
            set(0,'CurrentFigure',figs(3));
            title(['ROMS Difference (1-2): ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
            ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_z',num2str(z),'_roms_diff'],figsQuality);
            close(figs(3));
        end
    end
    % Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
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
			if lats(l) > 0
				title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			else
				title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name,': ',num2str(abs(lats(l))),'$^oS$'],'Interpreter','Latex');
			end
            ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lat',num2str(l),'_roms1'],figsQuality);
            close(figs(1));
            % ROM2 figure
            set(0,'CurrentFigure',figs(2));
			if lats(l)>0
				title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name,': ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			else
				title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name,': ',num2str(abs(lats(l))),'$^oS$'],'Interpreter','Latex');
			end
            ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lat',num2str(l),'_roms2'],figsQuality);
            close(figs(2));
            % Difference figure
            set(0,'CurrentFigure',figs(3));
			if lats(l)>0
				title(['ROMS Difference (1-2): ',num2str(lats(l)),'$^oN$'],'Interpreter','Latex');
			else
				title(['ROMS Difference (1-2): ',num2str(abs(lats(l))),'$^oS$'],'Interpreter','Latex');
			end
            ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lat',num2str(l),'_roms_diff'],figsQuality);
            close(figs(3))
        end
    end
    % Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
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
			if lons(l)>180
				title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name,': ',num2str(abs(lons(l)-360)),'$^oW$'],'Interpreter','Latex');
			else
				title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name,': ',num2str(lons(l)),'$^oE$'],'Interpreter','Latex');
			end
            ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lon',num2str(l),'_roms1'],figsQuality);
            close(figs(1));
            % ROM2 figure
            set(0,'CurrentFigure',figs(2));
			if lons(l)>180
				title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name,': ',num2str(abs(lons(l)-360)),'$^oW$'],'Interpreter','Latex');
			else
				title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name,': ',num2str(lons(l)),'$^oE$'],'Interpreter','Latex');
			end
            ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lon',num2str(l),'_roms2'],figsQuality);
            close(figs(2));
            % Difference figure
            set(0,'CurrentFigure',figs(3));
			if lons(l)>180
				title(['ROMS Difference (1-2): ',num2str(abs(lons(l)-360)),'$^oW$'],'Interpreter','Latex');
			else
				title(['ROMS Difference (1-2): ',num2str(lons(l)),'$^oE$'],'Interpreter','Latex');
			end
            ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
            export_fig(figsFormat,[plotsDir,vars{v},'_lon',num2str(l),'_roms_diff'],figsQuality);
            close(figs(3))
        end
    end
    % Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
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
        title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
        ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
        export_fig(figsFormat,[plotsDir,vars{v},'_roms1'],figsQuality);
        close(figs(1));
        set(0,'CurrentFigure',figs(2));
        title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
        ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
        export_fig(figsFormat,[plotsDir,vars{v},'_roms2'],figsQuality);
        close(figs(2));
        % Diff figure
        set(0,'CurrentFigure',figs(3));
        title(['ROMS Difference (1-2)'],'Interpreter','Latex');
        ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
        export_fig(figsFormat,[plotsDir,vars{v},'_roms_diff'],figsQuality);
        close(figs(3));
    end
    % Clear data
    for i = 1:length(obj)
        obj(i) = clearROMS(obj(i));
    end
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
        title(['ROMS(1) ',obj(1).data.avg.(vars{v}).name],'Interpreter','Latex');
        ylabel(cbs(1),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
        cbs(1).XTick      = absLevs([1:10:101]);
        cbs(1).XTickLabel = absLbls([1:10:101]); 
        export_fig(figsFormat,[plotsDir,vars{v},'_roms1'],figsQuality);
        close(figs(1));
        set(0,'CurrentFigure',figs(2));
        title(['ROMS(2) ',obj(2).data.avg.(vars{v}).name],'Interpreter','Latex');
        ylabel(cbs(2),obj(2).data.avg.(vars{v}).units,'Interpreter','Latex');
        cbs(2).XTick      = absLevs([1:10:101]);
        cbs(2).XTickLabel = absLbls([1:10:101]); 
        export_fig(figsFormat,[plotsDir,vars{v},'_roms2'],figsQuality);
        close(figs(2));
        % Diff figure
        set(0,'CurrentFigure',figs(3));
        title(['ROMS Difference (1-2)'],'Interpreter','Latex');
        ylabel(cbs(3),obj(1).data.avg.(vars{v}).units,'Interpreter','Latex');
        export_fig(figsFormat,[plotsDir,vars{v},'_roms_diff'],figsQuality);
        close(figs(3));
        % Clear data
        for i = 1:length(obj)
            obj(i) = clearROMS(obj(i));
        end
    catch
         disp(['...skipping ',vars{1},'...']);
    end    
end
