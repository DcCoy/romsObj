function obj = romsDiag(obj,file,plotchoice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to automatically generate diagnostic plots for individual ROMS simulations 
% Expects 'avg' output from ROMS simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
% - obj = romsDiag(obj,file,plotchoice)
%
% Inputs:
% - obj       = ROMS object via romsObj
% - file      = File(s) to diagnose 
% - pltchoice = Diagnostic to produce, see below 
%
%    ------------------------------
%    ---- STANDARD DIAGNOSTICS ----
%    ------------------------------
%
%    1  == Shallow z-slice comparisons
%    2  == transect comparisons
%    3  == 2D field comparisons
%    4  == OMZ thickness comparisons (0, 5, 10, 20, 50uM)
%    5  == POC flux comparisons
%    6  == Gridded N cycle tracers obs
%    7  == Deep z-slice comparisons 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace and opened figures
close all

% Load options
romsOpt;

% Clear objects
obj = clearROMS(obj);

% Initialize plotchoice array
tmpchoice = zeros(1,100);
tmpchoice(plotchoice) = 1;
plotchoice = tmpchoice;

% Process plotchoice(s)
for i = 1:7
    plots(i).on = plotchoice(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load diagnostic options
warning off;
run([diagPath,'/diag_options.m']);

% Change to simulation diagDir and load overrides options
diagDir  = [diagPath,obj.info.simName];
run([diagDir,'/diag_overrides.m']);

% Make plot directory
mkdir([diagDir,'/plots']);
plotsDir = [diagDir,'/plots/',obj.info.runName,'/']; 
mkdir(plotsDir);

% Start plot counter
pltcnt = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D diagnostics (shallow zslices)
% P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Variable loop
    for v = 1:length(vars)
        if ~opt(v)
            disp(['    ...skipping ',vars{v},'...']);
            continue
        end
        % Get data
        obj = clearROMS(obj);
        try
            obj = zslice(obj,vars(v),zdeps,file);
        catch
            disp(['    ...skipping ',vars{v},'...']);
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
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zshallow_',num2str(z),'_roms'],figsQuality);
                close(figs(1));
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                title([obj.diag.(vars{v})(d).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zshallow_',num2str(z),'_diag',num2str(d)],figsQuality);
                close(figs(2));   
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zshallow_',num2str(z),'_diff',num2str(d)],figsQuality);   
            end
        end
    end
    % Clear data
    obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transects
% P2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Variable loop
    for v = 1:length(vars);
        if ~opt(v) 
            disp(['    ...skipping ',vars{v},'...']);
            continue
        end
        % Transect loop
        for l = 1:length(degs)
            % Get ROMS data
            obj = clearROMS(obj);
            %try
                obj = sliceROMS(obj,vars(v),choice{l},degs(l),file,...
                   'zdep',obj.grid.z_avg_dep(obj.grid.z_avg_dep<=zlims(end)));
            %catch
            %    disp(['    ...skipping ',vars{v},'...']);
            %    continue
            %end
            % Get diag data
            obj = sliceDiag(obj,vars(v),choice{l},degs(l),'zlim',zlims(end));
            for d = 1:length(obj.diag.(vars{v}))
                % Extract data
                romsdat = obj.data.avg.(vars{v}).slice;
                diagdat = obj.diag.(vars{v})(d).slice;
                % Get levels
                levs  = linspace(lims{v}(1),lims{v}(2),40);
                dlevs = linspace(dlims{v}(1),dlims{v}(2),42);
                % Make figures
                [figs,cbs] = sliceCmp(obj,romsdat,diagdat,0,'cmap',cmaps{v},'xlims',xlims(l,:),'zlims',zlims(l,:),...
                   'figdim',0.5,'levels',levs,'difflevels',dlevs);
                % ROMS figure
                set(0,'CurrentFigure',figs(1));
                if strcmp(choice{l},'lat')
                    if degs(l)<0
                        title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(abs(degs(l))),'$^oS$'],'Interpreter','Latex');
                    else
                        title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(degs(l)),'$^oN$'],'Interpreter','Latex');
                    end
                elseif strcmp(choice{l},'lon')
                    if degs(l)>180
                        title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(abs(degs(l)-360)),'$^oW$'],'Interpreter','Latex');
                    else
                        title(['ROMS ',obj.data.avg.(vars{v}).name,': ',num2str(degs(l)),'$^oE$'],'Interpreter','Latex');
                    end
                end
                ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_transect_',num2str(l),'_roms'],figsQuality);
                close(figs(1))
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                if strcmp(choice{l},'lat')
                    if degs(l)<0
                        title([obj.diag.(vars{v})(d).name,': ',num2str(abs(degs(l))),'$^oS$'],'Interpreter','Latex');
                    else
                        title([obj.diag.(vars{v})(d).name,': ',num2str(degs(l)),'$^oN$'],'Interpreter','Latex');
                    end
                elseif strcmp(choice{l},'lon')
                    if degs(l)>180
                        title([obj.diag.(vars{v})(d).name,': ',num2str(abs(degs(l)-360)),'$^oW$'],'Interpreter','Latex');
                    else
                        title([obj.diag.(vars{v})(d).name,': ',num2str(degs(l)),'$^oE$'],'Interpreter','Latex');
                    end
                end
                ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_transect_',num2str(l),'_diag',num2str(d)],figsQuality);
                close(figs(2))
                % Difference figure
                set(0,'CurrentFigure',figs(3));
                if strcmp(choice{l},'lat')
                    if degs(l)<0
                        title(['Difference: ',num2str(abs(degs(l))),'$^oS$'],'Interpreter','Latex');
                    else
                        title(['Difference: ',num2str(degs(l)),'$^oN$'],'Interpreter','Latex');
                    end
                elseif strcmp(choice{l},'lon')
                    if degs(l)>180
                        title(['Difference: ',num2str(abs(degs(l)-360)),'$^oW$'],'Interpreter','Latex');
                    else
                        title(['Difference: ',num2str(degs(l)),'$^oE$'],'Interpreter','Latex');
                    end
                end
                ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_transect_',num2str(l),'_diff',num2str(d)],figsQuality);
                close(figs(3))
            end
        end
    end

    % Plot transect map
    close all
    [fig] = quickMap(obj);
    hold on

    % Find and plot lon choices
    idx = find(ismember('lon',choice)==1);
    if ~isempty(idx)
        lons = degs(idx);
        for i = 1:length(lons)
            lony{i} = [-90:0.1:90];
            lonx{i} = [lons(i)*ones(size(lony{i}))];
            [in,~] = inpolygon(lonx{i},lony{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
            m_plot(lonx{i}(in==1),lony{i}(in==1),'--k');
        end
    end

    % Find and plot lat choices
    idx = find(ismember('lat',choice)==1);
    if ~isempty(idx)
        lats = degs(idx);
        for i = 1:length(lats)
            latx{i} = [0:0.1:360];
            laty{i} = [lats(i)*ones(size(latx{i}))];
            [in,~] = inpolygon(latx{i},laty{i},obj.grid.polygon(:,1),obj.grid.polygon(:,2));
            m_plot(latx{i}(in==1),laty{i}(in==1),'--k');
        end
    end
    title(['Location of Transects'],'Interpreter','Latex');
    export_fig(figsFormat,[plotsDir,'trans_locations'],figsQuality);
    % Clear data
    obj = clearROMS(obj);
   % Clear data
   obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface diagnostics
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
   unpackStruct(plots(pltcnt));
    % Variable loop
    for v = 1:length(vars);
        if ~opt(v)
            disp(['    ...skipping ',vars{v},'...']);
            continue
        end
        % Get data
        obj = clearROMS(obj);
        try
            obj = loadData(obj,vars(v),file);
            obj = loadDiag(obj,vars(v),0);
            % Diagnostic loop (sometimes more than 1 available)
            for d = 1:length(obj.diag.(vars{v}))
                close all
                % Get levels
                levs  = linspace(lims{v}(1),lims{v}(2),41);
                dlevs = linspace(dlims{v}(1),dlims{v}(2),43);
                % Extract data
                romsdat    = nanmean(obj.data.avg.(vars{v}).data,3);
                diagdat    = nanmean(obj.diag.(vars{v})(d).slice,3);
                % Make specific corrections
                if strcmp(vars{v},'sustr');
                    romsdat = romsObj.u2rho(romsdat); % convert to rho points
                elseif strcmp(vars{v},'svstr'); % convert to rho points     
                    romsdat = romsObj.v2rho(romsdat);
                elseif strcmp(vars{v},'FG_N2O');
                    romsdat = -romsdat; % loss = positive
                elseif strcmp(vars{v},'SFC_CHL')
                    romsdat(romsdat<0) = 0;
                    romsdat = sqrt(romsdat); % for nonlinear colorbar
                    diagdat = sqrt(diagdat); % for nonlinear colorbar
                    levs    = sqrt(levs);    % for nonlinear levels
                    dlevs   = [dlims{v}(1):0.025:dlims{v}(2)];   % make difference levels simpler
                    dlevs(dlevs>0) = sqrt(dlevs(dlevs>0));       % for nonlinear difference levels
                    dlevs(dlevs<0) = -sqrt(abs(dlevs(dlevs<0))); % for nonlienar difference levels
                end
                % Make figures
                [figs,cbs] = mapCmp(obj,romsdat,diagdat,'cmap',cmaps{v},'levels',levs,'difflevels',dlevs);
                % ROMS figure
                set(0,'CurrentFigure',figs(1));
                title(['ROMS ',obj.data.avg.(vars{v}).name],'Interpreter','Latex');
                ylabel(cbs(1),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                if strcmp(vars{v},'SFC_CHL')
                    cbs(1).XTick      = levs(1:4:end);    % set nonlinear colorbar ticks
                    cbs(1).XTickLabel = levs(1:4:end).^2; % set nonlinear colorbar labels
                end
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_2D_roms'],figsQuality);
                close(figs(1));
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                title([obj.diag.(vars{v})(d).name],'Interpreter','Latex');
                ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                if strcmp(vars{v},'SFC_CHL')
                    cbs(2).XTick      = levs(1:4:end);    % set nonlinear colorbar ticks
                    cbs(2).XTickLabel = levs(1:4:end).^2; % set nonlinear colorbar labels
                end
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_2D_diag',num2str(d)],figsQuality);
                close(figs(2));   
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['Difference'],'Interpreter','Latex');
                ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                if strcmp(vars{v},'SFC_CHL')
                    cbs(3).XTick = dlevs(1:4:end);    % set nonlinear colorbar ticks
                    labels = dlevs(1:4:end).^2;
                    labels(dlevs(1:4:end)<0) = -labels(dlevs(1:4:end)<0);
                    cbs(3).XTickLabel = labels; % set nonlinear colorbar labels
                end
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_2D_diff',num2str(d)],figsQuality);
                close(figs(3));
            end
        catch
             disp(['    ...skipping ',vars{v},'...']);
             continue
        end
    end
    % Clear data
    obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OMZ thickness
% P4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Get OMZ thickness
    if ~opt(1);
        disp(['    ...skipping ',vars{1},'... (no variable)']);
    else
        try
            obj = computeVar(obj,{'OMZ'},file,'thresh',omzthresh);
            obj = loadDiag(obj,{'OMZ'},omzthresh);
            % Make comparison plots
            for d = 1:length(obj.diag.OMZ);
                for th = 1:length(omzthresh)
                    levs  = linspace(lims{th}(1),lims{th}(2),40);
                    dlevs = linspace(dlims{th}(1),dlims{th}(2),42);
                    if ndims(obj.data.avg.OMZ.int)==4
                        romsdat = nanmean(squeeze(obj.data.avg.OMZ.int(:,:,th,:)),3) .* obj.grid.mask_rho;
                    else
                        romsdat = squeeze(obj.data.avg.OMZ.int(:,:,th)) .* obj.grid.mask_rho;
                    end
                    diagdat = squeeze(obj.diag.OMZ(d).slice(:,:,th));
                    [figs,cbs] = mapCmp(obj,romsdat,diagdat,'levels',levs,'difflevels',dlevs,'cmap',cmaps{1});
                    % ROMS figure
                    set(0,'CurrentFigure',figs(1));
                    title(['ROMS: ',obj.data.avg.OMZ.name,'(O$_2$ $<$ ',num2str(omzthresh(th)),' mmol $m^{-3}$)'],'Interpreter','Latex');
                    ylabel(cbs(1),obj.data.avg.OMZ.units,'Interpreter','Latex')
                    export_fig(figsFormat,[plotsDir,'OMZ_roms_th',num2str(omzthresh(th))],figsQuality);
                    % Diag figure
                    set(0,'CurrentFigure',figs(2));
                    title([obj.diag.OMZ(d).name],'Interpreter','Latex');
                    ylabel(cbs(2),obj.data.avg.OMZ.units,'Interpreter','Latex')
                    export_fig(figsFormat,[plotsDir,'OMZ_diag_',num2str(d),'_th',num2str(omzthresh(th))],figsQuality);
                    % Diff figure
                    set(0,'CurrentFigure',figs(3));
                    title(['Difference'],'Interpreter','Latex');
                    ylabel(cbs(3),obj.data.avg.OMZ.units,'Interpreter','Latex')
                    export_fig(figsFormat,[plotsDir,'OMZ_diff_',num2str(d),'_th',num2str(omzthresh(th))],figsQuality);
                    close all
                end
            end
        catch
             disp(['    ...skipping OMZ thickness...']);
        end    
    end
    % Clear data
    obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slice maps of POC_FLUX_IN
% P5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Get POC FLUX IN 
    if ~opt(1);
        disp(['    ...skipping POC_FLUX_IN...']);
    else
        try
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
                export_fig(figsFormat,[plotsDir,name{1},'_roms'],figsQuality);
                close(figs(1));
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                title([obj.diag.POC_FLUX_IN(d).name,': Euphotic'],'Interpreter','Latex');
                ylabel(cbs(2),obj.data.avg.POC_FLUX_IN.units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,name{1},'_diag',num2str(d)],figsQuality);
                close(figs(2));
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['Difference'],'Interpreter','Latex');
                ylabel(cbs(3),obj.data.avg.POC_FLUX_IN.units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,name{1},'_diff',num2str(d)],figsQuality);
            end
        catch
             disp(['    ...skipping POC_FLUX_IN...']);
        end
    end
    % Clear data
    obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gridded obs vs ROMS (including ratios)
% P6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Get filename/path to gridded product
    fname = [valiPath,'/ncycle/Global_gridded_3d_2x2.mat'];
    % Call extract_3d_obs
    extract_3d_obs(obj,fname,plotsDir,rootPath);
    % Call extract_3d_roms
    extract_3d_roms(obj,file,fname,plotsDir,rootPath);
    % Clear data
    obj = clearROMS(obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D deep diagnostics (zslices)
% P7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pltcnt = pltcnt + 1;
if plots(pltcnt).on;
    unpackStruct(plots(pltcnt));
    % Variable loop
    for v = 1:length(vars)
        if ~opt(v);
            disp(['    ...skipping ',vars{v},'...']);
            continue
        end
        % Get data
        obj = clearROMS(obj);
        try    
            obj = zslice(obj,vars(v),zdeps,file);
        catch
            disp(['    ...skipping ',vars{v},'...']);
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
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zdeep_',num2str(z),'_roms'],figsQuality);
                close(figs(1));
                % Diag figure
                set(0,'CurrentFigure',figs(2));
                title([obj.diag.(vars{v})(d).name,': ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                ylabel(cbs(2),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zdeep_',num2str(z),'_diag',num2str(d)],figsQuality);
                close(figs(2));   
                % Diff figure
                set(0,'CurrentFigure',figs(3));
                title(['Difference: ',num2str(zdeps(z)),'m'],'Interpreter','Latex');
                ylabel(cbs(3),obj.data.avg.(vars{v}).units,'Interpreter','Latex');
                export_fig(figsFormat,[plotsDir,sprintf('%02d',v),'_',name{v},'_zdeep_',num2str(z),'_diff',num2str(d)],figsQuality);   
            end
        end
    end
    % Clear data
    obj = clearROMS(obj);
end

% Turn back on warnings
warning on
