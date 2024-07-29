%------------------------------------------------------------------------------------------------
classdef romsObj
%-------------------
% Matlab class used to load ROMS data and make various plots,
% comparison figures against validation data, and many
% other actions.
%
% Requires NCO tools (ncread, ncinfo, ncdump, etc)
%
% Documentation exists at:
% https://github.com/DcCoy/romsObj
%
% Basic instructions:
% - simName = 'pacmed_0p25';
% - runName = 'test';
% - obj = romsObj;                       % initialize object
% - obj = initROMS(obj,simName,runName); % initialize for simName using runName
%
% To view available routines:
% - methods(obj)
%
% For help on object or static methods:
% - help obj.(method)
%------------------
    %----------------------------------------------------------------------------------------
    properties
        info                      % struct containing simulation and variable information (via romsInfo)
        grid                      % struct containing grid data and dimensions (via loadGrid)
        slice                     % struct containing slice coordinates (via sliceROMS/sliceDiag)
        profile                   % struct containing profile coordinates (via getProfile)
        budget                    % struct containing BGC tracer budget results (via getBudg)
        data                      % struct containing ROMS output (via loadData, computeVar, sliceROMS)
        diag                      % struct containing validation data for comparions (via loadDiag, sliceDiag)
        paths = getDiagPaths;     % struct containing paths to data, directories, diagnostic data, etc (via initROMS)
        constants = getConstants; % struct containing useful constants and conversions
    end % end properties
    %----------------------------------------------------------------------------------------
    methods
        %--------------------------------------------------------------------------------
        function obj = initROMS(obj,simName,runName,varargin)
            % -------------------
            % Initialization method: gathers paths and coordinate variables 
            %
            % Usage: 
            % - obj = initROMS(obj,simName,runName,varargin)
            % 
            % Inputs:
            % - simName: ROMS simulation (e.g., peru_chile_0p1, pacmed_0p25)
            % - runName:  Simulation-specific name (matches output directory)
            %
            % Optional Inputs (varargin):
            % - domain:   [xi_start xi_end eta_start eta_end] to analyze a subdomain 
            % - coast:    [km] to mask data within #km of coast
            % - depth:    [depth_top depth_bottom] to mask points outside of depth range 
            % - bathy:    [m] to mask data where bathymetry is less than #m
            %
            % Example:
            % - obj = initROMS(romsObj,'peru_chile_0p1','VKV4_tune2')                         <-- if obj undefined
            % - obj = initROMS(obj,'peru_chile_0p1','VKV4_tune2')                             <-- if obj defined
            % - obj = initROMS(obj,'peru_chile_0p1','VKV4_tune2','domain',[100 200 100 200]); <-- define subdomain
            % -------------------

            %  Begin
            disp(' ');
            disp('---------------------------------');
            disp('---------------------------------');
            disp('            romsObj              ');
            disp('---------------------------------');
            disp('---------------------------------');
            disp(' ');                   
            disp('initROMS: Initializing simulation paths');

            % Load options
            romsOpt
            addpath(scriptsPath)

            % Process inputs 
            A.domain = [];
            A.coast  = [];
            A.depth  = [];
            A.bathy  = [];
            A        = romsObj.parse_pv_pairs(A,varargin);

            % Clear any loaded ROMS data
            if ~isempty(obj.data) | ~isempty(obj.diag)
                obj = clearROMS(obj);
            end

            % Get simPath and grid file
            obj.paths.simPath   = [simPath,simName,'/'];
            tmp = dir([obj.paths.simPath,'grid/']);
            obj.paths.grid      = [obj.paths.simPath,'grid/',tmp(3).name];

            % Get runPath and paths to output
            obj.paths.runPath    = [obj.paths.simPath,runName,'/']; % run directory
            obj.paths.phy_avg    = [obj.paths.runPath,'phy/avg/'];  % physical tracer averages
            obj.paths.phy_his    = [obj.paths.runPath,'phy/his/'];  % physical tracer snapshots
            obj.paths.bgc_avg    = [obj.paths.runPath,'bgc/avg/'];  % biogeochemical tracer averages
            obj.paths.bgc_his    = [obj.paths.runPath,'bgc/his/'];  % biogoechemical tracer snapshots
            obj.paths.dia_avg    = [obj.paths.runPath,'dia/avg/'];  % diagnostic averages
            obj.paths.dia_his    = [obj.paths.runPath,'dia/his/'];  % diagnostics snapshots
            obj.paths.flux_avg   = [obj.paths.runPath,'flux/avg/']; % advective/diffusive tracer flux averages
            obj.paths.flux_his   = [obj.paths.runPath,'flux/his/']; % advective/diffusive tracer flux snapshots
            obj.paths.flx_avg    = [obj.paths.runPath,'flx/avg/'];  % interface flux averages
            obj.paths.flx_his    = [obj.paths.runPath,'flx/his/'];  % interface flux snapshots
            obj.paths.plots.figs = [obj.paths.runPath,'Figures/'];  % run Figures

            % initiate directories if they dont exist
            warning off
            mkdir([obj.paths.runPath,'phy']);
            mkdir([obj.paths.runPath,'bgc']);
            mkdir([obj.paths.runPath,'dia']);
            mkdir([obj.paths.runPath,'phy/avg/']);
            mkdir([obj.paths.runPath,'phy/his/']);
            mkdir([obj.paths.runPath,'bgc/avg/']);
            mkdir([obj.paths.runPath,'bgc/his/']);
            mkdir([obj.paths.runPath,'dia/avg/']);
            mkdir([obj.paths.runPath,'dia/his/']);
            mkdir([obj.paths.runPath,'flux/avg/']);
            mkdir([obj.paths.runPath,'flux/his/']);
            mkdir([obj.paths.runPath,'flx/avg/']);
            mkdir([obj.paths.runPath,'flx/his/']);
            mkdir([obj.paths.plots.figs]);
            warning on

            % Initialize info for subregion
            obj.grid.region.lon_lim   = [];
            obj.grid.region.lat_lim   = [];
            obj.grid.region.coast_lim = [];
            obj.grid.region.dep_lim   = [-inf inf];
            obj.grid.region.bathy_lim = [];
            if ~isempty(A.domain)
                obj.grid.region.lon_lim = [A.domain(1:2)];
                obj.grid.region.lat_lim = [A.domain(3:4)];
            end
            if ~isempty(A.coast)
                obj.grid.region.coast_lim = [A.coast];
            end
            if ~isempty(A.depth)
                obj.grid.region.dep_lim = [A.depth];
            end
            if ~isempty(A.bathy)
                obj.grid.region.bathy_lim = [A.bathy];
            end

            % Initialize info for ROMS variables
            obj.info = A;
            obj.info.simName = simName;
            obj.info.runName = runName;

            % Get info for files
            obj = romsInfo(obj);

            % Load grid
            obj = loadGrid(obj);
        end % end method initROMS

        %--------------------------------------------------------------------------------
        function [obj] = romsInfo(obj)
            % ----------------------
            % Obtains info from sample roms files in run directory
            % Automatically called during initROMS
            % Assumes all file output is identical between files
            %
            % Usage:
            % - [obj] = romsInfo(obj) 
            % ----------------------
            disp('romsInfo: Loading file meta data (ncinfo) into obj.info');

            % List file_types
            file_types  = {'phy','bgc','dia','flux','flx'};
            file_dirs   = {'avg','his'};

            % Get info for each file type
            good_dir = 0;
            for i = 1:length(file_types)
                for j = 1:length(file_dirs)
                    % Get files and simplify paths
                    tmp = dir(obj.paths.([file_types{i},'_',file_dirs{j}]));
                    tmpfiles = {tmp.name};
                    clear tmp

                    % Remove first 2 entries ({.},{..})
                    tmpfiles(1:2) = [];

                    % Skip if no files exist
                    if isempty(tmpfiles)
                        disp(['    NOTE(romsInfo): No ',[file_types{i},'_',file_dirs{j}],' files found...skipping']);
                        obj.info.([file_types{i},'_',file_dirs{j}]) = [];
                        continue
                    else
                        good_dir = 1;
                    end

                    % Get sample file (for ncinfo)
                    file_path = obj.paths.([file_types{i},'_',file_dirs{j}]);

                    % Grab ncinfo for each file
                    for k = 1:length(tmpfiles) 

                        % Call ncinfo (once, assumes all files are identical)
                        if k == 1
                            tmp = ncinfo([file_path,tmpfiles{k}]);
                            % Remove useless fields from obj.info
                            tmp = rmfield(tmp,'Name');  
                            tmp = rmfield(tmp,'Groups');
                            tmp = rmfield(tmp,'Format');
                            % Extract dimensions for easier access
                            dims = {tmp.Dimensions.Name};
                            for l = 1:length(dims)
                                tmp.(dims{l}) = tmp.Dimensions(l).Length;
                            end
                            % Apply any regional overrides
                            if ~isempty(obj.grid.region.lon_lim); 
                                tmp.xi_rho = [obj.grid.region.lon_lim];
                                tmp.xi_u   = [obj.grid.region.lon_lim(1) obj.grid.region.lon_lim(2)-1];
                            else
                                tmp.xi_rho = [1 tmp.xi_rho];
                                try;
                                    tmp.xi_u = [1 tmp.xi_u];
                                catch; 
                                    tmp.xi_u = [1 tmp.xi_rho(end)-1];
                                end
                            end
                            if ~isempty(obj.grid.region.lat_lim); 
                                tmp.eta_rho = [obj.grid.region.lat_lim];
                                tmp.eta_v   = [obj.grid.region.lat_lim(1) obj.grid.region.lat_lim(2)-1];
                            else
                                tmp.eta_rho = [1 tmp.eta_rho];
                                try; 
                                    tmp.eta_v = [1 tmp.eta_v];
                                catch;
                                    tmp.eta_v = [1 tmp.eta_rho(end)-1];
                                end
                            end
                        end
                            
                        % Save into struct
                        obj.info.([file_types{i},'_',file_dirs{j}])(k) = tmp;

                        % Override filename
                        filename = obj.info.([file_types{i},'_',file_dirs{j}])(k).Filename;
                        fidx = strfind(filename,tmpfiles{1});
                        filename(fidx:end) = tmpfiles{k};
                        obj.info.([file_types{i},'_',file_dirs{j}])(k).Filename = filename;
                    end
                end
            end

            % Kill object if no files are found
            if good_dir == 0
                disp('    ERROR(romsInfo): No ROMS files found, check initROMS inputs or fix output directories');
                kill
            end

            % Load other parameter values
            params = {'rho0','theta_s','theta_b','hc','dt','ntimes','nwrt','navg'};
            for i = 1:length(params)
                if ~isempty(obj.info.phy_avg)
                    idx = find(strcmp(params{i},{obj.info.phy_avg(1).Attributes.Name})==1);    
                    if ~isempty(idx)
                        obj.info.params.(params{i}) = obj.info.phy_avg(1).Attributes(idx).Value;
                    end
                end
            end
        end % end method romsInfo

        %--------------------------------------------------------------------------------
        function [obj] = loadGrid(obj)
            % ----------------------
            % Loads 2D grid information into obj.grid
            % Automatically called during initROMS
            %
            % Usage:
            % - obj = loadGrid(obj)
            % ----------------------
            disp('loadGrid: Loading grid data into obj.grid');

            % Requires romsInfo
            try; obj.info;
            catch; obj = romsInfo(obj); 
            end

            % Grab z_avg_deps (WOA-18 depths, useful for diagnostics/validation)
            romsOpt
            tmp = load([valiPath,'wcoord.mat']);
            obj.grid.z_avg_dep = tmp.wcoord0p25.depth;

            % Extract dimensions from file
            file_types = {'phy_avg','phy_his','bgc_avg','bgc_his','dia_avg','dia_his','flx_avg','flx_his'};
            for i = 1:length(file_types)
                if isempty(obj.info.(file_types{i}));
                    continue
                else
                    % Save grid dimensions (do it only once)
                    if ~isfield(obj.grid,'xi_rho') | ~isfield(obj.grid,'s_rho')
                        obj.grid.xi_rho  = obj.info.(file_types{i})(1).xi_rho;
                        obj.grid.eta_rho = obj.info.(file_types{i})(1).eta_rho;
                        obj.grid.xi_u    = obj.info.(file_types{i})(1).xi_u;
                        obj.grid.eta_v   = obj.info.(file_types{i})(1).eta_v;
                        obj.grid.s_rho   = obj.info.(file_types{i})(1).s_rho;
                        % Create regional indices
                        % rho points
                        obj.grid.rho2D = [obj.grid.xi_rho(1) obj.grid.eta_rho(1);...
                            diff(obj.grid.xi_rho)+1 diff(obj.grid.eta_rho)+1];
                        obj.grid.rho3D = [obj.grid.xi_rho(1) obj.grid.eta_rho(1) 1;...
                            diff(obj.grid.xi_rho)+1 diff(obj.grid.eta_rho)+1 inf];
                        obj.grid.rho4D = [obj.grid.xi_rho(1) obj.grid.eta_rho(1) 1 1;...
                            diff(obj.grid.xi_rho)+1 diff(obj.grid.eta_rho)+1 inf inf];
                        % u points
                        obj.grid.u2D = obj.grid.rho2D; obj.grid.u2D(2,1) = obj.grid.u2D(2,1)-1;
                        obj.grid.u3D = obj.grid.rho3D; obj.grid.u3D(2,1) = obj.grid.u3D(2,1)-1;
                        obj.grid.u4D = obj.grid.rho4D; obj.grid.u4D(2,1) = obj.grid.u4D(2,1)-1;
                        % v points
                        obj.grid.v2D = obj.grid.rho2D; obj.grid.v2D(2,2) = obj.grid.v2D(2,2)-1;
                        obj.grid.v3D = obj.grid.rho3D; obj.grid.v3D(2,2) = obj.grid.v3D(2,2)-1;
                        obj.grid.v4D = obj.grid.rho4D; obj.grid.v4D(2,2) = obj.grid.v4D(2,2)-1;
                    end
                end
            end

            % Load basic grid
            fields = {'lon_rho','lat_rho','mask_rho','h','pm','pn','angle'};
            for i = 1:length(fields)
                obj.grid.(fields{i}) = double(ncread(obj.paths.grid,fields{i},...
                                       [obj.grid.rho2D(1,:)],[obj.grid.rho2D(2,:)]));
            end

            % lon360 fix and mask_rho fix
            obj.grid.lon_rho(obj.grid.lon_rho<0) = obj.grid.lon_rho(obj.grid.lon_rho<0)+360;
            obj.grid.mask_rho(obj.grid.mask_rho==0) = NaN;

            % Basic dimensions
            obj.grid.nx = diff(obj.grid.xi_rho)+1;
            obj.grid.ny = diff(obj.grid.eta_rho)+1; 
            obj.grid.nz = obj.grid.s_rho; 
            obj.grid.ndim_xy    = [obj.grid.nx,obj.grid.ny];
            obj.grid.ndim_xyz   = [obj.grid.nx,obj.grid.ny,obj.grid.nz];
            obj.grid.minlon_rho = nanmin(obj.grid.lon_rho(:));
            obj.grid.maxlon_rho = nanmax(obj.grid.lon_rho(:));
            obj.grid.minlat_rho = nanmin(obj.grid.lat_rho(:));
            obj.grid.maxlat_rho = nanmax(obj.grid.lat_rho(:));

            % Get grid polygon
            obj.grid.polygon(:,1) = [obj.grid.lon_rho(1,:) obj.grid.lon_rho(:,end)'...
                                     fliplr(obj.grid.lon_rho(end,:)) flipud(obj.grid.lon_rho(:,1))'];
            obj.grid.polygon(:,2) = [obj.grid.lat_rho(1,:) obj.grid.lat_rho(:,end)'...
                                     fliplr(obj.grid.lat_rho(end,:)) flipud(obj.grid.lat_rho(:,1))'];

            % Apply mask
            % Try loading u/v fields
            try
                obj.grid.lon_u  = double(ncread(obj.paths.grid,'lon_u',...
                                  [obj.grid.u2D(1,:)],[obj.grid.u2D(2,:)]));
                obj.grid.lat_u  = double(ncread(obj.paths.grid,'lat_u',...
                                  [obj.grid.u2D(1,:)],[obj.grid.u2D(2,:)]));
                obj.grid.mask_u = double(ncread(obj.paths.grid,'mask_u',...
                                  [obj.grid.u2D(1,:)],[obj.grid.u2D(2,:)]));
                obj.grid.lon_v  = double(ncread(obj.paths.grid,'lon_v',...
                                  [obj.grid.v2D(1,:)],[obj.grid.v2D(2,:)]));
                obj.grid.lat_v  = double(ncread(obj.paths.grid,'lat_v',...
                                  [obj.grid.v2D(1,:)],[obj.grid.v2D(2,:)]));
                obj.grid.mask_v = double(ncread(obj.paths.grid,'mask_v',...
                                  [obj.grid.v2D(1,:)],[obj.grid.v2D(2,:)]));
            catch
                obj.grid.lon_u  = (obj.grid.lon_rho(1:end-1,:) + obj.grid.lon_rho(2:end,:))./2;
                obj.grid.lat_u  = (obj.grid.lat_rho(1:end-1,:) + obj.grid.lat_rho(2:end,:))./2;
                obj.grid.mask_u = obj.grid.mask_rho(1:end-1,:);
                obj.grid.lon_v  = (obj.grid.lon_rho(:,1:end-1) + obj.grid.lon_rho(:,2:end))./2;
                obj.grid.lat_v  = (obj.grid.lat_rho(:,1:end-1) + obj.grid.lat_rho(:,2:end))./2;
                obj.grid.mask_v = obj.grid.mask_rho(:,1:end-1);
            end

            % Get area dimensions
            obj.grid.area     = (1./(obj.grid.pm .* obj.grid.pn));
            obj.grid.dx       = 1./(obj.grid.pm);
            obj.grid.dy       = 1./(obj.grid.pn);
            obj.grid.area_rho = double(obj.grid.dx.*obj.grid.dy);

            % Fix masks
            obj.grid.mask_u(obj.grid.mask_u==0) = NaN;
            obj.grid.mask_v(obj.grid.mask_v==0) = NaN;

            % Make grid data all single
            obj.grid = romsObj.struct2double(obj.grid);

            % If 'coast' is called during initROMS 
            if obj.grid.region.coast_lim>(-inf)
                % Calculate distance-from-coast (m)
                if isfield(obj.grid,'coastdist') == 0 
                    disp('    Refining mask via Dist2Coast');
                    obj = Dist2Coast(obj);
                end
                % Update mask
                obj.grid.mask_rho(obj.grid.coastdist < obj.grid.region.coast_lim*1000) = NaN;
            end

            % If 'bathy' is called during initROMS
            if ~isempty(obj.grid.region.bathy_lim)
                % Update mask to exclude shallow points
                disp('    Refining mask to exclude shallow points');
                obj.grid.mask_rho(obj.grid.h <= obj.grid.region.bathy_lim) = NaN;
            end
        end % end method loadGrid

        %--------------------------------------------------------------------------------
        function [obj] = loadDepth(obj,file,varargin) 
            % ----------------------
            % Load Z-grid info (z_r,z_w), which are output-file dependent.
            % Also updates 3D mask
            %
            % Usage:
            % - obj = loadDepth(obj,file,varargin) 
            %
            % Inputs:
            % - file = (#s) load specific files, see obj.info 
            %
            % Optional:
            % - type = file type (his, avg (default), phys_flux)
            % - full = if (1), also loads z_u, z_v, z_w
            % - mean = if (1), average across all files
            %
            % Example:
            % - obj = loadDepth(obj,1,'type','avg');
            % ----------------------
            disp('loadDepth: Grabbing z_r, z_w, Hz');

            % Call loadGrid
            try; obj.grid.h;
            catch; obj = loadGrid(obj);
            end
            
            % Process optional inputs
            A.type = 'avg'; % file type
            A.full = 0;     % full switch
            A.mean = 0;     % mean switch
            A      = romsObj.parse_pv_pairs(A,varargin);

            % Check for bad inputs
            if isempty(obj.info.(['phy_',A.type]))
                disp(['    ERROR(loadDepth): No physical output found in phy_',A.type,'! Check type or runPath directories']);
                kill
            end

            % Set file type
            file_type = ['phy_',A.type];

            % Grab files
            files = {obj.info.(file_type).Filename};

            % Initialize matrices to fill
            files = files(file);
            nt    = sum(cell2mat({obj.info.(file_type)(file).time}));
            obj.grid.(A.type).Hz  = nan(obj.grid.nx,obj.grid.ny,obj.grid.nz,nt);
            obj.grid.(A.type).z_r = nan(obj.grid.nx,obj.grid.ny,obj.grid.nz,nt);
            if A.full
                obj.grid.(A.type).z_u = nan(obj.grid.nx-1,obj.grid.ny,obj.grid.nz,nt);
                obj.grid.(A.type).z_v = nan(obj.grid.nx,obj.grid.ny-1,obj.grid.nz,nt);
                obj.grid.(A.type).z_w = nan(obj.grid.nx,obj.grid.ny,obj.grid.nz-1,nt);
            end
            
            % Go through file(s) and load data (or calculate it via zlevs4)
            for ff = 1:length(files)
                fname = [files{ff}];

                % Load SLA
                tmp.zeta = ncread(fname,'zeta',[obj.grid.rho3D(1,:)],[obj.grid.rho3D(2,:)]);

                % Average for more accurate depth calc
                tmp.h = obj.grid.h;

                % Initialize
                tmp.z_r = nan(size(tmp.zeta,1),size(tmp.zeta,2),obj.grid.s_rho,size(tmp.zeta,3));
                tmp.z_w = nan(size(tmp.zeta,1),size(tmp.zeta,2),obj.grid.s_rho+1,size(tmp.zeta,3));

                % Get depths
                for i = 1:size(tmp.zeta,3);
                    z_r = zlevs4(tmp.h,tmp.zeta(:,:,i),...
                        obj.info.params.theta_s,obj.info.params.theta_b,...
                        obj.info.params.hc,obj.grid.s_rho,'r','new2012');
                    z_w = zlevs4(tmp.h,tmp.zeta(:,:,i),...
                        obj.info.params.theta_s,obj.info.params.theta_b,...
                        obj.info.params.hc,obj.grid.s_rho,'w','new2012');
                    z_r = permute(z_r,[2 3 1]);
                    if A.full
                        z_u = (z_r(1:end-1,:,:) + z_r(2:end,:,:))./2;
                        z_v = (z_r(:,1:end-1,:) + z_r(:,2:end,:))./2;
                    end
                    z_w = permute(z_w,[2 3 1]);
                    tmp.z_r(:,:,:,i) = z_r;
                    if A.full
                        tmp.z_u(:,:,:,i) = z_u;
                        tmp.z_v(:,:,:,i) = z_v;
                    end
                    tmp.z_w(:,:,:,i) = z_w;
                end
                tmp.Hz = diff(tmp.z_w,1,3);

                % Save output
                if A.mean
                    if ff == 1 
                        out.Hz  = [];
                        out.z_r = [];
                        if A.full
                            out.z_u = [];
                            out.z_v = [];
                            out.z_w = [];
                        end
                    end
                    out.Hz  = [out.Hz + (tmp.Hz)];
                    out.z_r = [out.z_r + (tmp.z_r)];
                    if A.full
                        out.z_u = [out.z_u + (tmp.z_u)];
                        out.z_v = [out.z_v + (tmp.z_v)];
                        out.z_w = [out.z_w + (tmp.z_w)];
                    end
                else
                    out.Hz{ff} = tmp.Hz;
                    out.z_r{ff} = tmp.z_r;
                    if A.full
                        out.z_u{ff} = tmp.z_u;
                        out.z_v{ff} = tmp.z_v;
                        out.z_w{ff} = tmp.z_w;
                    end
                end
                clear tmp
            end

            % If the average of depths across multiple files are requested...
            if A.mean
                obj.grid.(A.type).Hz  = (out.Hz ./ ff);
                obj.grid.(A.type).z_r = (out.z_r ./ ff);
                if A.full
                    obj.grid.(A.type).z_u = (out.z_u ./ ff);
                    obj.grid.(A.type).z_v = (out.z_v ./ ff);
                    obj.grid.(A.type).z_w = (out.z_w ./ ff);
                end
            % Else, concatenate
            else
                obj.grid.(A.type).Hz  = cat(4,out.Hz{:});
                obj.grid.(A.type).z_r = cat(4,out.z_r{:});
                if A.full
                    obj.grid.(A.type).z_u = cat(4,out.z_u{:});
                    obj.grid.(A.type).z_v = cat(4,out.z_v{:});
                    obj.grid.(A.type).z_w = cat(4,out.z_w{:});
                end
            end

            % Make 3D area_rho, grid_mask with depth limits applied
            if A.mean
                obj.grid.(A.type).mask_rho3d = repmat(obj.grid.mask_rho,[1 1 obj.grid.nz]);
                if A.full
                    try
                        obj.grid.(A.type).mask_u3d = repmat(obj.grid.mask_u,[1 1 obj.grid.nz]);
                        obj.grid.(A.type).mask_v3d = repmat(obj.grid.mask_v,[1 1 obj.grid.nz]);
                    catch
                    end
                    obj.grid.(A.type).area3d = repmat(obj.grid.area_rho,[1 1 obj.grid.nz]);
                    obj.grid.(A.type).volume = obj.grid.(A.type).area3d .* obj.grid.(A.type).Hz;                
                end
            else
                obj.grid.(A.type).mask_rho3d = repmat(obj.grid.mask_rho,[1 1 obj.grid.nz nt]);
                if A.full
                    try
                        obj.grid.(A.type).mask_u3d = repmat(obj.grid.mask_u,[1 1 obj.grid.nz nt]);
                        obj.grid.(A.type).mask_v3d = repmat(obj.grid.mask_v,[1 1 obj.grid.nz nt]);
                    catch
                    end
                    obj.grid.(A.type).area3d = repmat(obj.grid.area_rho,[1 1 obj.grid.nz nt]);
                    obj.grid.(A.type).volume = obj.grid.(A.type).area3d .* obj.grid.(A.type).Hz;
                end
            end

            % Apply depth limits to 3D mask (if 
            if ~isempty(obj.grid.region.dep_lim);            
                obj.grid.(A.type).mask_rho3d(obj.grid.(A.type).z_r<obj.grid.region.dep_lim(1)) = NaN;
                obj.grid.(A.type).mask_rho3d(obj.grid.(A.type).z_r>obj.grid.region.dep_lim(2)) = NaN;
                if A.full
                    try
                        obj.grid.(A.type).mask_u3d(obj.grid.(A.type).z_u<obj.grid.region.dep_lim(1)) = NaN;
                        obj.grid.(A.type).mask_u3d(obj.grid.(A.type).z_u>obj.grid.region.dep_lim(2)) = NaN;
                        obj.grid.(A.type).mask_v3d(obj.grid.(A.type).z_v<obj.grid.region.dep_lim(1)) = NaN;
                        obj.grid.(A.type).mask_v3d(obj.grid.(A.type).z_v>obj.grid.region.dep_lim(2)) = NaN;
                    catch
                    end
                    obj.grid.(A.type).area3d(obj.grid.(A.type).z_r<obj.grid.region.dep_lim(1)) = NaN;
                    obj.grid.(A.type).area3d(obj.grid.(A.type).z_r>obj.grid.region.dep_lim(2)) = NaN;
                end
            end
        end % end method loadDepth

        %--------------------------------------------------------------------------------
        function obj = loadData(obj,vars,file,varargin)
            % -------------------
            % Method to Load ROMS data.
            % Use dispVars or obj.info to list variables
            %
            % Usage:
            % - obj = loadData(obj,vars,file,varargin)
            %
            % Inputs:
            % - vars   = variable(s) to load, as a cell array
            % - file   = (#s) files to load (time dimension will be concatenated)
            %
            % Optional inputs (varargin):
            % - type    = 'avg','his','rst','bgc','phy','phys_flux' 
            % - mean    = (1) to average across multiple files
            % 
            % Examples:
            % - obj = loadData(obj)
            % - obj = loadData(obj,{'temp','salt'},1,'type','avg');
            % -------------------

            % Defaults for optional arguments
            romsOpt
            A.type   = ['avg']; % file type
            A.mean   = [0];     % averaging switch
            A        = romsObj.parse_pv_pairs(A,varargin);

            % Restrict file types to 'avg' or 'his'
            file_types = {'phy_avg','phy_his',...
                          'bgc_avg','bgc_his',...
                          'dia_avg','dia_his',...
                          'flux_avg','flux_his',...
                          'flx_avg','flx_his'};
            idx = contains(file_types,A.type);
            file_types(idx==0) = [];
    
            % Ignore empty file types
            idx = ones(size(file_types));;
            for i = 1:length(file_types)
                if isempty(obj.info.(file_types{i}));
                    idx(i) = 0;
                end
            end
            file_types(idx==0) = [];

            % Check inputs
            if ~isempty(vars) % vars
                for i = 1:length(vars)
                    if ~iscell(vars(i));
                        disp('    ERROR(loadData): vars must be a cell array');
                        return
                    end
                end
            end

            % Load each variable into structure
            for i = 1:length(vars)
                disp(['loadData: Loading ',vars{i},' into struct']);
                % Find location
                file_dir = [];
                for j = 1:length(file_types)
                    all_vars = {obj.info.(file_types{j})(file(1)).Variables.Name};
                    if ismember(vars{i},all_vars);
                        file_dir = j;
                        break
                    end
                end
                if isempty(file_dir)
                    % Call computeVar to compute variables (if they don't exist in output)
                    try; obj = computeVar(obj,vars(i),file,'type',A.type,'mean',A.mean); continue
                    catch; disp(['    ERROR(loadData): ',vars{i},' is not an variable, and isnt available via computeVar']);
                    end
                end

                % Check for dimensions, grab index
                idx      = find(strcmp(vars{i},{obj.info.(file_types{file_dir})(file(1)).Variables.Name})==1);
                tmp_size = squeeze(obj.info.(file_types{file_dir})(file(1)).Variables(idx).Size); 
                % 4D data (x,y,z,t)
                if length(tmp_size)==4
                    ind = obj.grid.rho4D;    
                    if ismember(vars{i},{'u','sustr','sustr_avg'});
                        ind = obj.grid.u4D;
                    elseif ismember(vars{i},{'v','svstr','svstr_avg'});
                        ind = obj.grid.v4D;
                    end
                % 3D data (x,y,z) or (x,y,t)
                elseif length(tmp_size)==3
                    ind = obj.grid.rho3D;
                    if ismember(vars{i},{'u','sustr','sustr_avg'});
                        ind = obj.grid.u3D;
                    elseif ismember(vars{i},{'v','svstr','svstr_avg'});
                        ind = obj.grid.v3D;
                    end
                % 2D data (x,y)
                elseif length(tmp_size)==2
                    ind = obj.grid.rho2D;
                    if ismember(vars{i},{'u','sustr','sustr_avg'});
                        ind = obj.grid.u2D;
                    elseif ismember(vars{i},{'v','svstr','svstr_avg'});
                        ind = obj.grid.v2D;
                    end
                else
                    ind = [1 inf];
                end

                % Grab variable info from output file
                atts = {obj.info.(file_types{file_dir})(file(1)).Variables(idx).Attributes.Name};
                idxname  = find(strcmp('long_name',atts)==1);
                idxunits = find(strcmp('units',atts)==1);
                % Try loading output name
                if ~isempty(idxname)
                    obj.data.(A.type).(vars{i}).name  = obj.info.(file_types{file_dir})(file(1)).Variables(idx).Attributes(idxname).Value;
                else
                    obj.data.(A.type).(vars{i}).name = vars{i};
                end
                % Try loading output units
                if ~isempty(idxunits)
                    obj.data.(A.type).(vars{i}).units = obj.info.(file_types{file_dir})(file(1)).Variables(idx).Attributes(idxunits).Value;
                else
                    obj.data.(A.type).(vars{i}).units = 'N/A';
                end

                % Go through all files and load data
                for ff = 1:length(file)

                    % Load data into temporary struct
                    tmp.data = ncread([obj.info.(file_types{file_dir})(file(ff)).Filename],vars{i},[ind(1,:)],[ind(2,:)]);
                    tmp.dims = {obj.info.(file_types{file_dir})(file(1)).Variables(idx).Dimensions.Name}; 

                    % Check for averaging switch
                    if A.mean
                        if ff == 1
                            out.data = [];
                        end
                        if length(ind)>=2
                            % Apply mask
                            if ismember(vars{i},{'u','sustr','sustr_avg'});
                                out.data = [out.data + (tmp.data .* obj.grid.mask_u)];
                            elseif ismember(vars{i},{'v','svstr','svstr_avg'});
                                out.data = [out.data + (tmp.data .* obj.grid.mask_v)];
                            else
                                out.data = [out.data + (tmp.data .* obj.grid.mask_rho)];
                            end
                        else
                            out.data = [out.data + tmp.data];
                        end
                    % Else, set up cell array for concatenation
                    else    
                        if length(ind)>=2
                            % Apply mask
                            if ismember(vars{i},{'u','sustr','sustr_avg'});
                                out.data{ff} = tmp.data .* obj.grid.mask_u;
                            elseif ismember(vars{i},{'v','svstr','svstr_avg'});
                                out.data{ff} = tmp.data .* obj.grid.mask_v;
                            else
                                out.data{ff} = tmp.data .* obj.grid.mask_rho;
                            end
                        else
                            out.data{ff} = tmp.data;
                        end
                    end
                end

                % If average of multiple files is requested...
                if A.mean
                    obj.data.(A.type).(vars{i}).data = (out.data ./ ff);
                % Else, concatenate
                else
                    obj.data.(A.type).(vars{i}).data = cat(length(ind),out.data{:});
                end

                % If loading 'rho', correct for density anomaly (rho0)
                if strcmp(vars{i},'rho');
                    % Remove 'anomaly'
                    disp('    NOTE(loadData): Correcting density anomaly (rho0) --> density');
                    obj.data.(A.type).rho.data = obj.data.(A.type).rho.data + (obj.info.params.rho0 - 1000);
                    obj.data.(A.type).rho.name = 'averaged density';
                end

                % Save dimensions in structure
                obj.data.(A.type).(vars{i}).dims = tmp.dims;
                clear out

                % Replace strings
                for j = 1:length(strings_to_replace)
                    if strcmp(strings_to_replace{j,1},obj.data.(A.type).(vars{i}).units)==1;
                        obj.data.(A.type).(vars{i}).units = strings_to_replace{j,2};
                    end
                end

                % Replace instances of underscores
                if contains(obj.data.(A.type).(vars{i}).name,'_');
                    obj.data.(A.type).(vars{i}).name = replace(obj.data.(A.type).(vars{i}).name,'_',' ');
                end

                % Replace avg with averaged
                if contains(obj.data.(A.type).(vars{i}).name,'avg');
                    obj.data.(A.type).(vars{i}).name = replace(obj.data.(A.type).(vars{i}).name,'avg','averaged');
                end
            end
        end % end method loadData

        %--------------------------------------------------------------------------------
        function obj = computeVar(obj,vars,file,varargin)
            % ------------------
            % Computes additional fields like wind stress, Okubo-Weiss, AOU, etc. 
            % 
            % List of available variables:
            % - pres     (pressure)
            % - rho      (density)
            % - bvf      (Brunt Vaisala frequency)
            % - pv       (potential vorticity)
            % - SSH      (sea-surface height, corrected via zeta) 
            % - ws       (wind-stress)
            % - wsc      (wind-stress curl)
            % - NPP      (net primary production)
            % - nstar    (N* via NO3 - 16*PO4 + 2.9);
            % - MLD      (mixed layer depth)
            % - SFC_CHL  (surface integrated chla, to compare against satellite)
            % - OW       (Okubo-Weiss)
            % - vort     (relative vorticity)
            % - sN       (normal component of strain)
            % - sS       (shear component of strain)
            % - Jden_N2O (SMS of N2O via denitrification)
            % - Jnit_N2O (SMS of N2O via nitrification)
            % - AOU      (Apparent oxygen utilization)
            % - DeltaN2O (Supersaturated N2O)
            % - OMZ      (OMZ thickness)
            % - NH4vNO2  (NH4 vs NO2 ratio)
            % - NO2vNH4  (NO2 vs NH4 ratio)
            % - NOX      (NO2+NO3)
            % - N2O_BRY  (N2O_BOU + N2O_ATM)
            % - NO2AMMOX (AMMOX - 2.*(N2OAMMOX))
            % - sustr    (zonal wind stress, via old naming style 'sustr_avg')
            % - svstr    (meridional wind stress, via old naming style 'svstr_avg')
            %
            % Usage:
            % - obj = computeVar(obj,vars,file,varargin)
            %
            % Inputs:
            % - vars = variable(s) to compute, as a cell array    
            % - file = file number to extract data from
            %
            % Inputs (varargin)
            % - type   = 'avg' (default), 'his'
            % - dep/ip = depths/isopycnal to slice Okubo-Weiss fields on
            % - thresh = thresholds, if needed
            % - mean   = (1) to average across multiple fields
            %
            % Example:
            % - obj = computeVar(obj,{'MLD'},1,'type','avg');
            % ------------------

            % process inputs
            A.type    = 'avg'; % file type
            A.ip      = [];    % isopycnal for Okubo_Weiss calc
            A.dep     = [];    % depth for Okubo_Weiss calc
            A.thresh  = [];    % uM O2 for OMZ thickness calc
            A.mean    = 0;     % averaging switch
            A         = romsObj.parse_pv_pairs(A,varargin);

            % dims for vertical computations
            nx = obj.grid.nx;
            ny = obj.grid.ny;
            nz = obj.grid.nz;

            % process requests
            for i = 1:length(vars)
                % pressure calc
                if strcmp(vars{i},'pres') % & ~isfield(obj.data,'pres');
                    disp('computeVar: Calculating sea pressure');
                    nt = sum(cell2mat({obj.info.(['phy_',A.type])(file).time}));
                    if isempty(obj.grid.(A.type))
                        obj = loadDepth(obj,file);
                    end
                    obj.data.(A.type).pres.data = sw_pres(-obj.grid.(A.type).z_r,repmat(obj.grid.lat_rho,[1 1 nz nt]));
                    obj.data.(A.type).pres.data = obj.data.(A.type).pres.data .* obj.grid.(A.type).mask_rho3d;
                    obj.data.(A.type).pres.name  = 'averaged Pressure';
                    obj.data.(A.type).pres.units = 'dbar';
                    if length(size(obj.grid.(A.type).z_r))==3
                        obj.data.(A.type).pres.dims = {'xi_rho','eta_rho','s_rho'};
                    elseif length(size(obj.grid.(A.type).z_r))==4
                        obj.data.(A.type).pres.dims = {'xi_rho','eta_rho','s_rho','time'};
                    end
                % rho calc
                elseif strcmp(vars{i},'rho');
                    disp('computeVar: Calculating sea water density');
                    obj = loadData(obj,{'temp','salt'},file,'type',A.type,'mean',A.mean);
                    tmp = sw_dens0(obj.data.(A.type).salt.data,obj.data.(A.type).temp.data);
                    if isempty(obj.grid.(A.type))
                        obj = loadDepth(obj,file);
                    end
                    obj.data.(A.type).rho.data = tmp .* obj.grid.(A.type).mask_rho3d;
                    obj.data.(A.type).rho.name = 'averaged density';
                    obj.data.(A.type).rho.units = 'kg m$^{-3}$';
                    if length(size(obj.grid.(A.type).z_r))==3
                        obj.data.(A.type).rho.dims = {'xi_rho','eta_rho','s_rho'};
                    elseif length(size(obj.grid.(A.type).z_r))==4
                        obj.data.(A.type).rho.dims = {'xi_rho','eta_rho','s_rho','time'};
                    end                
                % bvf calc
                elseif strcmp(vars{i},'bvf') | strcmp(vars{i},'pv');
                    disp('computeVar: Calculating bvf (Brunt Vaisala)');
                    nt = sum(cell2mat({obj.info.(['phy_',A.type])(file).time}));
                    obj = computeVar(obj,{'pres'},file,'type',A.type);
                    if isfield(obj.data.(A.type),'temp');
                        tt = 1;
                        origtemp = obj.data.(A.type).temp;
                    else
                        tt = 0;
                    end
                    if isfield(obj.data.(A.type),'salt');
                        ss = 1;
                        origsalt = obj.data.(A.type).salt;
                    else
                        ss = 0;
                    end
                    obj = loadData(obj,{'temp','salt'},file,'type',A.type,'mean',A.mean);
                    tmptemp = reshape(sw_temp(obj.data.(A.type).salt.data(:),obj.data.(A.type).temp.data(:),...
                                              obj.data.(A.type).pres.data(:),0),size(obj.data.(A.type).salt.data)); 
                    tmplat  = repmat(obj.grid.lat_rho,1,1,nz,nt);
                    % Permute to get depth in front
                    tmptemp = permute(obj.data.(A.type).temp.data,[3 1 2 4]);
                    tmpsalt = permute(obj.data.(A.type).salt.data,[3 1 2 4]);
                    tmppres = permute(obj.data.(A.type).pres.data,[3 1 2 4]);
                    tmplat  = permute(tmplat,[3 1 2 4]);
                    % Reshape to MxN
                    dims = size(tmptemp);
                    tmptemp = reshape(tmptemp,dims(1),dims(2)*dims(3)*dims(4)); 
                    tmpsalt = reshape(tmpsalt,dims(1),dims(2)*dims(3)*dims(4));    
                    tmppres = reshape(tmppres,dims(1),dims(2)*dims(3)*dims(4));
                    tmplat  = reshape(tmplat,dims(1),dims(2)*dims(3)*dims(4));
                    [tmpbvf,tmpq,tmppav] = sw_bfrq(tmpsalt,tmptemp,tmppres,tmplat);
                    % Interpolate
                    bvf_rho = 0.5*(tmpbvf(1:end-1,:) + tmpbvf(2:end,:));
                    pv_rho  = 0.5*(tmpq(1:end-1,:) + tmpq(2:end,:));
                    % Add NaN data above and below
                    bvf_rho = [nan(1,size(bvf_rho,2));bvf_rho;nan(1,size(bvf_rho,2))];
                    pv_rho  = [nan(1,size(pv_rho,2));pv_rho;nan(1,size(pv_rho,2))];
                    % Remake
                    bvf_rho = reshape(bvf_rho,dims(1),dims(2),dims(3),dims(4));
                    bvf_rho = permute(bvf_rho,[2 3 1 4]);
                    pv_rho  = reshape(pv_rho,dims(1),dims(2),dims(3),dims(4));
                    pv_rho  = permute(pv_rho,[2 3 1 4]);
                    obj.data.(A.type).bvf.data  = bvf_rho;
                    obj.data.(A.type).bvf.name  = 'averaged Brunt Vaisala frequency';
                    obj.data.(A.type).bvf.units = 's$^{-2}$';
                    obj.data.(A.type).bvf.dims  = obj.data.(A.type).temp.units;
                    obj.data.(A.type).pv.data   = pv_rho;
                    obj.data.(A.type).pv.name   = 'averaged Potential Vorticity';
                    obj.data.(A.type).pv.units  = 'm$^{-1}$ s$^{-1}$';
                    obj.data.(A.type).pv.dims  = obj.data.(A.type).temp.units;
                    if tt == 1
                        obj.data.(A.type).temp = origtemp;
                    else
                        obj.data.(A.type).temp = [];
                    end
                    if ss == 1
                        obj.data.(A.type).salt = origsalt;
                    else
                        obj.data.(A.type).salt = [];
                    end
                    obj.data.(A.type).bvf.data = obj.data.(A.type).bvf.data .* obj.grid.(A.type).mask_rho3d;
                    obj.data.(A.type).pv.data  = obj.data.(A.type).pv.data .* obj.grid.(A.type).mask_rho3d;
                % nstar
                elseif strcmp(vars{i},'nstar');
                    disp('computeVar: Calculating N* (Deutsch et al.)');
                    obj = loadData(obj,{'NO3','NO2','PO4'},file,'type',A.type,'mean',A.mean);
                    tmpnstar = (obj.data.(A.type).NO3.data + obj.data.(A.type).NO2.data) - 16.*obj.data.(A.type).PO4.data + 2.9;
                    obj.data.(A.type).nstar.data = tmpnstar;
                    try
                        obj.data.(A.type).nstar.data = obj.data.(A.type).nstar.data .* obj.grid.(A.type).mask_rho3d;
                    catch
                        obj.data.bgc.nstar.data = obj.data.bgc.nstar.data .* obj.grid.phy.mask_rho3d;
                    end
                    obj.data.(A.type).nstar.name = 'averaged N*';
                    obj.data.(A.type).nstar.units = 'mmol N m$^{-3}$';
                    obj.data.(A.type).nstar.dims  = obj.data.(A.type).NO3.dims;
                % NPP
                elseif strcmp(vars{i},'NPP') | strcmp(vars{i},'npp');
                    disp('computeVar: Calculating NPP');
                    obj     = loadData(obj,{'TOT_PROD'},file,'type',A.type,'mean',A.mean);
                    obj     = loadDepth(obj,file,'type',A.type);
                    tmpprod = obj.data.(A.type).TOT_PROD.data;
                    tmpHz   = obj.grid.(A.type).Hz; 
                    tmpNPP  = squeeze(nansum(tmpprod.*tmpHz,3)).*3600*24*12;
                    obj.data.(A.type).NPP.data  = tmpNPP;
                    obj.data.(A.type).NPP.data  = obj.data.(A.type).NPP.data .* obj.grid.mask_rho;
                    obj.data.(A.type).NPP.name  = 'averaged Net Primary Production (NPP)';
                    obj.data.(A.type).NPP.units = 'mg C m$^{-2}$ d$^{-1}$'; 
                    obj.data.(A.type).NPP.dims  = obj.data.(A.type).TOT_PROD.dims;
                % SSH    
                elseif strcmp(vars{i},'SSH') | strcmp(vars{i},'ssh');
                    disp('computeVar: Calculating SSH');
                    obj = loadData(obj,{'zeta'},file,'type',A.type,'mean',A.mean);
                    obj = loadDiag(obj,{'SSH'},0);
                    slacorr = nanmedian(obj.diag.SSH.slice(:)) - nanmedian(obj.data.(A.type).zeta.data(:));
                    disp(' '); disp(['    Adding correction of ',num2str(slacorr),'m to ROMS SSH']);
                    obj.data.(A.type).SSH.data  = obj.data.(A.type).zeta.data + slacorr;
                    obj.data.(A.type).SSH.name  = 'averaged sea-surface height';
                    obj.data.(A.type).SSH.units = 'm'; 
                    obj.data.(A.type).SSH.dims  = obj.data.(A.type).zeta.dims; 
                % Wind Stress or Wind Stress Curl
                elseif strcmp(vars{i},'WS') | strcmp(vars{i},'ws') | strcmp(vars{i},'WSC') | strcmp(vars{i},'wsc');
                    disp('computeVar: Calculating wind-stress fields');
                    obj    = loadData(obj,{'sustr'},file,'type',A.type,'mean',A.mean);
                    obj    = loadData(obj,{'svstr'},file,'type',A.type,'mean',A.mean);
                    tmpu   = obj.data.(A.type).sustr.data;
                    tmpv   = obj.data.(A.type).svstr.data;
                    tmpang = obj.grid.angle;
                    [tmpws,tmpwsc] = romsObj.WindStress(tmpu,tmpv,obj.grid.lon_rho,obj.grid.lat_rho,tmpang);
                    obj.data.(A.type).ws.data   = tmpws;
                    obj.data.(A.type).wsc.data  = tmpwsc;
                    obj.data.(A.type).ws.data   = obj.data.(A.type).ws.data .* obj.grid.mask_rho;
                    obj.data.(A.type).wsc.data  = obj.data.(A.type).wsc.data .* obj.grid.mask_rho .* 1000;
                    obj.data.(A.type).ws.name   = 'averaged wind-stress';
                    obj.data.(A.type).wsc.name  = 'averaged wind-stress curl';
                    obj.data.(A.type).ws.units  = 'm$^{2}$ s$^{-2}$';
                    obj.data.(A.type).wsc.units = 'N m$^{-3}$';
                    obj.data.(A.type).ws.dims   = obj.data.(A.type).sustr.dims;
                    obj.data.(A.type).wsc.dims  = obj.data.(A.type).sustr.dims;
                % Mixed-layer depth (MLD)
                elseif strcmp(vars{i},'MLD') | strcmp(vars{i},'mld');
                    disp('computeVar: Calculating MLD');
                    obj = zslice(obj,{'rho'},10,file,'type',A.type,'clear','off');
                    d10  = obj.data.(A.type).rho.slice + 0.03;
                    dat  = obj.data.(A.type).rho.data;
                    dep  = obj.grid.(A.type).z_r;
                    [nx,ny,nz,nt] = size(dat);
                    lvls = zeros(nx,ny,nz,nt);
                    lvls(dat>=d10) = 1;
                    lvls = nansum(lvls,3);
                    outdep = NaN(nx,ny,nt);
                    for i = 1:nx
                        for j = 1:ny
                            for k = 1:nt
                                if lvls(i,j,k) > 0
                                    outdep(i,j,k) = dep(i,j,(lvls(i,j,k)));
                                else
                                    outdep(i,j,k) = NaN;
                                end
                            end
                        end
                    end
                    % Save ROMS data (make positive)
                    obj.data.(A.type).MLD.data  = -outdep;    
                    obj.data.(A.type).MLD.data  = obj.data.(A.type).MLD.data;
                    obj.data.(A.type).MLD.name  = 'averaged mixed-layer depth';
                    obj.data.(A.type).MLD.units = 'm';
                    obj.data.(A.type).MLD.dims  = obj.data.(A.type).rho.dims;
                    ind = find(ismember(obj.data.(A.type).MLD.dims,'s_rho')==1);
                    if ~isempty(ind)
                        obj.data.(A.type).MLD.dims(ind) = [];
                    end
                % ChlA
                elseif strcmp(vars{i},'SFC_CHL') | strcmp(vars{i},'sfc_chl');
                    disp('computeVar: Calculating SFC_CHL');
                    obj     = zslice(obj,{'TOT_CHL'},obj.grid.z_avg_dep(obj.grid.z_avg_dep<=50),file,'type',A.type);
                    tmpchla = obj.data.(A.type).TOT_CHL.slice;
                    tmpchla = squeeze(nanmean(tmpchla,3));
                    obj.data.(A.type).SFC_CHL.data  = tmpchla .* obj.grid.mask_rho;
                    obj.data.(A.type).SFC_CHL.name  = 'averaged surface chlA';
                    obj.data.(A.type).SFC_CHL.units = 'mg chlA m$^{-3}$';
                    obj.data.(A.type).SFC_CHL.dims  = obj.data.(A.type).TOT_CHL.dims;
                    ind = find(ismember(obj.data.(A.type).SFC_CHL.dims,'s_rho')==1);
                    if ~isempty(ind)
                        obj.data.(A.type).SFC_CHL.dims(ind) = [];
                    end
                % Okubo-Weiss, vorticity, sN, or sS (normal, shear of strain)
                elseif strcmp(vars{i},'OW') | strcmp(vars{i},'vort') | strcmp(vars{i},'sN') | strcmp(vars{i},'sS');
                    disp('computeVar: Calculating Okubo-Weiss');
                    if isempty(A.ip) & isempty(A.dep)
                        disp('    ERROR(computeVar): Supply a depth input via ''ip'' or ''dep''');
                        return
                    elseif ~isempty(A.ip);
                        obj = ipslice(obj,{'u','v'},A.ip,file);
                    elseif ~isempty(A.dep);
                        if ~ismember(A.dep,obj.grid.z_avg_dep);
                            disp('    ERROR(computeVar): Choose a depth from obj.grid.z_avg_dep only');
                            return
                        end
                        obj = zslice(obj,{'u','v'},A.dep,file);
                    end    
                    [OW,vort,sS,sN] = romsObj.okubo_weiss(obj.data.(A.type).u.slice,obj.data.(A.type).v.slice,obj.grid.pm,obj.grid.pn);
                    obj.data.(A.type).OW.slice   = OW;
                    obj.data.(A.type).OW.name    = 'averaged Okubo-Weiss parameter';
                    obj.data.(A.type).OW.units   = 's$^{-2}$';
                    obj.data.(A.type).OW.dims    = obj.data.(A.type).u.dims;
                    obj.data.(A.type).vort.slice = vort;
                    obj.data.(A.type).vort.name  = 'averaged relative vorticity';
                    obj.data.(A.type).vort.units = 's$^{-1}$';
                    obj.data.(A.type).vort.dims  = obj.data.(A.type).u.dims;
                    obj.data.(A.type).sN.slice   = sN;
                    obj.data.(A.type).sN.name    = 'averaged normal component of strain';
                    obj.data.(A.type).sN.units   = 's$^{-1}$';
                    obj.data.(A.type).SN.dims    = obj.data.(A.type).u.dims;
                    obj.data.(A.type).sS.slice   = sS;
                    obj.data.(A.type).sS.name    = 'averaged shear component of strain';
                    obj.data.(A.type).sS.units   = 's$^{-1}$';
                    obj.data.(A.type).sS.dims    = obj.data.(A.type).u.dims;
                    % Override
                    obj.data.(A.type).OW.dims{1}   = 'xi_rho';
                    obj.data.(A.type).vort.dims{1} = 'xi_rho'; 
                    obj.data.(A.type).sN.dims{1}   = 'xi_rho';
                    obj.data.(A.type).sS.dims{1}   = 'xi_rho';
                % Apparent Oxygen Utiliziation (AOU)
                elseif strcmp(vars{i},'AOU');
                    disp('computeVar: Calculating AOU');
                    obj = loadData(obj,{'temp','salt','O2'},file,'type',A.type,'mean',A.mean);
                    o2_sat = romsObj.o2_sat(obj.data.(A.type).temp.data,obj.data.(A.type).salt.data);
                    obj.data.(A.type).AOU.data  = o2_sat - obj.data.(A.type).O2.data;
                    obj.data.(A.type).AOU.name  = 'averaged AOU';
                    obj.data.(A.type).AOU.units = 'mmol O$_2 m$^{-3}$';
                    obj.data.(A.type).AOU.dims  = obj.data.(A.type).O2.dims;
                % Sat/Delta N2O
                elseif strcmp(vars{i},'DeltaN2O');
                    disp('computeVar: Calculating DeltaN2O');
                    obj = loadData(obj,{'temp','salt','N2O'},file,'type',A.type,'mean',A.mean);
                    n2o_sat = romsObj.n2o_sat(obj.data.(A.type).temp.data,obj.data.(A.type).salt.data);
                    obj.data.(A.type).satN2O.data    = n2o_sat;
                    obj.data.(A.type).satN2O.name    = 'averaged N$_2$O$_{sat}$';
                    obj.data.(A.type).satN2O.units   = 'mmol N$_2$O m$^{-3}$';
                    obj.data.(A.type).satN2O.dims    = obj.data.(A.type).N2O.dims;
                    obj.data.(A.type).deltaN2O.data  = obj.data.(A.type).N2O.data - n2o_sat;
                    obj.data.(A.type).deltaN2O.name  = 'averaged $\Delta$N$_2$O';
                    obj.data.(A.type).deltaN2O.units = 'mmol N$_2$O m$^{-3}$';
                    obj.data.(A.type).deltaN2O.dims    = obj.data.(A.type).N2O.dims;
                % OMZ thickness
                elseif strcmp(vars{i},'OMZ');
                    disp('computeVar: Calculating OMZ thickness');
                    obj = loadDepth(obj,file);
                    obj = loadData(obj,{'O2'},file,'type',A.type,'mean',A.mean);
                    nt = sum(cell2mat({obj.info.(['bgc_',A.type])(file).time}));
                    tmp.OMZ = nan(obj.grid.nx,obj.grid.ny,length(A.thresh),nt);
                    for th = 1:length(A.thresh)
                        tmp.O2 = obj.data.(A.type).O2.data;    
                        tmp.Hz = obj.grid.(A.type).Hz;
                        tmp.Hz(tmp.O2>A.thresh(th)) = 0;
                        tmp.OMZ(:,:,th,:) = squeeze(nansum(tmp.Hz,3));
                    end
                    obj.data.(A.type).OMZ.int = tmp.OMZ;
                    for t = 1:nt;
                        for th = 1:length(A.thresh)
                            obj.data.(A.type).OMZ.tot(th,t) = ...
                                nansum(obj.data.(A.type).OMZ.int(:,:,th,t).*obj.grid.mask_rho.* obj.grid.area_rho,'all');
                        end
                    end
                    obj.data.(A.type).OMZ.name   = 'OMZ thickness';
                    obj.data.(A.type).OMZ.units  = 'm'; 
                    obj.data.(A.type).OMZ.thresh = A.thresh;
                    obj.data.(A.type).OMZ.dims   = obj.data.(A.type).O2.dims;
                    ind = find(ismember(obj.data.(A.type).O2.dims,'s_rho')==1);
                    if ~isempty(ind)
                        obj.data.(A.type).OMZ.dims(ind) = [];
                    end    
                % NH4 / NO2 ratio
                elseif strcmp(vars{i},'NH4vNO2') | strcmp(vars{i},'NO2vNH4');
                    disp('computeVar: Calculating NH4 / NO2 ratios');
                    obj = loadData(obj,{'NH4','NO2'},file,'type',A.type,'mean',A.mean);
                    tmpnh4 = obj.data.(A.type).NH4.data;
                    tmpnh4(tmpnh4<0) = 0;
                    tmpno2 = obj.data.(A.type).NO2.data;
                    tmpno2(tmpno2<0) = 0;
                    tmpratio = tmpnh4 ./ tmpno2;
                    tmpratio(tmpratio==0)   = NaN;
                    tmpratio(tmpratio==Inf) = NaN;
                    obj.data.(A.type).NH4vNO2.data   = tmpratio;
                    obj.data.(A.type).NH4vNO2.name   = 'NH$^{+}_4$ / NO$^{-}_2$';
                    obj.data.(A.type).NH4vNO2.units  = 'ratio ($\mu$M/$\mu$M)';
                    obj.data.(A.type).NH4vNO2.dims   = obj.data.(A.type).NH4.dims;
                    tmpratio = tmpno2 ./ tmpnh4;
                    tmpratio(tmpratio==0)   = NaN;
                    tmpratio(tmpratio==Inf) = NaN;
                    obj.data.(A.type).NO2vNH4.data   = tmpratio;
                    obj.data.(A.type).NO2vNH4.name   = 'NO$^{-}_2$ / NH$^{+}_4$';
                    obj.data.(A.type).NO2vNH4.units  = 'ratio ($\mu$M/$\mu$M)';
                    obj.data.(A.type).NO2vNH4.dims   = obj.data.(A.type).NH4.dims;
                % NOX (NO2 + NO3)
                elseif strcmp(vars{i},'NOX')
                    disp('computeVar: Calculating NOx (NO3 + NO2)');
                    obj = loadData(obj,{'NO3','NO2'},file,'type',A.type,'mean',A.mean);
                    tmpno3 = obj.data.(A.type).NO3.data;
                    tmpno2 = obj.data.(A.type).NO2.data;
                    tmpno3(tmpno3<0) = 0;
                    tmpno2(tmpno2<0) = 0;
                    tmpnox = tmpno3 + tmpno2;
                    obj.data.(A.type).NOX.data  = tmpnox;
                    obj.data.(A.type).NOX.name  = 'averaged NO$_x$';
                    obj.data.(A.type).NOX.units = obj.data.(A.type).NO3.units;
                    obj.data.(A.type).NOX.dims  = obj.data.(A.type).NO3.dims;
                    obj.data.(A.type).NO3 = [];
                    obj.data.(A.type).NO2 = [];
                % SUSTR (fix for old naming style)
                elseif strcmp(vars{i},'sustr');
                    disp('computeVar: Calculating sustr');
                    obj = loadData(obj,{'sustr_avg'},file,'type',A.type,'mean',A.mean);
                    obj.data.(A.type).sustr.data  = obj.data.(A.type).sustr_avg.data;
                    obj.data.(A.type).sustr.name  = obj.data.(A.type).sustr_avg.name;
                    obj.data.(A.type).sustr.units = obj.data.(A.type).sustr_avg.units;
                    obj.data.(A.type).sustr.dims  = obj.data.(A.type).sustr_avg.dims;
                    obj.data.(A.type).sustr_avg   = [];
                % SVSTR (fix for old naming style)
                elseif strcmp(vars{i},'svstr');
                    disp('computeVar: Calculating svstr');
                    obj = loadData(obj,{'svstr_avg'},file,'type',A.type,'mean',A.mean);
                    obj.data.(A.type).svstr.data  = obj.data.(A.type).svstr_avg.data;
                    obj.data.(A.type).svstr.name  = obj.data.(A.type).svstr_avg.name;
                    obj.data.(A.type).svstr.units = obj.data.(A.type).svstr_avg.units;
                    obj.data.(A.type).svstr.dims  = obj.data.(A.type).svstr_avg.dims;
                    obj.data.(A.type).svstr_avg   = [];
                else
                    disp(['    ERROR(computeVar): ',vars{i},' is already, or cant be, calculated']);
                end
            end
        end % end method computeVar

        %--------------------------------------------------------------------------------
        function obj = sliceROMS(obj,vars,choice,deg,file,varargin);
            % -------------------
            % Takes depth slice of ROMS along a given latitude, longitude, xi-index, or eta-index
            % 
            % Usage:
            % - obj = sliceROMS(obj,vars,choice,deg,file,varargin);
            % 
            % Inputs:
            % - vars   = ROMS variable(s) to slice, as a cell array
            % - choice = 'lon','lat','xi','eta' 
            %              (lat/lon slices along a given lat/lon degree)
            %              (NOTE: lon is in 0-360 format)
            %              (xi/eta slices along a given xi or eta index, use gridView(obj) to help make selection)
            % - deg    = lon/lat degree or x/y index
            % - file   = file number to slice
            %
            % Inputs (varargin):
            % - type       = 'avg' or 'his'
            % - zdep       = depth(s) to regrid ROMS output to before slicing
            %                 NOTE: if calling 'lon' or 'lat' slices, this will default to
            %                 obj.grid.z_avg_dep levels (feel free to override them!) 
            % 
            % Example
            % - obj = sliceROMS(obj,{'temp','salt'},'lon',0,1);
            % 
            % This will slice temp and salt data from file 1 at 0 degrees longitude
            % -------------------
            disp('sliceROMS: Grabbing transects');
            % Clear slice struct
            obj.slice = [];
    
            % Check inputs
            if nargin<5
                disp('    ERROR(sliceROMS): Incorrect number of inputs, see help sliceROMS');
                return
            end    

            % Grab user inputs
            A.type = ['avg'];
            A.zdep = [];
            % If calling lat/lon, force zslice
            if strcmp(choice,'lat') | strcmp(choice,'lon');
                A.zdep = obj.grid.z_avg_dep;    
            end
            A = romsObj.parse_pv_pairs(A,varargin);

            % Use raw output data
            nz = obj.grid.nz;
            nt = sum(cell2mat({obj.info.(['phy_',A.type])(file).time}));
            ns = length(deg);

            % Check if grid is loaded
            try; obj.grid.(A.type).z_r;
            catch; obj = loadDepth(obj,file,'type',A.type);
            end
            % Load variables
            for i = 1:length(vars)
                obj = loadData(obj,vars(i),file,'type',A.type);
            end
            nt = size(obj.grid.(A.type).Hz,4);

            % X slice
            if strcmp(choice,'xi')
                order = [1 2 4 3];
                for v = 1:length(vars)
                    obj.data.(A.type).(vars{v}).slice  = squeeze(obj.data.(A.type).(vars{v}).data(deg,:,:,:));
                    slicedepth = squeeze(obj.grid.(A.type).z_r(deg,:,:,:));
                end
                slicelon = squeeze(obj.grid.lon_rho(deg,:))';
                slicelat = squeeze(obj.grid.lat_rho(deg,:))';
            % Y slice
            elseif strcmp(choice,'eta');
                order = [1 2 4 3];
                for v = 1:length(vars)
                    obj.data.(A.type).(vars{v}).slice  = squeeze(obj.data.(A.type).(vars{v}).data(:,deg,:,:));
                    slicedepth = squeeze(obj.grid.(A.type).z_r(:,deg,:,:));
                end
                slicelon = squeeze(obj.grid.lon_rho(:,deg));
                slicelat = squeeze(obj.grid.lat_rho(:,deg));
            end

            % If xi-idx or eta-idx, organize slice output and end routine
            if strcmp(choice,'xi') | strcmp(choice,'eta');
                if ns>1
                    for v = 1:length(vars)
                        obj.data.(A.type).(vars{v}).slice = permute(obj.data.(A.type).(vars{v}).slice,order);
                    end
                    obj.slice.depth = permute(slicedepth,order);
                    obj.slice.lat   = permute(repmat(slicelat,[1 1 nz nt]),[1 3 2 4]);
                    obj.slice.lon   = permute(repmat(slicelon,[1 1 nz nt]),[1 3 2 4]);
                else
                    for v = 1:length(vars)
                        obj.data.(A.type).(vars{v}).slice = permute(obj.data.(A.type).(vars{v}).slice,order);
                    end
                    obj.slice.depth = permute(slicedepth,order);
                    obj.slice.lat   = repmat(slicelat,[1 nz ns nt]);
                    obj.slice.lon   = repmat(slicelon,[1 nz ns nt]);
                end
                if strcmp(choice,'xi');
                    obj.slice.deg = obj.slice.lat;
                    obj.slice.coord = 'longitude';
                elseif strcmp(choice,'eta');
                    obj.slice.deg = obj.slice.lon;
                    obj.slice.coord = 'latitude';
                end
                obj.slice.mask  = double(~isnan(obj.slice.depth));
                obj.slice.mask(obj.slice.mask==0) = NaN;
                % Kill routine    
                obj.data.(A.type).(vars{v}).data = [];
                return
            end

            % Latitude slice    
            if strcmp(choice,'lat');
                % Interpolate raw data to latitude
                for v = 1:length(vars)
                    % Go through all xi points
                    for i = 1:size(obj.grid.lon_rho,1);
                        tmplat  = squeeze(obj.grid.lat_rho(i,:));
                        tmplon  = squeeze(obj.grid.lon_rho(i,:));
                        % Find closest points
                        idx     = sum(tmplat < deg);
                        if idx == 0 | idx == length(tmplat)
                            disp('    ERROR(sliceROMS): Slice outside ROMS domain (check latitude)')
                            kill
                        end
                        tmplat  = tmplat(idx:idx+1); 
                        tmplon  = tmplon(idx:idx+1); 
                        % 2D data (x,y)
                        if ndims(obj.data.(A.type).(vars{v}).data)==2
                            tmpdat      = squeeze(obj.data.(A.type).(vars{v}).data(i,:));
                            tmpdat      = tmpdat(idx:idx+1);
                            outdat(i)   = interp1(tmplat,tmpdat,deg);
                            outdep(i)   = NaN;
                            outdeg(i)   = interp1(tmplat,tmplon,deg);
                        % 3D data
                        elseif ndims(obj.data.(A.type).(vars{v}).data)==3
                            tmpdat = squeeze(obj.data.(A.type).(vars{v}).data(i,:,:));
                            % (x,y,z)
                            if ismember('s_rho',obj.data.(A.type).(vars{v}).dims);
                                tmpz = squeeze(obj.grid.(A.type).z_r(i,:,:));
                                dat_raw = [];
                                z_raw   = [];
                                for j = 1:obj.grid.s_rho 
                                    TMPDAT     = tmpdat(idx:idx+1,j);    
                                    TMPZ       = tmpz(idx:idx+1,j);
                                    dat_raw(j) = interp1(tmplat,TMPDAT,deg);
                                    z_raw(j)   = interp1(tmplat,TMPZ,deg);
                                end
                                % Interpolate to standard depths    
                                if length(find(~isnan([dat_raw+z_raw]))==1) == 0
                                    outdat(i,:) = nan(size(A.zdep));
                                else
                                    outdat(i,:) = interp1(-z_raw,dat_raw,A.zdep);
                                end
                                outdeg(i,:) = interp1(tmplat,tmplon,deg) .* ones(1,length(A.zdep));
                                outdep(i,:) = A.zdep; 
                            % (x,y,t)
                            else
                                for j = 1:nt
                                    TMPDAT      = tmpdat(idx:idx+1,j);    
                                    outdat(i,j) = interp1(tmplat,TMPDAT,deg);
                                    outdeg(i,j) = interp1(tmplat,tmplon,deg);
                                    outdep(i,j) = NaN;
                                end
                            end
                        % 4D data (x,y,z,t)
                        elseif ndims(obj.data.(A.type).(vars{v}).data)==4
                            for j = 1:nt
                                tmpdat = squeeze(obj.data.(A.type).(vars{v}).data(i,:,:,j));
                                tmpz   = squeeze(obj.grid.(A.type).z_r(i,:,:,j));
                                for k = 1:obj.grid.s_rho
                                    TMPDAT     = tmpdat(idx:idx+1,k);
                                    TMPZ       = tmpz(idx:idx+1,k);
                                    dat_raw(k) = interp1(tmplat,TMPDAT,deg); 
                                    z_raw(k)   = interp1(tmplat,TMPZ,deg);
                                end
                                % Interpolate to standard depths
                                if length(find(~isnan([dat_raw+z_raw]))==1) == 0
                                    outdat(i,:,j) = nan(size(A.zdep));
                                else
                                    outdat(i,:,j) = interp1(-z_raw,dat_raw,A.zdep);
                                end
                                outdep(i,:,j) = A.zdep;
                                outdeg(i,:,j) = interp1(tmplat,tmplon,deg) .* ones(1,length(A.zdep)); 
                            end
                        end
                    end 
                    % Save output
                    obj.data.(A.type).(vars{v}).slice = outdat;
                end
                % Save slice coordinates
                obj.slice.depth = outdep; 
                obj.slice.sect  = deg;
                obj.slice.coord = 'latitude';
                obj.slice.deg   = outdeg; 
            % Longitude slice
            elseif strcmp(choice,'lon');
                % Interpolate raw data to longitude
                for v = 1:length(vars)
                    % Go through all xi points
                    for i = 1:size(obj.grid.lat_rho,2);
                        tmplon  = squeeze(obj.grid.lon_rho(:,i));    
                        tmplat  = squeeze(obj.grid.lat_rho(:,i));
                        % Find closest points
                        idx     = sum(tmplon < deg);
                        if idx == 0 | idx == length(tmplon)
                            disp('    ERROR(sliceROMS): Slice outside ROMS domain (check longitude)')
                            kill
                        end
                        tmplon  = tmplon(idx:idx+1); 
                        tmplat  = tmplat(idx:idx+1); 
                        % 2D data (x,y)
                        if ndims(obj.data.(A.type).(vars{v}).data)==2
                            tmpdat      = squeeze(obj.data.(A.type).(vars{v}).data(:,i));
                            tmpdat      = tmpdat(idx:idx+1);
                            outdat(i)   = interp1(tmplon,tmpdat,deg);
                            outdep(i)   = NaN;
                            outdeg(i)   = interp1(tmplon,tmplat,deg);
                        % 3D data
                        elseif ndims(obj.data.(A.type).(vars{v}).data)==3
                            tmpdat = squeeze(obj.data.(A.type).(vars{v}).data(:,i,:));
                            % (x,y,z)
                            if ismember('s_rho',obj.data.(A.type).(vars{v}).dims);
                                tmpz = squeeze(obj.grid.(A.type).z_r(:,i,:));
                                dat_raw = [];
                                z_raw   = [];
                                for j = 1:obj.grid.s_rho 
                                    TMPDAT     = tmpdat(idx:idx+1,j);    
                                    TMPZ       = tmpz(idx:idx+1,j);
                                    dat_raw(j) = interp1(tmplon,TMPDAT,deg);
                                    z_raw(j)   = interp1(tmplon,TMPZ,deg);
                                end
                                % Interpolate to standard depths    
                                if length(find(~isnan([dat_raw+z_raw]))==1) == 0
                                    outdat(i,:) = nan(size(A.zdep)); 
                                else
                                    outdat(i,:) = interp1(-z_raw,dat_raw,A.zdep);
                                end
                                outdeg(i,:) = interp1(tmplon,tmplat,deg) .* ones(1,length(A.zdep));
                                outdep(i,:) = A.zdep; 
                            % (x,y,t)
                            else
                                for j = 1:nt
                                    TMPDAT      = tmpdat(idx:idx+1,j);    
                                    outdat(i,j) = interp1(tmplon,TMPDAT,deg);
                                    outdeg(i,j) = interp1(tmplon,tmplat,deg);
                                    outdep(i,j) = NaN;
                                end
                            end
                        % 4D data (x,y,z,t)
                        elseif ndims(obj.data.(A.type).(vars{v}).data)==4
                            for j = 1:nt
                                tmpdat = squeeze(obj.data.(A.type).(vars{v}).data(:,i,:,j));
                                tmpz   = squeeze(obj.grid.(A.type).z_r(:,i,:,j));
                                for k = 1:obj.grid.s_rho
                                    TMPDAT     = tmpdat(idx:idx+1,k);
                                    TMPZ       = tmpz(idx:idx+1,k);
                                    dat_raw(k) = interp1(tmplon,TMPDAT,deg); 
                                    z_raw(k)   = interp1(tmplon,TMPZ,deg);
                                end
                                % Interpolate to standard depths
                                if length(find(~isnan([dat_raw+z_raw]))==1) == 0
                                    outdat(i,:,j) = nan(size(A.zdep)); 
                                else
                                    outdat(i,:,j) = interp1(-z_raw,dat_raw,A.zdep);
                                end
                                outdep(i,:,j) = A.zdep;
                                outdeg(i,:,j) = interp1(tmplon,tmplat,deg) .* ones(1,length(A.zdep));
                            end
                        end
                    end 
                    % Save output
                    obj.data.(A.type).(vars{v}).slice = outdat;
                end
                % Save slice coordinates
                obj.slice.depth = outdep; 
                obj.slice.sect  = deg;
                obj.slice.coord = 'longitude';
                obj.slice.deg   = outdeg;     
            end
        end % end method sliceROMS

        %--------------------------------------------------------------------------------
        function obj = zslice(obj,vars,zdep,file,varargin)
            % ------------------
            % Slices ROMS along a given depth 
            % Use obj.grid.z_avg_dep for WOA18 depths
            % 
            % Usage:
            % - obj = zslice(obj,vars,zdep,file);
            % 
            % Inputs:
            % - vars = cell array of ROMS variables to slice
            % - zdep = depth(s) to slice along
            % - file = file to slice
            %
            % Optional:
            % - type  = 'avg','his','rst',etc
            % - clear = 'on' to clear raw data after slicing ('off' otherwise) 
            %
            % Example:
            % - obj = zslice(obj,{'O2'},[0 100 200],file);
            % -------------------    
            disp('zslice: Grabbing zsliced variables');
            
            % Process optional inputs
            A.type  = 'avg'; % file type
            A.clear = 'on';  % if on, clear raw data after calculation
            A       = romsObj.parse_pv_pairs(A,varargin);

            % Compute depth ('full' if requesting u/v points) 
            if ismember('u',vars) | ismember('v',vars);
                try; obj.grid.(['phy_',A.type]).z_u;
                catch; obj = loadDepth(obj,file,'type',A.type,'full',1);
                end
            else
                try; obj.grid.(['phy_',A.type]).z_r;
                catch; obj = loadDepth(obj,file,'type',A.type);
                end
            end

            % Change '0' to 5
            zdep(zdep==0) = 5;

            % Load variables
            for i = 1:length(vars)
                try
                    obj = loadData(obj,vars(i),file,'type',A.type);
                catch
                    obj = computeVar(obj,vars(i),file,'type',A.type);
                end
            end

            % Get time dimension
            nt = sum(cell2mat({obj.info.(['phy_',A.type])(file).time}));

            % Interpolate variable to zdep level
            for i = 1:length(vars)
                for t = 1:nt
                    vnew   = nan(obj.grid.nx,obj.grid.ny,length(zdep));
                    tmpvar = obj.data.(A.type).(vars{i}).data(:,:,:,t);
                    if strcmp(vars{i},'u');
                        tmpz = -obj.grid.(A.type).z_u(:,:,:,t);
                        vnew = vnew(1:end-1,:,:);
                    elseif strcmp(vars{i},'v');
                        tmpz = -obj.grid.(A.type).z_v(:,:,:,t);
                        vnew = vnew(:,1:end-1,:);
                    else
                        tmpz = -obj.grid.(A.type).z_r(:,:,:,t);
                    end
                    for z = 1:length(zdep)        
                        % Get sizes
                        [lx,ly,lz] = size(tmpvar);
                        N = lx*ly;
            
                        % Convert to 2D                
                        tmp.z = permute(tmpz,[1 2 3]);
                        tmp.var = permute(tmpvar,[1 2 3]);
                        tmp.z = reshape(tmp.z,[lx*ly lz]);
                        tmp.var = reshape(tmp.var,[lx*ly lz]);

                        % Find level where z_r > z
                        a=tmp.z>zdep(z);
                        levs=squeeze(sum(a,2));
                        levs(levs==lz)=lz;
                        levs(levs==0)=1;
                        mask=levs./levs;
                        mask(isnan(mask)) = 0;
                        mask(levs==1) = 0;
                        mask(levs==lz) = 0;
                        tmp.vnew = nan(size(levs));

                        % Go through and interpolate
                        for n = 1:(lx*ly)
                            if mask(n)==0
                                tmp.vnew(n) = NaN;
                                continue;
                            else
                                % Grab upper and lower levels
                                z1 = tmp.z(n,levs(n)+1);
                                z2 = tmp.z(n,levs(n));
                                v1 = tmp.var(n,levs(n)+1);
                                v2 = tmp.var(n,levs(n));
                                % Interpolate
                                tmp.vnew(n,1) = [v2.*(z1-zdep(z)) + v1.*(zdep(z)-z2)]./[(z1-z2)];
                            end
                        end
                        vnew(:,:,z) = reshape(tmp.vnew,[lx,ly]);
                    end
                    obj.data.(A.type).(vars{i}).slice(:,:,:,t)  = squeeze(vnew);
                end
                % Replace data
                if strcmp(A.clear,'on');
                    obj.data.(A.type).(vars{i}).data   = [];
                end
            end

            % Save into grid
            obj.slice.zdep  = zdep;
            if nt == 1
                obj.slice.depth = permute(repmat(zdep,[1 lx ly]),[2 3 1]);
            else
                [a,b] = size(zdep);
                if a > b
                    obj.slice.depth = permute(repmat(zdep,[1 lx ly nt]),[2 3 1 4]);
                else
                    obj.slice.depth = permute(repmat(zdep',[1 lx ly nt]),[2 3 1 4]);
                end
            end
        end % end method zslice

        %--------------------------------------------------------------------------------
        function obj = ipslice(obj,vars,ip,file,varargin)
            % ------------------
            % Slices ROMS along a given isopycnal
            % 
            % Usage:
            % - obj = ipslice(obj,vars,ip,file);
            % 
            % Inputs:
            % - vars = cell array of ROMS variables to slice
            % - ip   = isopycnal(s) to slice along
            % - file = file to slice
            %
            % Optional:
            % - type = 'avg','his','rst',etc
            %
            % Example:
            % - obj = ipslice(obj,{'O2'},26.5,file);
            % -------------------    
            disp('ipslice: Grabbing ipsliced variables');

            % Process optional inputs
            A.type = 'avg'; % file type
            A      = romsObj.parse_pv_pairs(A,varargin);

            % Load rho
            obj = loadData(obj,{'rho'},file,'type',A.type);

            % Load variables
            for i = 1:length(vars)
                try
                    obj = loadData(obj,vars(i),file,'type',A.type);
                catch
                    obj = computeVar(obj,vars(i),file,'type',A.type);
                end
            end

            % Get time dimension
            nt = sum(cell2mat({obj.info.(['phy_',A.type])(file).time}));

            % Interpolate variable to ip level
            for i = 1:length(vars)
                for t = 1:nt
                    vnew   = nan(obj.grid.nx,obj.grid.ny,length(ip));
                    tmpvar = obj.data.(A.type).(vars{i}).data(:,:,:,t);
                    tmpip  = obj.data.(A.type).rho.data(:,:,:,t);
                    if strcmp(vars{i},'u')
                        tmpip = 0.5.*(tmpip(1:end-1,:,:) + tmpip(2:end,:,:));
                        vnew  = vnew(1:end-1,:);
                    elseif strcmp(vars{i},'v');
                        tmpip = 0.5.*(tmpip(:,1:end-1,:) + tmpip(:,2:end,:));
                        vnew  = vnew(:,1:end-1);
                    end
                    for z = 1:length(ip);
                        % Get sizes
                        [lx,ly,lz] = size(tmpvar);
                        N = lx*ly;

                        % Convert to 2D
                        tmp.ip = permute(tmpip,[1 2 3]);
                        tmp.var = permute(tmpvar,[1 2 3]);
                        tmp.ip = reshape(tmp.ip,[lx*ly lz]);
                        tmp.var = reshape(tmp.var,[lx*ly lz]);

                        % Find level where ip > ip(z)
                        a=tmp.ip>ip(z);
                        levs=squeeze(sum(a,2));
                        levs(levs==lz)=lz;
                        levs(levs==0)=1;
                        mask=levs./levs;
                        mask(isnan(mask)) = 0;
                        mask(levs==1) = 0;
                        mask(levs==lz) = 0;
                        tmp.new = nan(size(levs));

                        % Go through and interpolate    
                        for n = 1:(lx*ly)
                            if mask(n)==0
                                tmp.vnew(n) = NaN;
                                continue    
                            else                        
                                % Grab upper and lower levels
                                z1 = tmp.ip(n,levs(n)+1);
                                z2 = tmp.ip(n,levs(n));
                                v1 = tmp.var(n,levs(n)+1);
                                v2 = tmp.var(n,levs(n));
                                % Interpolate
                                tmp.vnew(n) = [v2.*(z1-ip(z)) + v1.*(ip(z)-z2)]./[(z1-z2)];
                            end
                        end
                        vnew(:,:,z) = reshape(tmp.vnew,[lx,ly]);
                    end
                    % Replace data
                    obj.data.(A.type).(vars{i}).slice(:,:,:,t) = squeeze(vnew);
                end
                obj.data.(A.type).(vars{i}).data  = [];
            end

            % Save into grid
            obj.slice.ip = ip;
            if nt == 1
                obj.slice.rho = repmat(ip,[lx ly]);
            else
                obj.slice.rho = repmat(ip,[lx ly nt]);
            end
        end % end method ipslice

        %--------------------------------------------------------------------------------
        function obj = intVar(obj,vars,varargin)
            % ------------------
            % Vertically integrate 3D variable(s) 
            %
            % Usage:
            % - obj = intVar(obj,vars,varargin)
            %
            % Inputs:
            % - vars = 3D variables to integrate, as a cell array
            %
            % Optional:
            % - type = file type ('avg','his', etc)
            %
            % Example:
            % - obj = loadData(obj,{'NO3'},file);
            % - obj = intVar(obj,{'NO3'});
            % ------------------
            disp('intVar: Integrating 3D variables');

            % Process optional inputs
            A.type = 'avg'; % file type
            A      = romsObj.parse_pv_pairs(A,varargin);

            % Load depths?
            try; obj.grid.(A.type).Hz;
            catch; disp('    ERROR(intVar): call loadDepth first'); return
            end

            % Go through each 3D rate, integrate vertically (and totally)
            for i = 1:length(vars)        
                tmpdata = obj.data.(A.type).(vars{i}).data .* obj.grid.(A.type).Hz .* obj.grid.(A.type).mask_rho3d; % mmol/m3/s --> mmol/m2/s
                tmpdata = squeeze(nansum(tmpdata,3));
                tmpdata = tmpdata .* obj.grid.mask_rho;
                obj.data.(A.type).(vars{i}).int = tmpdata; 
                obj.data.(A.type).(vars{i}).intunits = strrep(obj.data.(A.type).(vars{i}).units,'m$^{-3}$','m$^{-2}$');
                for t = 1:size(obj.grid.(A.type).Hz,4);
                    obj.data.(A.type).(vars{i}).tot(t) = nansum(obj.data.(A.type).(vars{i}).int(:,:,t) .*obj.grid.area_rho,'all');
                end
                obj.data.(A.type).(vars{i}).totunits = strrep(obj.data.(A.type).(vars{i}).units,' m$^{-3}$','');
            end
        end % end method intVar

        %--------------------------------------------------------------------------------
        function obj = getProfile(obj,vars,lon,lat,file,varargin)
            % ------------------
            % Loads profile data at the nearest lon/lat point 
            %
            % Usage:
            %    - obj = getProfile(obj,vars,lon,lat,file);
            %
            % Inputs:
            %    - vars  = variables to load, as a cell array
            %    - lon   = vector of longitude points
            %    - lat   = vector of latitude points
            %    - file  = input file number
            %
            % Inputs (varargin):
            %    - type     = 'avg' or 'z_avg' (raw cant be used with yr_range)
            %    - yr_range = range of years to load and average (i.e. 2045:2049)
            %
            % Example:
            %    - obj = getProfile(obj,{'temp','salt','O2','NO3'},[250 250],[-15 -20]);
            % -------------------
            disp('getProfile: Grabbing profile(s)');

            % process inputs
            A.type     = 'avg';
            A          = romsObj.parse_pv_pairs(A,varargin);    

            % Get indices of nearest point
            lon_idx = [];
            lat_idx = [];
            for i = 1:length(lon);
                % Get index of nearest lat/lon grid cell
                all_dist                = distance(lat(i),lon(i),obj.grid.lat_rho,obj.grid.lon_rho);
                [lon_idx(i),lat_idx(i)] = find(all_dist == min(all_dist(:)));
            end
               
            % Load data
            for i = 1:length(vars)
                try
                    obj = loadData(obj,vars(i),file,'type',A.type);
                catch
                    obj = computeVar(obj,vars(i),file,'type',A.type);
                end
            end
            % Load depth
            if isempty(obj.grid.(A.type))
                obj = loadDepth(obj,file);
            end
                
            % Go through all variables and lon/lats
            for v = 1:length(vars)
                tmp.data = obj.data.(A.type).(vars{v}).data;    
                % Save profile data
                for i = 1:length(lon);
                    obj.data.(A.type).(vars{v}).profile(i,:,:)  = squeeze(tmp.data(lon_idx(i),lat_idx(i),:,:));
                    obj.profile.depth(i,:,:) = squeeze(obj.grid.(A.type).z_r(lon_idx(i),lat_idx(i),:,:));
                    obj.profile.lon(i)  = obj.grid.lon_rho(lon_idx(i),lat_idx(i));
                    obj.profile.lat(i)  = obj.grid.lat_rho(lon_idx(i),lat_idx(i));
                    obj.profile.xidx(i) = lon_idx(i);
                    obj.profile.yidx(i) = lat_idx(i);
                end
            end
        end % end method getProfile

        %--------------------------------------------------------------------------------
        function dispVars(obj,file_type,display)
            % ----------------------
            % Lists output variables in command window
            % Assumes all files are similar to the first file
            %
            % Usage:
            % dispVars(obj,file_type,display)    
            %
            % Inputs:
            % - type    = 'phy_avg','bgc_his','dia_avg', etc. 
            % - display = any number will show full ncdisp output (default is 'min')
            %
            % Example:
            % - dispVars(obj,'type','his');
            % ----------------------
            if nargin<2 
                disp('    ERROR(dispVars): Not enough inputs, see help dispVars');
                disp('    Choose from: ');
                if ~isempty(obj.info.phy_avg);
                    disp('    phy_avg');
                end
                if ~isempty(obj.info.phy_his);
                    disp('    phy_his');
                end
                if ~isempty(obj.info.bgc_avg);
                    disp('    bgc_avg');
                end
                if ~isempty(obj.info.bgc_his);
                    disp('    bgc_his');
                end
                if ~isempty(obj.info.dia_avg);
                    disp('    dia_avg');
                end
                if ~isempty(obj.info.dia_his);
                    disp('    dia_his');
                end
                if ~isempty(obj.info.flux_avg);
                    disp('    flux_avg');
                end
                if ~isempty(obj.info.flux_his);
                    disp('    flux_his');
                end
                if ~isempty(obj.info.flx_avg);
                    disp('    flx_avg');
                end
                if ~isempty(obj.info.flx_his);
                    disp('    flx_his');
                end
                return
            end
 
            % Call ncdisp    
            if nargin < 3
                ncdisp([obj.info.(file_type)(1).Filename],'/','min');    
            else
                ncdisp([obj.info.(file_type)(1).Filename]);    
            end
        end % end method dispVars

        %--------------------------------------------------------------------------------
        function obj = getBudg(obj,vars,file,varargin)
            % --------------------
            % Main method to perform budget analysis on a ROMS tracer (vars).
            %
            % Usage:
            % - obj = getBudg(obj,vars,file,varargin)
            %
            % Inputs:
            % - vars = cell array of budget(s) to close 
            % - file = (#s) load specific time-averages 
            %        = 0 load all files and average (for monthly)
            %
            % Optional Inputs:
            % - hisfile  = alternate file#s for history files 
            % - physfile = alternate file#s for flux files 
            % - srcsnk   = (1) to call sourcesSinks (default 0)
            % - int      = (1) to call intBudg (default 0)
            % - clean    = (0) to keep all data loaded (default 1)
            %
            % Example:
            % - obj = getBudg(obj,'N2O',1)
            % --------------------

            % Process optional inputs (varargin)
            A.int      = [0];
            A.srcsnk   = [0];
            A.clean    = [1];
            A.hisfile  = [];
            A.physfile = [];
            A          = romsObj.parse_pv_pairs(A,varargin);

            % Check inputs
            if isempty(vars)
                disp('vars must be defined, see help getBudg')
                return
            end
            if isempty(A.hisfile)
                hisfile = file;
            else
                hisfile = A.hisfile;
            end
            if isempty(A.physfile)
                physfile = file;
            else
                physfile = A.physfile;
            end

            % Get budget
            romsOpt
            for i = 1:length(vars)
                disp(['getBudg: Computing budget for ',vars{i}]);

                % Grab budget settings
                if ~isempty(budget.(vars{i}))
                    obj.budget.(vars{i}).info = budget.(vars{i});
                else
                    disp('    ERROR(getBudg): budget settings not set in romsOpt.m');
                end
            
                % Load all 3D output terms
                terms = [vars{i},obj.budget.(vars{i}).info.rates,obj.budget.(vars{i}).info.fluxes];
                terms = [terms(find(~cellfun(@isempty,terms)))];
                obj   = loadDepth(obj,file,'full',1);
                obj   = loadData(obj,terms,file,'type','avg');

                % Integrate 3D variables and rates vertically
                terms = [vars{i},obj.budget.(vars{i}).info.rates];
                terms = [terms(find(~cellfun(@isempty,terms)))];
                obj   = intVar(obj,terms);

                % Load and process 2D fluxes
                obj = getFluxes(obj,vars(i),...
                    obj.budget.(vars{i}).info.fluxes,...
                    obj.budget.(vars{i}).info.lvls,...
                    obj.budget.(vars{i}).info.feq);

                % Get dCdt, advection, sms, net terms
                obj = computeDcDt(obj,vars(i),hisfile);
                obj = computeAdv(obj,vars(i),physfile);
                obj = computeSMS(obj,vars(i),...
                    obj.budget.(vars{i}).info.rates,...
                    obj.budget.(vars{i}).info.smseq,...
                    obj.budget.(vars{i}).info.rtits);
                obj = computeNet(obj,vars(i));

                % OPTIONAL SWITCHES
                % Integrate vertically (and horizontally)
                if A.int == 1
                    obj = intBudg(obj,vars(i)); 
                end
                % Split sources and sinks
                if A.srcsnk == 1
                    obj = sourcesSinks(obj,vars(i));
                end
                % Save average, clear loaded data
                if A.clean == 1
                    tmp = obj.data.avg.(vars{i});
                    obj.data = [];
                    obj.data.avg.(vars{i}) = tmp;
                end
            end
        end % end method getBudg

        %--------------------------------------------------------------------------------
        function obj = getFluxes(obj,vars,fluxes,lvls,feq,varargin);
            % -------------------
            % Grab 2D fluxes for budget, convert to 3D 
            % Called in getBudg
            % Fluxes in mmol/m2/s
            %
            % Usage:
            % - obj = getFluxes(obj,vars,fluxes,lvls,feq)
            %
            % Inputs:
            % - vars    = BGC budget that you are closing (i.e. 'NO2')
            % - fluxes  = 2D fluxes (air-sea, sediment, etc)
            % - lvls    = levels to apply 2D flux ('sfc','sed', as a cell array)
            % - feq     = factors to multiply fluxes by
            %
            % Optional inputs:
            % - type    = file type (avg, his, rst)
            %
            % Example:
            % - obj = getFluxes(obj,'N2O',{'FG_N2O'},{'sfc'},[1])
            % -------------------
            disp('getFluxes: Get 2D interface fluxes, convert to 3D')

            % Process optional inputs
            A.type = 'avg'; % file type
            A      = romsObj.parse_pv_pairs(A,varargin);
 
            % Initialize matrices-to-fill
            obj.budget.(vars{1}).fg  = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
            obj.budget.(vars{1}).sed = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;           

            % Kill if no fluxes
            if isempty(fluxes{1});
                return
            end

            % Convert 2D flux to 3D based on lvls
            for i = 1:length(fluxes)
            
                % Apply 2D mask to 2D data, apply any factors(feq)
                obj.data.(A.type).(fluxes{i}).data = obj.data.(A.type).(fluxes{i}).data .* feq(i) .* obj.grid.mask_rho;
                
                % Apply 2D flux to correct z-level to make 3D
                if strcmp(lvls{i},'sfc')
                    tmpfg = zeros(size(obj.grid.(A.type).mask_rho3d));
                    % Apply value into 3D grid
                    tmpfg(:,:,obj.grid.nz,:) = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
                    % Divide by z, save as 3D rate
                    obj.budget.(vars{1}).fg = tmpfg ./ obj.grid.(A.type).Hz;
                    % Mask padding (due to XYZ flux output)
                    obj.budget.(vars{1}).fg(end,:,:,:) = nan;
                    obj.budget.(vars{1}).fg(:,end,:,:) = nan;
                elseif strcmp(lvls{i},'sed')
                    tmpsed = zeros(size(obj.grid.(A.type).mask_rho3d));
                    % Apply value into 3D grid
                    tmpsed(:,:,1,:) = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
                    % Divide by z, save as 3D rate
                    obj.budget.(vars{1}).sed = tmpsed ./ obj.grid.(A.type).Hz;
                    % Mask padding (due to XYZ flux output)
                    obj.budget.(vars{1}).sed(end,:,:,:) = nan;
                    obj.budget.(vars{1}).sed(:,end,:,:) = nan;
                end
            end
        end % end method getFluxes

        %--------------------------------------------------------------------------------
        function obj = computeDcDt(obj,vars,file)
            % -------------------
            % Compute change in concentration with time
            % Called in getBudg
            % End result is mmol/m3/s
            %
            % Usage:
            % - obj = computeDcDt(obj,vars,file)
            %
            % Inputs:
            % - vars = variable to calculated dC/dt for
            % - file = (#s) load specific time-averages 
            %        = 0 load all files and average (for monthly)
            %
            % Example:
            % - obj = computeDcDt(obj,'NO2',1);
            % -------------------
            disp('computeDcDt: Grabbing tracer tendencies over averaging period')

            % Compute dt according to output frequency
            dt1 = [];
            dt2 = [];
            try
                dt1 = double(obj.info.params.dt)*double(obj.info.params.nwrt);
            catch
                dt1 = 0; % ERROR
            end
            try
                dt2 = double(obj.info.params.dt)*double(obj.info.params.navg);
            catch
                dt2 = 0; % ERROR
            end
            if dt1==dt2
                dt = dt1;
            elseif dt1>dt2
                dt = dt2; % Sometimes one is written to fill_values
            elseif dt2>dt1
                dt = dt1; % Sometimes one is written to fill_values
            end
            if dt == 0
                disp('    ERROR(computeDcDt): Cant find number of timesteps used in averaging')
                kill
            end

            % Load data on either side of 'history' (first,last snapshot)
            for v = 1:length(vars)
                % Get start/finish from history
                obj  = loadData(obj,vars(v),file,'type','his');        
                idx  = find(strcmp('time',obj.data.his.(vars{v}).dims)==1);
                if idx == 3 % 2D field
                    obj.budget.(vars{v}).dcdt = ((obj.data.his.(vars{v}).data(:,:,2:end) - ...
                                                 obj.data.his.(vars{v}).data(:,:,1:end-1))./dt) .* obj.grid.mask_rho;
                elseif idx == 4 % 3D field
                    obj.budget.(vars{v}).dcdt = ((obj.data.his.(vars{v}).data(:,:,:,2:end) - ...
                                                 obj.data.his.(vars{v}).data(:,:,:,1:end-1))./dt) .* obj.grid.avg.mask_rho3d;
                end

                % Mask padding (for XYZ output)
                obj.budget.(vars{v}).dcdt(end,:,:,:) = nan;
                obj.budget.(vars{v}).dcdt(:,end,:,:) = nan;
            end 
        end % end method computeDcDt

        %--------------------------------------------------------------------------------
        function obj = computeAdv(obj,vars,file)
            % -------------------
            % Compute HorXAdvFlux, HorYAdvFlux, and top/bottom advection
            % Also get diffusion, if it is available
            % Called in getBudg
            % End result is mmol/m3/s
            %
            % Usage:
            % - obj = computeAdv(obj,vars)
            %
            % Inputs:
            % - vars = variable to get flux terms for, as a cell array
            %
            % Example:
            % - obj = computeAdv(obj,{'NO2'});
            % -------------------
            disp('computeAdv: Get advective/diffusion terms in budget equation')
        
            % Cycle through and load XYZfluxes
            for i = 1:length(vars)
                % Load XYZ advection 
                obj = loadData(obj,{['HorXAdvFlux_',vars{i}]},file);  % units of mmol/s
                obj = loadData(obj,{['HorYAdvFlux_',vars{i}]},file);  % units of mmol/s
                obj = loadData(obj,{['VertAdvFlux_',vars{i}]},file);  % units of mmol/m2s
                obj = loadData(obj,{['VertDiffFlux_',vars{i}]},file); % units of mmol/m2s

                % Get temporary arrays
                tmp.X  = obj.data.avg.(['HorXAdvFlux_',vars{i}]).data;
                tmp.Y  = obj.data.avg.(['HorYAdvFlux_',vars{i}]).data;
                tmp.Z  = obj.data.avg.(['VertAdvFlux_',vars{i}]).data;
                tmp.Zd = obj.data.avg.(['VertDiffFlux_',vars{i}]).data;
                   
                % Add extra lons/lats for calc
                pad_x = nan(1,size(tmp.X,2),size(tmp.X,3),size(tmp.X,4));
                pad_y = nan(size(tmp.Y,1),1,size(tmp.Y,3),size(tmp.Y,4));
                tmp.X = cat(1,tmp.X,pad_x);
                tmp.Y = cat(2,tmp.Y,pad_y);

                % Initiate matrices-to-fill, compute adv X and Y
                % Converts from mmol/s to mmol/m3/s by dividing by grid area and height of cell (volume)
                % Get dimensions
                nx = obj.grid.nx;
                ny = obj.grid.ny;
                nz = obj.grid.nz;
                idx = find(strcmp('time',obj.data.avg.(['HorXAdvFlux_',vars{i}]).dims)==1);
                if isempty(idx)
                    nt = 1;
                else
                    nt = size(obj.data.avg.(['HorXAdvFlux_',vars{i}]).data,idx);
                end
                
                % X advection
                adx(1:nx,1:ny,1:nz,1:nt) = NaN;
                adx(1:nx,:,:,1:nt)       = (tmp.X(1:nx,:,:,:) - tmp.X(2:nx+1,:,:,:)) ./ ...
                                            obj.grid.avg.area3d(1:nx,:,:,:) ./ obj.grid.avg.Hz(1:nx,:,:,:);
                
                % Y advection
                ady(1:nx,1:ny,1:nz,1:nt) = NaN;         
                ady(:,1:ny,:,:)          = (tmp.Y(:,1:ny,:,:) - tmp.Y(:,2:ny+1,:,:)) ./ ...
                                            obj.grid.avg.area3d(:,1:ny,:,:) ./ obj.grid.avg.Hz(:,1:ny,:,:);
                
                % Z advection
                adz(1:nx,1:ny,1:nz,1:nt) = NaN;         
                adz(:,:,:,:)             = (tmp.Z(:,:,1:nz,:) - tmp.Z(:,:,2:nz+1,:)) ./ obj.grid.avg.Hz(:,:,:,:);
                
                % Z diffusion
                dfz(1:nx,1:ny,1:nz,1:nt) = NaN;         
                dfz(:,:,:,:)             = (tmp.Zd(:,:,1:nz,:) - tmp.Zd(:,:,2:nz+1,:)) ./ obj.grid.avg.Hz(:,:,:,:);
                
                % Apply 3D mask
                obj.budget.(vars{i}).adx = adx .* obj.grid.avg.mask_rho3d(:,:,:,:);
                obj.budget.(vars{i}).ady = ady .* obj.grid.avg.mask_rho3d(:,:,:,:);
                obj.budget.(vars{i}).adz = adz .* obj.grid.avg.mask_rho3d(:,:,:,:);
                obj.budget.(vars{i}).dfz = dfz .* obj.grid.avg.mask_rho3d(:,:,:,:);

                % Mask padding (greatest x/y indicies are forced to be NaN)
                obj.budget.(vars{i}).adx(end,:,:,:) = nan; % mask final 'x' indices
                obj.budget.(vars{i}).ady(end,:,:,:) = nan; % mask final 'x' indices
                obj.budget.(vars{i}).adz(end,:,:,:) = nan; % mask final 'x' indices
                obj.budget.(vars{i}).dfz(end,:,:,:) = nan; % mask final 'x' indices
                obj.budget.(vars{i}).adx(:,end,:,:) = nan; % mask final 'y' indices
                obj.budget.(vars{i}).ady(:,end,:,:) = nan; % mask final 'y' indices
                obj.budget.(vars{i}).adz(:,end,:,:) = nan; % mask final 'y' indices
                obj.budget.(vars{i}).dfz(:,end,:,:) = nan; % mask final 'y' indices

                % Also calculate horizontal and total advection
                obj.budget.(vars{i}).adxy = obj.budget.(vars{i}).adx + obj.budget.(vars{i}).ady;
                obj.budget.(vars{i}).adv  = obj.budget.(vars{i}).adx + obj.budget.(vars{i}).ady + obj.budget.(vars{i}).adz;

                % Clear object data structure memory 
                obj.data.avg.(['HorXAdvFlux_',vars{i}])  = [];
                obj.data.avg.(['HorYAdvFlux_',vars{i}])  = [];
                obj.data.avg.(['VertAdvFlux_',vars{i}])  = [];
                obj.data.avg.(['VertDiffFlux_',vars{i}]) = [];
            end
        end % end method computeXYZFlux

        %--------------------------------------------------------------------------------
        function obj = computeSMS(obj,vars,rates,smseq,titles,varargin)
            % ---------------------
            % Gathers sources and sinks
            % Called in getBudg
            % End result is mmol/m3/s
            %
            % Usage:
            % - obj = computeSMS(obj,rates,smseq)
            %
            % Inputs:
            % - vars   = budget variable (i.e. 'NO2')
            % - rates  = BGC rates, set in getBudg
            % - smseq  = S-M-S equation, set in getBudg
            % - titles = rate titles, set in getBudg
            %
            % Optional:
            % - type = file type, defaults to 'avg' for budget calcs
            %
            % Example:
            % - vars   = {'N2O_AO1'};
            % - rates  = {'N2OAMMOX','N2OAO1_CONS');
            % - smseq  = [(1) (-1)];
            % - titles = {'NH$^{+}_4$ oxidation to N$_2$O','N$_2$O$_{nit}$ reduction'};
            % - obj = computeSMS(obj,vars,rates,smseq);
            % ---------------------
            disp('computeSMS: Computing sources-minus-sinks (SMS)');

            % Process optional inputs
            A.type = 'avg'; % file type
            A      = romsObj.parse_pv_pairs(A,varargin);

            % Write equation
            eq = [];
            smseq_orig = smseq;
            for i = 1:length(smseq)
                if smseq(i) > 0
                    sym = ' + ';
                elseif smseq(i) < 0 
                    sym = ' - ';
                    smseq(i) = -smseq(i);
                end
                if i == 1
                    eq = [sym,'(',num2str(smseq(i)),')*',rates{i}];
                else
                    eq = [eq,sym,'(',num2str(smseq(i)),')*',rates{i}];
                end
            end
            smseq = smseq_orig;
            obj.budget.(vars{1}).info.sms_eq = eq;
            obj.budget.(vars{1}).info.sms_factors = smseq;
            obj.budget.(vars{1}).info.sms_titles = titles;

            % Get nt
            idx = find(strcmp('time',obj.data.avg.(vars{1}).dims)==1);
            if isempty(idx)
                nt = 1
            else
                nt = size(obj.data.avg.(vars{1}).data,idx);
            end

            % Get SMS 
            dims = ndims(obj.data.(A.type).(rates{1}).data);
            eq = [];
            for i = 1:length(rates)
                eq{i} = [(smseq(i).*obj.data.(A.type).(rates{i}).data)];
            end
            tmpsms = sum(cat(dims+1,eq{:}),dims+1);
            obj.budget.(vars{1}).sms = tmpsms;

            % Get production 
            ind = find(smseq>0);
            if ~isempty(ind);
                tmprates = rates(ind);
                obj.budget.(vars{1}).info.prod_sources = tmprates;
                obj.budget.(vars{1}).info.prod_factors = smseq(ind);
                obj.budget.(vars{1}).info.prod_titles  = titles(ind);
                eq = [];
                for i = 1:length(tmprates)
                    eq{i} = [(smseq(ind(i)).*obj.data.(A.type).(tmprates{i}).data)];
                end
                tmpprod = sum(cat(dims+1,eq{:}),dims+1);
                obj.budget.(vars{1}).prod = tmpprod;
            else
                obj.budget.(vars{1}).info.prod_sources = [];
                obj.budget.(vars{1}).info.prod_factors = []; 
                obj.budget.(vars{1}).info.prod_titles  = [];
                obj.budget.(vars{1}).prod = zeros([obj.grid.ndim_xyz nt]);
            end

            % Get consumption 
            ind = find(smseq<0);
            if ~isempty(ind)
                tmprates = rates(ind);
                obj.budget.(vars{1}).info.cons_sources = tmprates;
                obj.budget.(vars{1}).info.cons_factors = smseq(ind);
                obj.budget.(vars{1}).info.cons_titles  = titles(ind);
                eq = [];
                for i = 1:length(tmprates)
                    eq{i} = [(smseq(ind(i)).*obj.data.(A.type).(tmprates{i}).data)];
                end
                tmpcons = sum(cat(dims+1,eq{:}),dims+1);
                obj.budget.(vars{1}).cons = tmpcons;
            else
                obj.budget.(vars{1}).info.cons_sources = [];
                obj.budget.(vars{1}).info.cons_factors = []; 
                obj.budget.(vars{1}).info.cons_titles  = []; 
                obj.budget.(vars{1}).cons = zeros([obj.grid.ndim_xyz nt]);
            end

            % mask padding (due to XYZ flux output)
            obj.budget.(vars{1}).sms(end,:,:,:) = nan;
            obj.budget.(vars{1}).prod(end,:,:,:) = nan;
            obj.budget.(vars{1}).cons(end,:,:,:) = nan;
            obj.budget.(vars{1}).sms(:,end,:,:) = nan;
            obj.budget.(vars{1}).prod(:,end,:,:) = nan;
            obj.budget.(vars{1}).cons(:,end,:,:) = nan;
        end % end method computeSMS

        %--------------------------------------------------------------------------------
        function [obj] = sourcesSinks(obj,vars)
            % ------------------------------------------------
            % Function to extract individual sources and sinks
            % Rates are in the correct units (see obj.budget.(vars).info)
            %
            % Usage: 
            % - [obj] = sourcesSinks(obj,vars) 
            % 
            % Inputs:
            % - obj:    romsObj object
            % - vars:   budget variable
            %
            % Example:
            % - [obj] = sourcesSinks(obj,{'NH4'});
            % ------------------------------------------------

            % Get sources
            for i = 1:length(vars)
                for j = 1:length(obj.budget.(vars{i}).info.prod_sources);
                    this_var = obj.budget.(vars{i}).info.prod_sources{j};
                    this_factor = obj.budget.(vars{i}).info.prod_factors(j);
                    tmp.sources.(this_var) = obj.data.avg.(this_var).data.*this_factor;
                end
                for j = 1:length(obj.budget.(vars{i}).info.cons_sources);
                    this_var = obj.budget.(vars{i}).info.cons_sources{j};
                    this_factor = obj.budget.(vars{i}).info.cons_factors(j);
                    tmp.sinks.(this_var) = obj.data.avg.(this_var).data.*this_factor;
                end
                try
                    obj.budget.(vars{i}).sources = tmp.sources;
                catch
                    obj.budget.(vars{i}).sources = [];
                end
                try
                    obj.budget.(vars{i}).sinks = tmp.sinks;
                catch
                    obj.budget.(vars{i}).sinks = [];
                end
            end
        end % end method sourcesSinks

        %--------------------------------------------------------------------------------
        function obj = computeNet(obj,vars)
            % ---------------------
            % Computes remainder (net) from budget equation
            % dC/dt = adv + dfz + sms + fg
            % Called in getBudg
            % End result is mmol/m3/s
            %
            % Usage:
            % - obj = computeNet(obj,vars,varargin)
            %
            % Inputs:
            % - vars = budget to close (i.e 'NO2');
            %
            % Example:
            % - obj = computeNet(obj,'NO2');
            % ---------------------
            disp('computeNet: Computing budget remainder (net), hopefully very small');

            % Calculate remainder (net)
            obj.budget.(vars{1}).net = obj.budget.(vars{1}).dcdt - (obj.budget.(vars{1}).adx + ...
                                       obj.budget.(vars{1}).ady  +  obj.budget.(vars{1}).adz + ...
                                       obj.budget.(vars{1}).dfz  +  obj.budget.(vars{1}).sms + ...
                                       obj.budget.(vars{1}).fg   +  obj.budget.(vars{1}).sed);
        end % end method computeNet

        %--------------------------------------------------------------------------------
        function obj = intBudg(obj,vars)
            % ------------------
            % Vertically integrate budget terms
            % Called in getBudg
            % Terms are all in mmol/m3/s, get it in mmol/m2/s by using dz
            %
            % Usage:
            % - obj = intBudg(obj,vars,terms)
            %
            % Inputs:
            % - vars = budget variable (i.e. 'NO2');
            %
            % Example:
            % - obj = intBudg(obj,{'NO2'}); 
            % ------------------
            disp('intBudg: Integrating budget results');
    
            % Grab terms, including remainder
            terms = {'dcdt','adx','ady','adz','dfz','adxy','adv','sms','prod','cons','fg','sed','net'};
            for i = 1:length(terms)
                if ~isempty(obj.budget.(vars{1}).(terms{i}))
                    eval([terms{i},' = obj.budget.(vars{1}).',terms{i},' .* obj.grid.avg.mask_rho3d;']);
                end
            end

            % Integrate vertically (fg term should match real flux)
            % ...mmol/m3/s to mmol/m2/s
            for i = 1:length(terms)
                if ~isempty(obj.budget.(vars{1}).(terms{i}))
                    eval(['obj.budget.(vars{1}).int',terms{i},' = squeeze(nansum(',terms{i},'.*obj.grid.avg.Hz,3));']);
                    obj.budget.(vars{1}).(['int',terms{i}]) = obj.budget.(vars{1}).(['int',terms{i}]) .* obj.grid.mask_rho;
                end
            end

            % ...mmol/m2/s to mmol/s
            for i = 1:length(terms);
                for t = 1:size(obj.budget.(vars{1}).sms,4);
                    eval(['obj.budget.(vars{1}).tot',terms{i},'(t)  = nansum(obj.budget.(vars{1}).int',terms{i},...
                          '(:,:,t) .*obj.grid.area_rho,''all'');']); 
                end
            end
        end % end method intBudg

        %--------------------------------------------------------------------------------
        function plotIntBudg(obj,vars,varargin);
            % ------------------
            % Plot maps of the vertically integrated budget terms
            %
            % Usage:
            % - plotIntBudg(obj,vars,varargin)
            %
            % Inputs (varargin):
            % - time = time record to plot (if length > 1, plot the average)
            % - prc  = percentile to limit colorbar (default = 2)
            %       ...if prc == 2, then caxis = [2nd percentile - 98th percentile]
            %
            % Example:
            % - plotIntBudg(obj,{'N2O'},'time',10,'prc',2)
            % ------------------

            % defaults for optional  arguments
            romsOpt;
            A.time      = [];
            A.prc       = [2];
            A.lonbounds = [floor(min(obj.grid.lon_rho(:))) ceil(max(obj.grid.lon_rho(:)))];
            A.latbounds = [floor(min(obj.grid.lat_rho(:))) ceil(max(obj.grid.lat_rho(:)))];
            A.fontsize  = fontsize;
            A.coast     = coast;
            A.ticks     = ticks;
            A.polygon   = polygon;
            A           = romsObj.parse_pv_pairs(A,varargin); % parse method arguments to A
            
            % Check for integrated budget
            try; obj.budget.(vars{1}).intdcdt;
            catch; obj = intBudg(obj,vars);
            end
            % Process inputs
            if isempty(A.time)
                A.time = 1:size(obj.budget.(vars{1}).intdcdt,3);
            end

            % First generate uniform color bars for each axis    
            % Go through each variables
            terms = {'dcdt','adv','adxy','adx','ady','adz','dfz','sms','fg','sed','net'};
            tits  = {'d$C$/dt','$T_{u,v,w}$','$T_{u,v}$','$T_{u}$','$T_{v}$','T$_{w}$','$D$','$J$','$\Phi$','Sed','Net'};
            for i = 1:length(vars)
                % Get caxis limits that span the range of values
                lims = [];
                for j = 1:length(terms)
                    if ~strcmp(terms{j},'adv')
                        tmp.(terms{j}) = obj.budget.(vars{i}).(['int',terms{j}]);
                    else
                        tmp.adv = [obj.budget.(vars{i}).intadx + ...
                                   obj.budget.(vars{i}).intady + ... 
                                   obj.budget.(vars{i}).intadz];    
                    end
                    % Reduce to time
                    if length(A.time)==1
                        tmp.(terms{j}) = squeeze(tmp.(terms{j})(:,:,A.time));
                    else
                        tmp.(terms{j}) = nanmean(tmp.(terms{j})(:,:,A.time),3);
                    end
                    tmplims(1:2) = romsObj.prclims(tmp.(terms{j}),'prc',A.prc,'bal',1);
                    lims(:,j) = [-(max(abs(tmplims(:)))) (max(abs(tmplims(:))))];
                end
                lims = [-(max(abs(lims(:)))) (max(abs(lims(:))))];
                for j = 1:length(terms);
                    % Plot results
                    dat  = tmp.(terms{j});
                    [fig,cb] = mapPlot(obj,dat,...
                        'levels',linspace(lims(1),lims(end),56),...
                        'caxis',[lims(1) lims(end)],...
                        'cmap',cmocean('balance',55),...
                        'ticks',A.ticks,...
                        'fontsize',A.fontsize,...
                        'coast',A.coast,...
                        'polygon',A.polygon);
                    caxis([lims(1) lims(end)]);
                    hold on
                    title([obj.budget.(vars{i}).info.tits{1},': ',tits{j}],'Interpreter', 'Latex','FontSize',A.fontsize+2);
                    ylabel(cb,obj.budget.(vars{1}).info.units{2},'Interpreter','Latex');
                    fname = [vars{i},'_int',terms{j}];
                    if exist([obj.paths.plots.figs,fname,'.',figsFormat]) == 2
                        cmd = ['rm ',obj.paths.plots.budget,fname,'.',figsFormat];
                        system(cmd);
                    end
                    romsOpt;
                    export_fig(figsFormat,[obj.paths.plots.figs,fname],figsQuality);
                    close(fig); 
                end
            end
        end % end method plotIntBudg

        %--------------------------------------------------------------------------------
        function plotTotBudg(obj,vars,varargin);
            % ------------------
            % Plot the totally integrated budget terms
            %
            % Usage:
            % - plotTotBudg(obj,vars,varargin)
            %
            % Inputs (varargin):
            % - time = time record to plot (if length > 1, plot the average)
            %
            % Example:
            % - plotTotBudg(obj,{'N2O'},'time',1)
            % ------------------

            % defaults for optional  arguments
            A.time      = [1:size(obj.budget.(vars{1}).adx,4)];
            A           = romsObj.parse_pv_pairs(A,varargin); % parse method arguments to A

            % Check for integrated budget
            try; obj.budget.(vars{1}).totdcdt;
            catch; obj = intBudg(obj,vars);
            end

            % Get terms to plot
            terms = {'dcdt','adv','dfz','sms','fg','sed','net'};
            tits  = {'d$C$/dt','$T_{u,v,w}$','$D$','$J$','$\Phi$','Sed','Net'};
            clrs  = colormix(length(terms)+1,'w');
            
            % Loop through vars
            for i = 1:length(vars)

                % Grab mean and standard deviation from all times
                for j = 1:length(terms)
                    if length(A.time)==1
                        tmp.(terms{j}).mean = obj.budget.(vars{i}).(['tot',terms{j}])(A.time);
                        tmp.(terms{j}).std  = 0;
                    else
                        tmp.(terms{j}).mean = nanmean(obj.budget.(vars{i}).(['tot',terms{j}]));
                        tmp.(terms{j}).std  = nanstd(obj.budget.(vars{i}).(['tot',terms{j}]));
                    end
                end

                % Plot total terms
                romsOpt;
                fig = romsObj.piofigs(figtype,1);
                for j = 1:length(terms)
                    y(j)   = [tmp.(terms{j}).mean];
                    ye(j)  = [tmp.(terms{j}).std];
                    % Get bars
                    b(j) = bar(j,y(j)); hold on
                    e(j) = errorbar(j,y(j),ye(j));
                    % Fix
                    b(j).BarWidth  = 1;
                    b(j).FaceColor = clrs(j,:);
                    b(j).FaceAlpha = 0.7;
                    e(j).LineStyle = 'none';
                    e(j).Color     = clrs(j,:);
                    e(j).LineWidth = 1;
                    e(j).CapSize   = 15;
                end
                set(gca,'Color','None');
                xlim([0.5 length(terms)+0.5]);
                set(gca,'XTick',1:length(terms));
                set(gca,'XTickLabel',tits);
                set(gca,'TickLabelInterpreter','Latex');
                title([obj.budget.(vars{i}).info.tits{1},' budget'],'Interpreter', 'Latex');
                ylabel(obj.budget.(vars{i}).info.units{3},'Interpreter','Latex');
                fname = [vars{i},'_budget'];
                if exist([obj.paths.plots.figs,fname,'.',figsFormat]) == 2
                    cmd = ['rm ',obj.paths.plots.figs,fname,'.',figsFormat];
                    system(cmd);
                end
                % Save figure
                export_fig(figsFormat,[obj.paths.plots.figs,fname],figsQuality);
            end
        end % end method plotTotBudg

        %--------------------------------------------------------------------------------
        function [fig,ax] = gridView(obj,varargin)
            % ----------------------
            % Plots the xi/eta indices of a grid file
            % Useful for identifying transect, or lat/lon indices
            %
            % Usage:
            % - [fig,ax] = gridView(obj,varargin)
            %
            % Inputs (varargin):
            % - full     = 1 (default), view entire grid. Use 0 for regional grid
            % - dx       = plot lon lines separated by dx (default = 20)
            % - dy       = plot lat lines separated by dy (default = 20)
            % - ticks    = 0 (no lon/lat labels), 1 (yes), 2 (fancy box)
            % - fontsize = tick font size (default set in romsOpt) 
            % - save     = 1 (print to Figures directory), default = 0 (no print)
            % - figtype  = 'sfig', 'mfig', or 'lfig' (see romsObj.piofigs)
            % - figdim   = (multiplier of width set by figtype)
            %
            % Examples:
            % - [fig,ax] = gridView(obj,'dx',20,'dy',20)  <-- to create figure and axes handles
            % - gridView(obj,'dx',20,'dy',20)             <-- to print to tmpfigs directory
            % ----------------------

            % Grab inputs (varargin)
            romsOpt;
            A.full       = 1;
            A.dx         = [20];
            A.dy         = [20];
            A.ticks      = ticks;
            A.fontsize   = fontsize;
            A.save       = [0];
            A.figtype    = figtype; 
            A.coastcolor = coastcolor;
            A.background = background;
            A.polygon    = polygon;
            A.figdim     = round((obj.grid.ny/obj.grid.nx)*10)/10;
            A            = romsObj.parse_pv_pairs(A,varargin);

            % Get grid
            if A.full == 1
                % Get original lon/lat
                tmp.lon_rho = double(ncread(obj.paths.grid,'lon_rho'));
                tmp.lat_rho = double(ncread(obj.paths.grid,'lat_rho'));
            else
                tmp.lon_rho = obj.grid.lon_rho;
                tmp.lat_rho = obj.grid.lat_rho;
            end

            % Get a blank map 
            warning off % briefly, for contourf warning
            dat = zeros(size(tmp.lon_rho));
            [fig,cb] = mapPlot(obj,dat,...
                'ticks',A.ticks,...
                'font',A.fontsize,...
                'figtype',A.figtype,...
                'figdim',A.figdim,...
                'background',A.background,...
                'coastcolor',A.coastcolor,...
                'polygon',A.polygon);
            warning on
            delete(cb);

            % Add xi/eta grid lines
            set(0,'CurrentFigure',fig);
            [a,b]  = size(tmp.lon_rho);
            for i = 1:(A.dx*2):a
                m_plot(tmp.lon_rho(i,:),tmp.lat_rho(i,:),'r');
                for j = 1:(A.dy*2):b
                    hold on
                    m_plot(tmp.lon_rho(:,j),tmp.lat_rho(:,j),'b');
                end
            end
            for i = A.dx:(A.dx*2):a
                m_plot(tmp.lon_rho(i,:),tmp.lat_rho(i,:),'r');
                for j = A.dy:(A.dy*2):b
                    hold on
                    m_plot(tmp.lon_rho(:,j),tmp.lat_rho(:,j),'b');
                end
            end

            % Add labels to indices
            for i = 1:(A.dx*2):a
                m_text(tmp.lon_rho(i,end),tmp.lat_rho(i,end),num2str(i),...
                    'fontsize',A.fontsize);
                for j = 1:(A.dy*2):b
                    hold on
                    m_text(tmp.lon_rho(end,j),tmp.lat_rho(end,j),num2str(j),...
                        'fontsize',A.fontsize);
                end
            end
            for i = A.dx:(A.dx*2):a
                m_text(tmp.lon_rho(i,1),tmp.lat_rho(i,1),num2str(i),...
                    'fontsize',A.fontsize);
                for j = A.dy:(A.dy*2):b
                    hold on
                    m_text(tmp.lon_rho(1,j),tmp.lat_rho(1,j),num2str(j),...
                        'fontsize',A.fontsize);
                end
            end

            % Save figure
            hold on
            fname = [obj.info.simName,'_grid'];
            if A.save == 1
                romsOpt
                export_fig(figsFormat,[obj.paths.plots.figs,fname]);
            end
            if nargout < 1
                romsObj.pltjpg(1);
                close(fig(1));
            end
        end % end method gridView
        
        %--------------------------------------------------------------------------------
        function [fig,ax] = regionView(obj,varargin)
            % ----------------------
            % Plots a ROMS grid map along with a regional grid
            %
            % Usage:
            % - [fig,ax] = regionView(obj,varargin)
            %
            % Inputs (varargin):
            % - ticks = 0 (no lon/lat labels), 1 (yes), 2 (fancy box)
            % - font  = tick font size (default = 12) 
            % - save  = 1 (print), 0 (no print)
            %
            % Example:
            % - [fig,ax] = regionView(obj);
            % ----------------------

            % Grab inputs (varargin)
            romsOpt;
            A.ticks     = ticks;
            A.fontsize  = fontsize;
            A.save      = [0];
            A           = romsObj.parse_pv_pairs(A,varargin);

            % Get original lon/lat
            tmp.lon_rho = ncread(obj.paths.grid,'lon_rho');
            tmp.lat_rho = ncread(obj.paths.grid,'lat_rho');
    
            % Map of region
            if ~isempty(obj.grid.lon_rho)
                % Generate whole map
                fig(1)  = romsObj.piofigs('mfig',1.5);
                set(0,'CurrentFigure',fig(1));
                [ax(1)] = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks,'font',A.fontsize);
                hold on
                % Plot grid box
                m_plot(tmp.lon_rho(1,:),  tmp.lat_rho(1,:),'k','linewidth',2);
                m_plot(tmp.lon_rho(:,1),  tmp.lat_rho(:,1),'k','linewidth',2);
                m_plot(tmp.lon_rho(end,:),tmp.lat_rho(end,:),'k','linewidth',2);
                m_plot(tmp.lon_rho(:,end),tmp.lat_rho(:,end),'k','linewidth',2);
                % Plot region box
                m_plot(obj.grid.lon_rho(1,:),obj.grid.lat_rho(1,:),'--k','linewidth',2);
                m_plot(obj.grid.lon_rho(:,1),obj.grid.lat_rho(:,1),'--k','linewidth',2);
                m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'--k','linewidth',2);
                m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'--k','linewidth',2);
                fname = [obj.info.simName,'_region'];
                if A.save == 1
                    romsOpt;
                    export_fig(figsFormat,[obj.paths.plots.figs,fname]);
                else
                    pp = 1;
                    romsObj.pltjpg(1);
                end
                if nargout < 1 & pp ~= 1;
                    romsObj.pltjpg(1);
                end
            end
        end % end method regionView

        %--------------------------------------------------------------------------------
        function obj = Dist2Coast(obj)
            % --------------------
            % This function is automatically called if 'coast' is specified during initROMS
            % Calculate each grid cell's distance-to-coast (in meters) using 'mask_rho' as a proxy for the coastline
            % Adapted from Pierre Damien's routines
            %
            % Usage:
            % - obj = Dist2Coast(obj) 
            % 
            % Output: 
            % - obj.grid.coastdist
            % --------------------
            disp('Dist2Coast: Calculating distance from coast');

            % Load grid?
            try; obj.grid.lon_rho;
            catch; obj = loadGrid(obj);
            end

            % Get grid info
            lon  = obj.grid.lon_rho;
            lat  = obj.grid.lat_rho;
            mask = obj.grid.mask_rho;
            [Mp, Lp] = size(obj.grid.lon_rho);

            % Conversion to radians
            d2r = pi/180;
            cdist = 0*lon + 1e10;
            ic  = 1; jc  = 1;
            ncx = 1; ncy = 1;
            i0 = 1 + (ic-1)*ceil(Lp/ncx);
            i1 = i0+ceil(Lp/ncx) + 20;
            j0 = 1 + (jc-1)*ceil(Mp/ncy);
            j1 = j0+ceil(Mp/ncy) + 20;
            i1 = min(i1,Lp); j1 = min(j1,Mp);
            lons = lon(j0:j1,i0:i1);
            lats = lat(j0:j1,i0:i1);
            masks= mask(j0:j1,i0:i1);
            lab = 0*masks + 1; 
            lab(2:end-1,2:end-1) = masks(1:end-2,2:end-1)+...
                                   masks(3:end,2:end-1)+...
                                   masks(2:end-1,1:end-2)+...
                                   masks(2:end-1,3:end);
            mlon = lons;mlon(masks>0|lab<1) = [];
            mlat = lats;mlat(masks>0|lab<1) = [];
            if mlon
                for j = j0:j1
                    [j j1];
                    for i = i0:i1
                        if mask(j,i) < 1
                            cdist(j,i) = 0;
                        else
                            dist = romsObj.gc_dist(lon(j,i)*d2r,lat(j,i)*d2r,mlon*d2r,mlat*d2r);
                            mdist = min(min(dist));
                            cdist(j,i) = min(mdist,cdist(j,i));
                        end
                    end
                end
            end

            % Save output
            obj.grid.coastdist = cdist;
        end % end method Dist2Coast

        %--------------------------------------------------------------------------------
        function [fig,cb] = mapPlot(obj,dat,varargin);
            % -----------------------
            % A way to quickly plot a 2D field
            %
            % Usage:
            % - [fig,cb] = mapPlot(obj,dat,varargin);
            %
            % Inputs:
            % - dat = 2D field to plot (prepare it before using the script)
            %
            % Varargin:
            % - meta          = option to include structure to name and units (e.g., obj.data.avg.temp)
            % - lonbounds     = x-boundaries (defaults to whole domain)
            % - latbounds     = y-boundaries (defaults to whole domain)
            % - ticks         = 2 = fancy, 1 = on, 0 = off (default set in romsOpt)
            % - background    = background color (default set in romsOpt)
            % - coastcolor    = coast color (default set in romsOpt);
            % - fontsize      = fontsize (default set in romsOpt)
            % - figtype       = 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm) (default set in romsOpt)
            % - figdim        = figtype multiplier for height (e.g., 1 sets same width and height) 
            % - levels        = hard-coded levels to plot (can't be used with A.prc)
            % - cmap          = colormap(default = thermal for no balance, balance for balance)
            % - caxis         = colorbar limits 
            % - prc           = percentage to limit colorbar axes (if no levels supplied)
            % - bal           = force a balanced colorbar around 0 (1 == yes, 0 == no)
            % - log           = log-scale (1), use with caxis to set limits
            % - XaxisLocation = Override x-axis ticklabels location (default = bottom) 
            % - YaxisLocation = Override y-axis ticklabels location (default = left) 
            % - coast         = m_map map quality (coast, crude, low, high, intermediate)
            % - polygon       = (0) to turn off boundary polygon (default set in romsOpt)
            % -----------------------
            
            % User-inputs
            romsOpt
            A.meta          = [];
            A.lonbounds     = [];
            A.latbounds     = [];
            A.lonticks      = [];
            A.latticks      = [];
            A.ticks         = ticks;
            A.background    = background;
            A.coastcolor    = coastcolor;
            A.fontsize      = fontsize;
            A.figtype       = figtype;
            A.polygon       = polygon;
            A.coast         = coast;
            A.figdim        = round((obj.grid.ny/obj.grid.nx)*10)/10;;
            A.prc           = 2;
            A.bal           = 0;
            A.levels        = [];
            A.cmap          = [];
            A.caxis         = [];
            A.log           = 0;
            A.XaxisLocation = 'bottom';
            A.YaxisLocation = 'left';
            A = romsObj.parse_pv_pairs(A,varargin);

            % Make double
            if size(dat)==size(obj.grid.lon_rho);
                lon = double(obj.grid.lon_rho);
                lat = double(obj.grid.lat_rho);
            elseif size(dat)==size(obj.grid.lon_u);
                lon = double(obj.grid.lon_u);
                lat = double(obj.grid.lat_u);
            elseif size(dat)==size(obj.grid.lon_v)
                lon = double(obj.grid.lon_v);
                lat = double(obj.grid.lat_v);
            else
                disp('    ERROR(mapPlot): Check dimensions of input');
                kill
            end

            % Get auto-bounds if empty
            if isempty(A.lonbounds) & isempty(A.latbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            elseif isempty(A.lonbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = A.latbounds;
            elseif isempty(A.latbounds);
                lonbounds = A.lonbounds;
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            else
                lonbounds = A.lonbounds;
                latbounds = A.latbounds;
            end

            % Set up ticks
            if isempty(A.lonticks)
                dx = round(diff(lonbounds)/60)*10;
                lonticks = (lonbounds(1):dx:lonbounds(2));    
                lonticks = round(lonticks);
            else
                lonticks = A.lonticks;
            end
            if isempty(A.latticks)
                dy = round(diff(latbounds)/60)*10;
                latticks = (latbounds(1):dy:latbounds(2));
                latticks = round(latticks);
            else
                latticks = A.latticks;
            end

            % Initiate figure
            fig = romsObj.piofigs(A.figtype,A.figdim);

            % Get colormap limits
            if isempty(A.levels)
                clims = romsObj.prclims(dat,'prc',A.prc,'bal',A.bal);
                clevs = linspace(clims(1),clims(2),31); 
            else
                clevs = A.levels;
                clims = [A.levels(1) A.levels(end)];
            end
            if ~isempty(A.caxis);
                clevs(1) = A.caxis(1);
                clevs(end) = A.caxis(2);
                clims = A.caxis;
            end
            if max(dat(:)) == 0 & min(dat(:)) == 0 | isnan(max(dat)) ==1 & isnan(min(dat(:)) == 1);
                dat = nan(size(dat));
                clevs = [0 1];
                clims = linspace(0,1,11);
            else
                dat(dat<clevs(1))   = clevs(1);
                dat(dat>clevs(end)) = clevs(end);
            end

            % Get colormap
            if ischar(A.cmap)
                A.cmap = cmocean(A.cmap,length(clevs)-1);
            elseif isempty(A.cmap)
                if min(clims(:))<0 & max(clims(:)) > 0 | A.bal == 1
                    A.cmap = cmocean('balance',length(clevs)-1);
                else
                    A.cmap = cmocean('thermal',length(clevs)-1);
                end
            end

            % Make map
            set(0,'CurrentFigure',fig);
            m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
            hold on
            if A.ticks == 0
                m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
                       'xticklabels',[],'yticklabels',[],'backgroundcolor',A.background); drawnow
            elseif A.ticks == 1
                m_grid('box','on','linestyle','none','xtick',lonticks,...
                       'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
            elseif A.ticks == 2
                m_grid('box','fancy','linestyle','none','xtick',lonticks,...
                       'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
            end
            hold on
            m_contourf(lon,lat,dat,clevs,'LineStyle','none');
            cb = colorbar; drawnow
            cb.FontSize = A.fontsize;
            try
                caxis([clims]);
            catch
                caxis([0 1]);
            end
            ax = get(gca);
            if strcmp(A.coast,'coast');
                m_coast('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'crude')
                m_gshhs_c('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'low')
                m_gshhs_l('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'high')
                m_gshhs_h('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'intermediate')
                m_gshhs_i('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'full')
                m_gshhs_f('patch',A.coastcolor,'edgecolor','k'); drawnow
            end
            colormap(gca,A.cmap);
            if A.log == 1 & ~isempty(A.caxis)
                set(gca,'ColorScale','log');
                caxis(A.caxis);
            end

            % Include meta data if available
            if ~isempty(A.meta)
                title(A.meta.name,'Interpreter','Latex','fontsize',A.fontsize);
                ylabel(cb,A.meta.units,'Interpreter','Latex','fontsize',A.fontsize);
            end
        
            % Show boundary polygon
            if A.polygon
                hold on
                m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'-k','linewidth',1);
            end
            
            % Print figure
            if nargout<1
                romsObj.pltjpg(1);
            end
        end % end method mapPlot

        %--------------------------------------------------------------------------------
        function [figs,cbs] = mapCmp(obj,dat1,dat2,varargin);
            % -----------------------
            % A way to quickly plot a 2D field comparison
            % Produces 3 plots: dat1 and dat2 fields with the same colorbar, 
            % and a 3rd plot of the difference between them (dat1 - dat2)
            %
            % Usage:
            % - [figs,cbs] = mapCmp(obj,dat1,dat2,varargin);
            %
            % Inputs:
            % - dat1 = 2D field to plot (prepare it before using the script)
            % - dat2 = same same but different 2D field to plot (prepare it before using the script)
            %
            % Varargin:
            % - lonbounds:  x-boundaries (defaults to whole domain)
            % - latbounds:  y-boundaries (defaults to whole domain)
            % - ticks:      2 = fancy, 1 = on, 0 = off
            % - background: background color (default 'LightGray');
            % - coastcolor: coast color (default 'DimGray');
            % - fontsize:   default 10
            % - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
            % - figdim:     default 1 (same height as width, its a multiplier)
            % - prc:        percentage to limit colorbar axes
            % - bal:        balance colorbar around 0 (1 == yes, 0 == no, 2 == set min to 0)
            % - levels:     hard-coded levels to plot (overrides A.prc, A.bal)
            % - difflevels: hard-coded difference levels to plot
            % - cmap:       colormap(default = thermal for no balance, balance for balanced colormap)
            % - dmap:       colormap for difference plot (default = balance)
            % - units:      string containing units for colorbar
            % ----------------------

            % User-inputs
            romsOpt;
            A.lonbounds  = [];
            A.latbounds  = [];
            A.ticks      = ticks;
            A.background = background;
            A.coastcolor = coastcolor;
            A.fontsize   = fontsize;
            A.figtype    = figtype;
            A.figdim     = round((obj.grid.ny/obj.grid.nx)*10)/10;
            A.prc        = 0.1;
            A.bal        = 0;
            A.levels     = [];
            A.difflevels = [];
            A.cmap       = [];
            A.dmap       = [];
            A.units      = [];
            A = romsObj.parse_pv_pairs(A,varargin);
            
            % Get automatic levels?
            if isempty(A.levels)
                all_dat = [dat1(:) dat2(:)];
                A.levels = romsMaster.prclims(all_dat,'prc',A.prc,'bal',A.bal);
                A.levels = linspace(A.levels(1),A.levels(2),20);
            end
            % Get differences
            diff_dat  = dat1 - dat2;
            if isempty(A.difflevels)
                A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1); 
                A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
            end

            % If string used for A.cmap or A.dmap, call cmocean
            % If empty, use defaults
            if ischar(A.cmap)
                A.cmap = cmocean(A.cmap,length(A.levels)-1);
            elseif isempty(A.cmap)
                if min(all_dat) < 0 & max(all_dat) > 0 | A.bal == 1
                    A.cmap = cmocean('balance',length(A.levels)-1);
                else
                    A.cmap = cmocean('thermal',length(A.levels)-1);
                end
            end
            if ischar(A.dmap)
                A.dmap = cmocean(A.dmap,length(A.difflevels)-1);
            elseif isempty(A.dmap)
                A.dmap = cmocean('balance',length(A.difflevels)-1);
            end

            % Make figs(1) and figs(2)
            [figs(1),cbs(1)] = mapPlot(obj,dat1,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',A.cmap);
            [figs(2),cbs(2)] = mapPlot(obj,dat2,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',A.cmap);
            
            % Make figs(3)
            [figs(3),cbs(3)] = mapPlot(obj,diff_dat,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels,'cmap',A.dmap);

            % Add units?
            if ~isempty(A.units)
                ylabel(cbs(1),A.units,'Interpreter','Latex','FontSize',A.fontsize);
                ylabel(cbs(2),A.units,'Interpreter','Latex','FontSize',A.fontsize);
                ylabel(cbs(3),A.units,'Interpreter','Latex','FontSize',A.fontsize);
            end

            % Auto-print if no output
            if nargout == 0
                set(0,'CurrentFigure',figs(1));
                romsObj.pltjpg(1);
                
                set(0,'CurrentFigure',figs(2));
                romsObj.pltjpg(2);

                set(0,'CurrentFigure',figs(3));
                romsObj.pltjpg(3);
            end
        end % end method mapCmp

        %--------------------------------------------------------------------------------
        function [fig,ax,cb] = axPlot(obj,fig,ax,dat,varargin);
            % -----------------------
            % A way to quickly plot a 2D field into a specific axis (i.e. subplot)
            %
            % Usage:
            % - [fig,ax,cb] = axPlot(obj,fig,ax,dat,varargin);
            %
            % Inputs:
            % - fig = figure handle
            % - ax  = axes handle
            % - dat = 2D field to plot (prepare it before using the script)
            %
            % Varargin:
            % - lonbounds:  x-boundaries (defaults to whole domain)
            % - latbounds:  y-boundaries (defaults to whole domain)
            % - ticks:      2 = fancy, 1 = on, 0 = off
            % - xtick/ytick 0 = off (to toggle only one tick on/off)
            % - background: background color (default 'LightGray');
            % - coastcolor: coast color (default 'DimGray');
            % - fontsize:   default 10
            % - prc:        percentage to limit colorbar axes
            % - bal:        balance colorbar around 0 (1 == yes, 0 == no)
            % - levels:     hard-coded levels to plot (can't be used with A.prc)
            % - cmap:       colormap(default = thermal for no balance, balance for balance)
            % - caxis       colorbar limits 
            % - log         log-scale (1), use with caxis to set limits
            % ---------------------
            
            % User-inputs
            romsOpt;
            A.lonbounds  = [];
            A.latbounds  = [];
            A.lonticks   = [];
            A.latticks   = [];
            A.ticks      = ticks;
            A.ytick      = 1;
            A.xtick      = 1;
            A.background = background;
            A.coastcolor = coastcolor;
            A.fontsize   = fontsize;
            A.prc        = 2;
            A.bal        = 0;
            A.levels     = [];
            A.cmap       = cmocean('thermal');
            A.caxis      = [];
            A.log        = 0;
            A = romsObj.parse_pv_pairs(A,varargin);

            % Balance override
            if A.bal == 1
                A.cmap = cmocean('balance');
            end

            % Make double
            lon = double(obj.grid.lon_rho);
            lat = double(obj.grid.lat_rho);

            % Get auto-bounds if empty
            if isempty(A.lonbounds) & isempty(A.latbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            elseif isempty(A.lonbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = A.latbounds;
            elseif isempty(A.latbounds);
                lonbounds = A.lonbounds;
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            else
                lonbounds = A.lonbounds;
                latbounds = A.latbounds;
            end

            % Set up ticks
            if abs(diff(lonbounds)) > 50
                dx = 10;
            else
                dx = 5;
            end
            if abs(diff(latbounds)) > 50
                dy = 10;
            else
                dy = 5;
            end
            if isempty(A.latticks);
                latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
            else
                latticks = A.latticks;
            end
            if isempty(A.lonticks);
                lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));
            else    
                lonticks = A.lonticks;
            end

            % Set current axes
            set(0,'CurrentFigure',fig);
            set(fig,'CurrentAxes',ax);

            % Get colormap limits
            if isempty(A.levels)
                clims = romsObj.prclims(dat,'prc',A.prc,'bal',A.bal);
                clevs = linspace(clims(1),clims(2),31); 
            else
                clevs = A.levels;
                clims = [A.levels(1) A.levels(end)];
            end
            if ~isempty(A.caxis);
                clevs(1) = A.caxis(1);
                clevs(end) = A.caxis(2);
                clims = A.caxis;
            end
            if max(dat(:)) == 0 & min(dat(:)) == 0 | isnan(max(dat)) ==1 & isnan(min(dat(:)) == 1);
                dat = nan(size(dat));
                clevs = [0 1];
                clims = linspace(0,1,11);
            else
                dat(dat<clevs(1))   = clevs(1);
                dat(dat>clevs(end)) = clevs(end);
            end

            % Make map
            m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
            hold on
            if A.ticks == 0
                m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
                       'xticklabels',[],'yticklabels',[],'backgroundcolor',A.background); drawnow
            elseif A.ticks == 1
                if A.ytick == 1 & A.xtick == 1
                    m_grid('box','on','linestyle','none','xtick',lonticks,...
                           'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
                elseif A.ytick ==1 & A.xtick == 0
                    m_grid('box','on','linestyle','none','xtick',0,...
                           'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',[]);
                elseif A.ytick == 0 & A.xtick == 1
                    m_grid('box','on','linestyle','none','xtick',lonticks,...
                           'ytick',0,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',[],'xticklabels',lonticks);
                end
            elseif A.ticks == 2
                m_grid('box','fancy','linestyle','none','xtick',lonticks,...
                       'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
            end
            hold on
            try
                m_contourf(lon,lat,dat,clevs,'LineStyle','none');
            catch
            end
            cb = colorbar; drawnow
            cb.FontSize = A.fontsize;
            try
                caxis([clims])
            catch
                caxis([0 1]);
            end
            m_coast('patch',A.coastcolor,'edgecolor','k'); drawnow
            colormap(gca,A.cmap);
            if A.log == 1 & ~isempty(A.caxis)
                set(gca,'ColorScale','log');
                caxis(A.caxis);
            end
        end % end method axPlot

        %--------------------------------------------------------------------------------
        function [fig] = quickMap(obj,varargin);
            % -----------------------
            % A way to quickly plot a ROMS map
            %
            % Usage:
            % - [fig] = quickMap(obj,varargin);
            %
            % Varargin:
            % - lonbounds:  x-boundaries (defaults to whole domain)
            % - latbounds:  y-boundaries (defaults to whole domain)
            % - ticks:      2 = fancy, 1 = on, 0 = off
            % - background: background color (default 'LightGray');
            % - coastcolor: coast color (default 'DimGray');
            % - fontsize:   default 10
            % - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
            % - figdim:     default 1 (same height as width, its a multiplier)
            % - poly:       default 1 (plot ROMS boundaries)
            % ---------------------

            % Reject if no output provided
            if nargout < 1
                disp('    ERROR(quickMap): No point using this without [fig] output, see help quickMap');
                return
            end

            % User-inputs
            romsOpt;
            A.lonbounds  = [obj.grid.minlon_rho obj.grid.maxlon_rho];
            A.latbounds  = [obj.grid.minlat_rho obj.grid.maxlat_rho];
            A.latticks   = [];
            A.lonticks   = [];
            A.ticks      = ticks;
            A.box        = 'on';
            A.background = background;
            A.coastcolor = coastcolor;
            A.fontsize   = fontsize;
            A.figtype    = figtype;
            A.figdim     = round((obj.grid.ny/obj.grid.nx)*10)/10;;
            A.coast      = coast;
            A.poly       = 1;
            A = romsObj.parse_pv_pairs(A,varargin);

            % Make double
            lon = double(obj.grid.lon_rho);
            lat = double(obj.grid.lat_rho);

            % Get auto-bounds if empty
            if isempty(A.lonbounds) & isempty(A.latbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            elseif isempty(A.lonbounds);
                lonbounds = [floor(min(lon(:))) ceil(max(lon(:)))];
                latbounds = A.latbounds;
            elseif isempty(A.latbounds);
                lonbounds = A.lonbounds;
                latbounds = [floor(min(lat(:))) ceil(max(lat(:)))];
            else
                lonbounds = A.lonbounds;
                latbounds = A.latbounds;
            end

            % Set up ticks
            if isempty(A.lonticks)
                dx = round(diff(lonbounds)/60)*10;
                lonticks = (lonbounds(1):dx:lonbounds(2));    
                lonticks = round(lonticks);
            else
                lonticks = A.lonticks;
            end
            if isempty(A.latticks)
                dy = round(diff(latbounds)/60)*10;
                latticks = (latbounds(1):dy:latbounds(2));
                latticks = round(latticks);
            else
                latticks = A.latticks;
            end

            % Make map
            fig = romsObj.piofigs(A.figtype,A.figdim);
            set(0,'CurrentFigure',fig);
            m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
            hold on
            if A.ticks == 0
                m_grid('box',A.box,'linestyle','none','xtick',0,'ytick',0,...
                       'xticklabels',[],'yticklabels',[],'backgroundcolor',A.background); drawnow
            elseif A.ticks == 1
                m_grid('box','on','linestyle','none','xtick',lonticks,...
                       'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
            elseif A.ticks == 2
                m_grid('box','fancy','linestyle','none','xtick',lonticks,...
                       'ytick',latticks,'backgroundcolor',A.background,'fontsize',A.fontsize,'yticklabels',latticks,'xticklabels',lonticks);
            end
            if strcmp(A.coast,'coast');
                m_coast('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'crude')
                m_gshhs_c('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'low')
                m_gshhs_l('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'high')
                m_gshhs_h('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'intermediate')
                m_gshhs_i('patch',A.coastcolor,'edgecolor','k'); drawnow
            elseif strcmp(A.coast,'full')
                m_gshhs_f('patch',A.coastcolor,'edgecolor','k'); drawnow
            end
    
            if A.poly == 1
                m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'k','linewidth',1);
            end
        end % end method quickMap

        %--------------------------------------------------------------------------------
        function [fig,cb] = slicePlot(obj,slicedata,t,varargin)
            % ------------------
            % Plots 2D sliced variables obtained from sliceROMS or sliceDiag
            %
            % Usage:
            % - [fig,cb] = slicePlot(obj,slicedata,t,varargin);
            %
            % Inputs:
            % - slicedata = sliced data from sliceROMS or sliceDiag (must be size(obj.slice.deg)) 
            % - t         = time to plot (0 == average)
            %
            % Optional inputs (varargin):
            % - xlims:      hard-coded x-limits (degrees)
            % - zlims:      hard-coded z-limits (meters)
            % - cmap:       colormaps (if used, must be a cell array of length (vars))
            % - fontsize:   default 10
            % - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
            % - figdim:     default 1 (same height as width, its a multiplier)
            % - prc:        percentage to limit colorbar axes
            % - bal:        balance colorbar around 0 (1 == yes, 0 == no)
            % - levels:     hard-coded levels to plot (can't be used with A.prc)
            % - background: background color (i.e. coast). Default = rgb('DimGray');
            %
            % Example:
            % - [fig,cb] = slicePlot(obj,tmpdata,0,'xlims',[140 260],'zlims',[-500 0]);
            % ------------------
            
            % User-inputs
            romsOpt;
            A.xlims      = [];
            A.zlims      = [];
            A.cmap       = [];
            A.fontsize   = fontsize;
            A.figtype    = figtype;
            A.figdim     = 0.33;
            A.prc        = 0.5;
            A.bal        = 0;
            A.levels     = [];
            A.background = coastcolor;
            A = romsObj.parse_pv_pairs(A,varargin);

            % Reduce data to average, or specified slice    
            if t == 0
                tmpdeg = nanmean(obj.slice.deg,3);
                tmpdep = nanmean(obj.slice.depth,3);
                tmpdat = nanmean(slicedata,3);
            else
                tmpdeg = squeeze(obj.slice.deg(:,:,t));
                tmpdep = squeeze(obj.slice.depth(:,:,t));
                tmpdat = squeeze(slicedata(:,:,t));
            end

            % Check for xlim or zlims
            if ~isempty(A.xlims);
                xind = A.xlims(1) <= tmpdeg & tmpdeg <= A.xlims(2);
            end
            if ~isempty(A.zlims);
                zind = A.zlims(1) <= tmpdep & tmpdep <= A.zlims(2);
            end

            % Get universal levels
            if isempty(A.levels)
                % Reduce data?
                if ~isempty(A.xlims);
                    tmpdat(xind==0) = NaN;
                end
                if ~isempty(A.zlims);
                    tmpdat(zind==0) = NaN;
                end
                A.levels = romsObj.prclims(tmpdat(:),'prc',A.prc,'bal',A.bal);
                A.levels = linspace(A.levels(1),A.levels(2),20);
            end
            tmpdat(tmpdat<A.levels(1)) = A.levels(1);
            tmpdat(tmpdat>A.levels(end)) = A.levels(end);
            
            % Generate figure(s)    
            fig = romsObj.piofigs(A.figtype,A.figdim);
            set(0,'CurrentFigure',fig);
            if strcmp(obj.slice.coord,'latitude') | strcmp(obj.slice.coord,'longitude');
                contourf(tmpdeg,tmpdep,tmpdat,A.levels,'linestyle','none');
            end
            if ~isempty(A.xlims);
                xlim(A.xlims);
            end
            if ~isempty(A.zlims);
                ylim(A.zlims);
            end
            if strcmp(obj.slice.coord,'longitude')
                xlabel('Latitude','Interpreter','Latex','FontSize',A.fontsize);
                xlbl = get(gca,'XTickLabel');
                for i = 1:length(xlbl)
                    if str2num(xlbl{i})<0
                        newlbl{i} = [num2str(-str2num(xlbl{i})),char(176),'S'];
                    else
                        newlbl{i} = [num2str(str2num(xlbl{i})),char(176),'N'];
                    end
                end
                set(gca,'XTickLabel',newlbl);    
            elseif strcmp(obj.slice.coord,'latitude')    
                xlabel('Longitude','Interpreter','Latex','FontSize',A.fontsize);
                xlbl = get(gca,'XTickLabel');
                for i = 1:length(xlbl)
                    if str2num(xlbl{i})>180
                        newlbl{i} = [num2str(str2num(xlbl{i})-360),char(176),'W'];
                    else
                        newlbl{i} = [num2str(str2num(xlbl{i})),char(176),'E'];
                    end
                end
                set(gca,'XTickLabel',newlbl);    
            end
            ylabel('Depth (m)','Interpreter','Latex','FontSize',A.fontsize);
            cb = colorbar('location','eastoutside');
            caxis([A.levels(1) A.levels(end)]);
            if ~isempty(A.cmap)
                if ischar(A.cmap)
                    set(gca,'Colormap',cmocean(A.cmap,length(A.levels)));
                else
                    set(gca,'Colormap',A.cmap);
                end
            end
            set(gca,'FontSize',A.fontsize);
            if nanmean(tmpdep(:)) > 0
                set(gca,'YDir','Reverse');
            end
            set(gca,'Color',A.background);
            set(gcf,'inverthardcopy','off');

            % Print if no output provided
            if nargout<1
                romsObj.pltjpg(1);
            end
        end % end method slicePlot

        %--------------------------------------------------------------------------------
        function obj = sliceDiag(obj,vars,choice,deg,varargin);
            % -------------------
            % Takes 2D depth slice of diagnostic data along a given latitude or longitude
            % Options are set by user
            % 
            % Usage:
            % - obj = sliceDiag(obj,vars,choice,deg,varargin);
            % 
            % Inputs:
            % - vars   = diagnostic variable(s) to slice, as a cell array
            %             (see obj.paths.diag)
            % - choice = 'lon','lat','xi', or 'eta' 
            %              (lat/lon slices along a given lat/lon degree)
            %              (x/y slices along a given x or y index (lon_rho or lat_rho))
            % - deg    = lon/lat degree or x/y index
            %
            % Inputs (varargin):
            % - zlim = depth limit to restrict output
            % 
            % Example
            % - obj = sliceDiag(obj,{'temp','salt'},'lon',0);
            % 
            % -------------------
            disp('sliceDiag: Get transects of validation data');

            % Get optional inputs
            A.zlim = inf;
            A = romsObj.parse_pv_pairs(A,varargin);

            % Grab longitude/latitude data
            lon  = obj.grid.lon_rho;
            lat  = obj.grid.lat_rho;

            % Choose latitude, longitude
            if strcmp(choice,'lat')
                lonlat = lat;
                dmsn   = obj.grid.nx; 
                nz     = length(obj.grid.z_avg_dep(obj.grid.z_avg_dep<=A.zlim)); 
                nt     = 12; 
                dstr   = ['latitude'];
            elseif strcmp(choice,'lon')
                lonlat = lon; 
                dmsn   = obj.grid.ny; 
                nz     = length(obj.grid.z_avg_dep(obj.grid.z_avg_dep<=A.zlim)); 
                nt     = 12; 
                dstr   = ['longitude'];
            end
            
            % Get lon/lat data (outside parfor)
            tmpdeg  = NaN(dmsn,1);
            for i = 1:dmsn
                if dmsn == obj.grid.nx;
                    tmpdeg(i)  = interp1(squeeze(lat(i,:)),squeeze(lon(i,:)),deg);
                elseif dmsn == obj.grid.ny;
                    tmpdeg(i)  = interp1(squeeze(lon(:,i)),squeeze(lat(:,i)),deg);
                end
            end
            tmpdeg = repmat(tmpdeg,1,1,nz); tmpdeg = permute(tmpdeg,[1 3 2]);

            % Dimensions to fill
            fillmat  = [dmsn nz];
            tmpdepth = obj.grid.z_avg_dep(obj.grid.z_avg_dep<=A.zlim)'.*ones(fillmat);

            % Perform slices of each variable
            for ff = 1:length(vars);
                obj.diag.(vars{ff}).slice = [];
                % get path and coords for current variable
                if ~isfield(obj.paths.diag,(vars{ff}));
                    disp(['    ERROR(sliceDiag): ',vars{ff},' is not a diagnostic variable']);
                    disp(['    Filling with NaN']);
                    obj.diag.(vars{ff}).slice = nan(dmsn,nz,nt);
                    obj.diag.(vars{ff}).name = 'null';
                    obj.diag.(vars{ff}).units = 'null';
                    continue
                end
                curVar  = obj.paths.diag.(vars{ff});
                for i = 1:length(curVar.file); 
                    if strcmp(curVar.type{i},'nc');
                        tmp.lon   = ncread(curVar.file{i},curVar.lon{i}); tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
                        tmp.lat   = ncread(curVar.file{i},curVar.lat{i});
                        tmp.depth = ncread(curVar.file{i},curVar.zvar{i});
                        if nanmean(tmp.depth(:)) < 0
                            tmp.depth = -tmp.depth;
                        end
                        tmp.data  = ncread(curVar.file{i},curVar.var{i});
                        % Apply any factors (eg nMol to uMol)    
                        if isfield(curVar,'factor');
                            tmp.data = tmp.data*curVar.factor{i};
                        end
                    elseif strcmp(curVar.type{i},'mat');
                        tmp.lon              = load(curVar.file{i},curVar.lon{i});
                        tmp.lat              = load(curVar.file{i},curVar.lat{i});
                        tmp.depth            = load(curVar.file{i},curVar.zvar{i});
                        tmp.data           = load(curVar.file{i},curVar.var{i});
                        tmp.lon               = tmp.lon.(curVar.lon{i});
                        tmp.lat            = tmp.lat.(curVar.lat{i});
                        tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
                        tmp.depth          = tmp.depth.(curVar.zvar{i});
                        tmp.data           = tmp.data.(curVar.var{i});
                        % Apply any factors (eg nMol to uMol)    
                        if isfield(curVar,'factor');
                            tmp.data = tmp.data*curVar.factor{i};
                        end
                    end

                    % make lon/lat data gridded if it isnt setup that way
                    if sum(size(tmp.lon(:,:))==1) > 0
                        [tmp.latg,tmp.long,tmp.depthg] = meshgrid(tmp.lat,tmp.lon,tmp.depth);
                    else
                        tmp.latg = tmp.lat;    
                        tmp.long = tmp.lon;
                        tmp.depthg = tmp.depth;
                        tmp.lat = squeeze(tmp.lat(1,:,1,1));
                        tmp.lon = squeeze(tmp.lon(:,1,1,1));
                        tmp.depth = squeeze(tmp.depth(1,1,:,1));
                    end
                    
                    % go through each time record
                    for rcrd = 1:nt
                        % clear variables to fill
                        tmp.latr   = []; 
                        tmp.lonr   = [];
                        tmp.depthr = [];
                        tmp.datar  = [];
                        % Get reduced matrix to feed into scatteredInterpolant
                        if strcmp(choice,'lat');
                            idx = find(abs([tmp.lat-deg]) == min(abs([tmp.lat-deg])));
                            tmp.latr   = tmp.latg(:,idx,:);
                            tmp.lonr   = tmp.long(:,idx,:);
                            tmp.depthr = tmp.depthg(:,idx,:);
                            try
                                tmp.datar = tmp.data(:,idx,:,rcrd);
                            catch
                                tmp.datar = tmp.data(:,idx,:);
                            end
                            if length(idx)==2
                                tmp.latr   = squeeze(nanmean(tmp.latr,2));
                                tmp.lonr   = squeeze(nanmean(tmp.lonr,2));
                                tmp.depthr = squeeze(nanmean(tmp.depthr,2));
                                tmp.datar  = squeeze(nanmean(tmp.datar,2));
                            else
                                tmp.latr   = squeeze(tmp.latr);
                                tmp.lonr   = squeeze(tmp.lonr);
                                tmp.depthr = squeeze(tmp.depthr);
                                tmp.datar  = squeeze(tmp.datar);
                            end
                        elseif strcmp(choice,'lon');
                        % Get reduced matrix to feed into scatteredInterpolant
                            idx = find(abs([tmp.lon-deg]) == min(abs([tmp.lon-deg])));
                            tmp.latr   = tmp.latg(idx,:,:);
                            tmp.lonr   = tmp.long(idx,:,:);
                            tmp.depthr = tmp.depthg(idx,:,:);
                            try
                                tmp.datar = tmp.data(idx,:,:,rcrd);
                            catch
                                tmp.datar = tmp.data(idx,:,:); 
                            end
                            if length(idx)==2
                                tmp.latr   = squeeze(nanmean(tmp.latr,1));
                                tmp.lonr   = squeeze(nanmean(tmp.lonr,1));
                                tmp.depthr = squeeze(nanmean(tmp.depthr,1));
                                tmp.datar  = squeeze(nanmean(tmp.datar,1));
                            else
                                tmp.latr   = squeeze(tmp.latr);
                                tmp.lonr   = squeeze(tmp.lonr);
                                tmp.depthr = squeeze(tmp.depthr);
                                tmp.datar  = squeeze(tmp.datar);
                            end
                        end
                        % build interpolant and interpolate
                        if strcmp(choice,'lat');
                            tmp.lonr      = double(tmp.lonr(:));
                            tmp.depthr    = double(tmp.depthr(:));
                            tmp.datar     = double(tmp.datar(:)); 
                            F             = scatteredInterpolant(tmp.lonr(~isnan(tmp.datar)),...
                                            tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
                            tmp.out{rcrd} = F(tmpdeg(:,:),tmpdepth(:,:));
                        elseif strcmp(choice,'lon');
                            tmp.latr      = double(tmp.latr(:));
                            tmp.depthr    = double(tmp.depthr(:));
                            tmp.datar     = double(tmp.datar(:));
                            F             = scatteredInterpolant(tmp.latr(~isnan(tmp.datar)),...
                                            tmp.depthr(~isnan(tmp.datar)),tmp.datar(~isnan(tmp.datar)));
                            tmp.out{rcrd} = F(tmpdeg(:,:),tmpdepth(:,:));
                        end
                    end % end rcrd-loop
                    % Save results
                    obj.diag.(vars{ff})(i).slice = cat(ndims(tmp.out{1})+1,tmp.out{:});
                    obj.diag.(vars{ff})(i).units = curVar.units{i};
                    obj.diag.(vars{ff})(i).name  = curVar.name{i};
                    % Check if slice is already filled
                    if isempty(obj.slice)
                        if dmsn == obj.grid.nx;
                            obj.slice.coord = 'latitude';
                        elseif dmsn == obj.grid.ny;
                            obj.slice.coord = 'longitude';
                        end
                        fill_slice = 1;
                    else
                        fill_slice = 0;
                    end
                    
                    if fill_slice    
                        % Organize slice output
                        obj.slice.deg   = permute(repmat(tmpdeg,[1 1 nt]),[1 2 3]);
                        obj.slice.depth = permute(repmat(tmpdepth,[1 1 nt]),[1 2 3]);
                        obj.slice.sect  = deg;
                    end
                end % end i-loop
            end % end var-loop
        end % end method sliceDiag

        %--------------------------------------------------------------------------------
        function [figs,cbs] = sliceCmp(obj,dat1,dat2,t,varargin)
            % ------------------
            % A way to quickly plot slice comparisons
            % Produces 3 plots: dat1 and dat2 fields with the same colorbar, 
            % and a 3rd plot of the difference between them (dat1 - dat2)
            %
            % Usage:
            % - [figs,cbs] = sliceCmp(obj,dat1,dat2,t,varargin);
            %
            % Inputs:
            % - dat1: 2D field to plot (prepare it before using the script)
            % - dat2: same same but different 2D field to plot (prepare it before using the script)
            % - t:    time dimension to plot (0 == average)    
            %
            % Options:
            % - slice: if several slices were made, specify which one you want to plot here
            % - xlims: hard-coded x-limits (degrees)
            % - zlims: hard-coded z-limits (meters)
            % - cmap:  colormaps (if used, must be a cell array of length (vars))
            % - fontsize:   default 10
            % - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
            % - figdim:     default 1 (same height as width, its a multiplier)
            % - prc:        percentage to limit colorbar axes
            % - bal:        balance colorbar around 0 (1 == yes, 0 == no)
            % - levels:     hard-coded levels to plot (can't be used with A.prc)
            %
            % Example:
            % - [figs,cbs] = sliceCmp(obj,dat1,dat2,0,'xlims',[140 260],'zlims',[-500 0]);
            % This will take the average of dat1 + dat2 and plot them as well as their differences
            % ------------------

            % User-inputs
            romsOpt;
            A.xlims      = [];
            A.zlims      = [];
            A.cmap       = [];
            A.fontsize   = fontsize;
            A.figtype    = figtype;
            A.figdim     = 0.33;
            A.prc        = 2;
            A.bal        = 0;
            A.levels     = [];
            A.difflevels = [];
            A = romsObj.parse_pv_pairs(A,varargin);

            % Set levels if not supplied
            if isempty(A.levels);
                alldat  = [dat1(:) dat2(:)];
                lvls     = romsMaster.prclims(alldat,'prc',A.prc,'bal',A.bal);
                A.levels = linspace(lvls(1),lvls(2),20);
            end

            % Make figs(1) and figs(2)
            [figs(1),cbs(1)] = slicePlot(obj,dat1,t,...
                'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);
            [figs(2),cbs(2)] = slicePlot(obj,dat2,t,...
                'xlims',A.xlims,'zlims',A.zlims,'cmap',A.cmap,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels);

            % Get differences
            diff_dat  = dat1 - dat2;
            if isempty(A.difflevels)
                A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1);
                A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
            end

            % Make figs(3), difference
            [figs(3),cbs(3)] = slicePlot(obj,diff_dat,t,...
                'xlims',A.xlims,'zlims',A.zlims,'cmap','balance','fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels);

        end % end method sliceCmp

        %--------------------------------------------------------------------------------
        function obj = loadDiag(obj,vars,depth,varargin)
            % ------------------
            % Loads and interpolates monthly 2D validation data to the ROMS grid
            % See obj.paths.diag for available fields
            %
            % Usage:
            % - obj = loadDiag(obj,vars,depth,varargin);
            %
            % Inputs:
            % - vars  = validation variable(s) to load, as a cell array; paths set in getDiagpaths 
            % - depth = depths to interpolate data, will round to nearest z_avg_dep
            %            for surface data, use depth = 0
            %
            % Optional inputs (varargin):
            % - outer = degrees of lon/lat to pad interpolant, default = 3
            %
            % Example:
            % - obj = loadDiag(obj,vars,depth);
            % -------------------
            disp('loadDiag: Loading validation data');

            % Optional arguments
            A.outer = [3]; % interpolant padding (degrees)
            A=romsObj.parse_pv_pairs(A,varargin);

            % Check inputs
            diagfields = fieldnames(obj.paths.diag); 
            for i = 1:length(vars)
                if ~strcmp(vars{i},diagfields) & ~strcmp(vars{i},upper(diagfields));
                    disp(['    ERROR(loadDiag): ',vars{i},' is not a diagnostic variable']);
                    disp(['    Filling with NaN']);
                    obj.diag.(vars{i}).slice = nan(obj.grid.nx,obj.grid.ny,length(depth),12);
                    obj.diag.(vars{i}).name  = ' ';
                    obj.diag.(vars{i}).units = ' ';
                    obj.diag.(vars{i}).depth = depth;
                    skip(i) = 1;
                else
                    skip(i) = 0;
                end
            end
            if min(depth) < 0 | max(depth) > max(obj.grid.z_avg_dep)
                disp('    ERROR(loadDiag): Check depth input');
                return
            end
            for i = 1:length(depth)
                diffd = abs(depth(i) - obj.grid.z_avg_dep);
                ind   = find(diffd == min(diffd));
                depth(i) = obj.grid.z_avg_dep(ind);
                disp(['    Closest depth = ',num2str(depth(i))])
            end
            
            % Get temporary output longitude and latitude
            tmp.outlon = obj.grid.lon_rho;
            tmp.outlat = obj.grid.lat_rho;

            % Process each variable
            for i = 1:length(vars)        
                % Skip empty data
                if skip(i) == 1
                    continue
                end
                % Get path and coords for current variable
                curVar  = obj.paths.diag.(vars{i});

                % Go through multiple files, if they exist
                for j = 1:length(curVar.file)
                    
                    % Get coordinates
                    if strcmp(curVar.type{j},'nc');
                        tmp.lon = romsObj.lon360(ncread(curVar.file{j},curVar.lon{j}));
                        tmp.lat = ncread(curVar.file{j},curVar.lat{j});
                        if ndims(tmp.lon)==1
                            [tmp.lat,tmp.lon] = meshgrid(tmp.lat,tmp.lon);
                        elseif ndims(tmp.lon)>2
                            tmp.lon = squeeze(tmp.lon(:,:,1));
                            tmp.lat = squeeze(tmp.lat(:,:,1));
                            if ismember(curVar.dim{j},{'yxt','yxz','yxzt'})
                                tmp.lon = tmp.lon';
                                tmp.lat = tmp.lat';
                            end
                        end
                    elseif strcmp(curVar.type{j},'mat');
                        tmp.lon = load(curVar.file{j},curVar.lon{j});
                        tmp.lat = load(curVar.file{j},curVar.lat{j});
                        tmp.lon = tmp.lon.(curVar.lon{j});
                        tmp.lat = tmp.lat.(curVar.lat{j});
                        tmp.lon(tmp.lon<0) = tmp.lon(tmp.lon<0)+360;
                        if ndims(tmp.lon)==3
                            tmp.lat = squeeze(tmp.lat(:,:,1));
                            tmp.lon = squeeze(tmp.lon(:,:,1));
                        end
                    end

                    % Make lon/lat data gridded if it isnt setup that way
                    if sum(size(tmp.lon(:,:))==1) > 0
                        [tmp.lat,tmp.lon] = meshgrid(tmp.lat, tmp.lon);
                    end;
                    tmp = romsObj.struct2double(tmp);
                
                    % Get indeces for reduced domain interpolation
                    % Note: this breaks if ROMS boundary longitude is close to 0 or 360
                    idx = find(tmp.lon(:) > obj.grid.minlon_rho-A.outer ...
                             & tmp.lon(:) < obj.grid.maxlon_rho+A.outer ...
                             & tmp.lat(:) > obj.grid.minlat_rho-A.outer ...
                             & tmp.lat(:) < obj.grid.maxlat_rho+A.outer);

                    % Initialize interpolated output
                    obj.diag.(vars{i})(j).slice = nan(obj.grid.nx,obj.grid.ny,length(depth),12);

                    % Interpolate data on ROMS coords for each month
                    lvls = length(depth);
                    skip(i) = 0;
                    for k = 1:12
                        tmpdata = []; tmpdepth = [];

                        % Get data
                        if strcmp(curVar.type{j},'nc')
                            if strcmp(curVar.dim{j},'xyt') | strcmp(curVar.dim{j},'yxt');
                                tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
                                            [1,1,k],[inf,inf,1]));
                                tmpdepth = 0;
                                if strcmp(curVar.dim{j},'yxt');
                                    tmpdata = permute(tmpdata,[2 1 3]);
                                end
                            elseif strcmp(curVar.dim{j},'txy');
                                tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
                                            [k,1,1],[1,inf,inf]));
                                tmpdepth = 0;
                            elseif strcmp(curVar.dim{j},'xyzt')
                                tmpdepth = ncread(curVar.file{j},curVar.zvar{j});
                                tmpdata  = squeeze(ncread(curVar.file{j},curVar.var{j},...
                                       [1,1,1,k],[inf,inf,inf,1]));
                            elseif strcmp(curVar.dim{j},'xyz') | strcmp(curVar.dim{j},'yxz')
                                skip(i) = 1;
                                tmpdata = ncread(curVar.file{j},curVar.var{j});
                                tmpdepth = ncread(curVar.file{j},curVar.zvar{j});
                                if strcmp(curVar.dim{j},'yxz')
                                    tmpdata = permute(tmpdata,[2 1 3]);
                                    tmpdepth = permute(tmpdepth,[2 1 3]);
                                end
                            end
                        elseif strcmp(curVar.type{j},'mat')
                            if strcmp(curVar.dim{j},'xyt') | strcmp(curVar.dim{j},'yxt');
                                tmpdata  = load(curVar.file{j},curVar.var{j});
                                tmpdata  = tmpdata.(curVar.var{j})(:,:,k);
                                tmpdepth = 0;
                                if strcmp(curVar.dim{j},'yxt');
                                    tmpdata = permute(tmpdata,[2 1 3]);
                                end
                            elseif strcmp(curVar.dim{j},'xyzt')
                                tmpdepth = load(curVar.file{j},curVar.zvar{j});
                                tmpdepth = tmpdepth.(curVar.zvar{j});
                                tmpdata  = load(curVar.file{j},curVar.var{j});
                                tmpdata  = tmpdata.(curVar.var{j})(:,:,:,k);
                            elseif strcmp(curVar.dim{j},'xyz') | strcmp(curVar.dim{j},'yxz');
                                skip(i) = 1;
                                tmpdepth = load(curVar.file{j},curVar.zvar{j});
                                tmpdepth = tmpdepth.(curVar.zvar{j});
                                tmpdata  = load(curVar.file{j},curVar.var{j});
                                tmpdata  = tmpdata.(curVar.var{j});
                                if strcmp(curVar.dim{j},'yxz');
                                    tmpdepth = permute(tmpdepth,[2 1 3]);
                                    tmpdata  = permute(tmpdata,[2 1 3]);
                                end
                            end
                        end

                        % Apply any factors (eg nMol to uMol)    
                        if isfield(curVar,'factor');
                            tmpdata = tmpdata*curVar.factor{j};
                        end
                
                        % Fix tmpdepth?
                        if ndims(tmpdepth)==3
                            tmpdepth = squeeze(tmpdepth(1,1,:));
                        end

                        % Build interpolant and interpolate over all depths
                        for z = 1:lvls
                            data  = tmpdata;
                            zdepth = tmpdepth; 
                            if nanmean(depth) < 0
                                depth = -depth;
                            end
                            zind  = find(abs([zdepth-depth(z)]) == min(abs([zdepth-depth(z)])));
                            if ~strcmp(curVar.dim{j},'xyt') | ~strcmp(curVar.dim{j},'yxt');
                                data  = data(:,:,zind);
                            end
                            if strcmp(curVar.dim{j},'xyt') & z > 1 | strcmp(curVar.dim{j},'yxt') & z > 1
                                disp(['    NOTE(loadDiag): No z-data for ',vars{i},', skipping']);
                                continue
                            end
                            % Interpolate
                            F = scatteredInterpolant(double(tmp.lon(idx)),double(tmp.lat(idx)),...
                                                     double(data(idx)),'linear','nearest');
                            tmpout{z,k} = F(double(tmp.outlon),double(tmp.outlat));
                            tmpout{z,k}(isnan(obj.grid.mask_rho)) = nan;
                        end % end z-loop
                    end % end k-loop
                    % Save interpolated data
                    for ll = 1:length(depth);
                        obj.diag.(vars{i})(j).slice(:,:,ll,:)  = cat(3, tmpout{ll,:});
                        obj.diag.(vars{i})(j).slice(:,:,ll,:)  = single(obj.diag.(vars{i})(j).slice(:,:,ll,:));
                    end
                    obj.diag.(vars{i})(j).slice = squeeze(obj.diag.(vars{i})(j).slice);
                    obj.diag.(vars{i})(j).depth = depth;
                    obj.diag.(vars{i})(j).units = curVar.units{j};
                    obj.diag.(vars{i})(j).name  = curVar.name{j};
                end % end j-loop
            end % end i-loop
        end % end method loadDiag

        %--------------------------------------------------------------------------------
        function [obj] = clearROMS(obj) 
            % ----------------------
            % Clears loaded data and coordinates
            %
            % Usage:
            % - [obj] = clearROMS(obj) 
            % ----------------------
            % Clear fields
            obj.slice = [];
            obj.profile = [];
            obj.diag = [];
            obj.data = [];

            % Clear file-specific grid data
            file_types = {'avg','his'};
            for ff = 1:length(file_types)
                obj.grid.(file_types{ff}) = [];
            end
        end % end method clearROMS
    end % end object methods

    %----------------------------------------------------------------------------------------
    % Utility functions (static methods)
    methods (Static)
        %--------------------------------------------------------------------------------
        function x = lon360(x)
            % ------------------
            % - corrects negative longitudes such that all lon are between 0:360
            % ------------------
            x(x<0) = x(x<0) + 360;
        end % end static method lon360

        %--------------------------------------------------------------------------------
        function x = lon180(x)
            % ------------------
            % - corrects longitudes such that all lon are between -180:180  
            % ------------------
            x(x>180) = x(x>180) - 360;
        end % end static method lon180

        %--------------------------------------------------------------------------------
        function x = struct2double(x)
            % ------------------
            % - converts structure fields to double  
            % ------------------
            ffields = fields(x);
            for ff = 1 : length(ffields)
                    if isa(x.(ffields{ff}),'double')
                        x.(ffields{ff}) = double(x.(ffields{ff}));
                    end
            end
        end % end static method struct2double

        %--------------------------------------------------------------------------------
        function x = struct2single(x)
            % ------------------
            % - converts structure fields to single  
            % ------------------
            ffields = fields(x);
            for ff = 1 : length(ffields)
                    if isa(x.(ffields{ff}),'double')
                        x.(ffields{ff}) = single(x.(ffields{ff}));
                    end
            end
        end % end static method struct2single

        %--------------------------------------------------------------------------------
        function [ws,wsc] = WindStress(u,v,lon,lat,ang)
            % ------------------
            % - calculates wind stress and wind stress curl from zonal/meridional wind stress 
            % ------------------
            % Get rho points 
            u = romsObj.u2rho(u);
            v = romsObj.v2rho(v);
            % - calculate wind stress
            ws = sqrt(u.*u + v.*v);
            % - check that lon/lat are gridded
            [a,b] = size(lon);
            if a == 1 | b == 1
                [lon,lat] = meshgrid(lon,lat);
            end
            % - check for yearly or monthly file
            [a,b,c] = size(u);
            % - calculate wind stress curl
            E   = cos(ang).*u - sin(ang).*v;
            N   = sin(ang).*u + cos(ang).*v;
            wsc = ws.*NaN;
            for i = 2:size(lon,1)-1;
                for j = 2:size(lon,2)-1;
                    dx = sw_dist([lat(i+1,j) lat(i-1,j)],...
                             [lon(i+1,j) lon(i-1,j)],'km')*1000;
                    dy = sw_dist([lat(i,j+1) lat(i,j-1)],...
                             [lon(i,j+1) lon(i,j-1)],'km')*1000;
                    for t = 1:c
                        wsc(i,j,t) = (N(i+1,j,t)-N(i-1,j,t))/dx - ...
                                 (E(i,j+1,t)-E(i,j-1,t))/dy ;
                    end
                end
            end
            end % end static method WindStress
        %--------------------------------------------------------------------------------
        function dis = gc_dist(lon1,lat1,lon2,lat2); 
            % -------------------
            % - Distance between 2 points along a great circle
            % - Lat/lon in Radians
            % - Jeroen Molemaker, UCLA 2008
            % -------------------
            dlat = lat2-lat1;
            dlon = lon2-lon1;
            dang = 2*asin( sqrt( sin(dlat/2).^2 + cos(lat2).*cos(lat1).*sin(dlon/2).^2 ) );  %% haversine function
            r_earth = 6371315.;
            dis = r_earth*dang;
        end % end static method gc_dist

        %--------------------------------------------------------------------------------
        function lims = prclims(data,varargin)
            % -------------------
            % Script to automatically create limits of data to show
            % i.e. xlim( ) or BinLimits( ) 
            %
            % Usage:
            %   lims = prclims(data,varargin)
            %
            % Inputs:
            %   prc -- percentile to limit data by (default 2% to 98%);
            %   bal -- balance limits (i.e. for histograms around 0)
            %          0 = auto-limits, 1 = evenly balanced around 0, 2 = set min to 0
            % -------------------
            % Default arguments
            A.prc = [2]; % default range 2 - 98%
            A.bal = [1]; 
            A     = romsObj.parse_pv_pairs(A,varargin);
            % Get limits
            low = prctile(data(:),[A.prc]);
            hih = prctile(data(:),[100-A.prc]);
            % If low < 0 and hih > 0, check for balance
            if low < 0 & hih > 0 & A.bal ~= 2
                lims = [-(max(abs([low hih]))) max(abs([low hih]))];
            elseif A.bal == 1
                lims = [-(max(abs([low hih]))) max(abs([low hih]))];
            elseif A.bal == 2
                lims = [0 hih];
            else
                lims = [low hih];
            end
        end % end static method prclims

        %--------------------------------------------------------------------------------
        function [outvar] = grid_to_grid(inlon,inlat,invar,outlon,outlat)
            % --------------------------------------------------------------------
            % Function to grid one meshgrid to another.
            % Only works if 'inlon/inlat' encapsulates 'outlon/outlat'
            % i.e. a subgrid inside a larger grid.
            %
            % Usage:
            % - outvar = grid_to_grid(inlon,inlat,invar,outlon,outlat)
            %
            % Inputs:
            % - inlon  = meshgrid of lon points (m x n)
            % - inlat  = meshgrid of lat points (m x n)
            % - invar  = meshgrid of data on surface (m x n x t)
            % - outlon = meshgrid of lat points (mm x nn)
            % - outlat = meshrgid of lat points (mm x nn)
            %
            % Outputs:
            % - outvar = meshgrid of data on surface (mm x nn x t)
            %
            % NOTES:
            % Cannot be used with 4D grids, reduce your grid to 3D in a loop first
            % --------------------------------------------------------------------
            % get indices for reduced domain interpolation
            minlon_rho = min(outlon(:));
            maxlon_rho = max(outlon(:));
            minlat_rho = min(outlat(:));
            maxlat_rho = max(outlat(:));
            outer      = 3;
            idx = find(inlon(:) > minlon_rho-outer ...
                 & inlon(:) < maxlon_rho+outer ...
                 & inlat(:) > minlat_rho-outer ...
                 & inlat(:) < maxlat_rho+outer);
            % interpolate over all time steps
            [x]  = size(invar);
            [xi] = length(x);
            if xi == 2 % no time
                tmpvar = invar;
                F      = scatteredInterpolant(double(inlon(idx)),double(inlat(idx)),...
                                              double(tmpvar(idx)),'linear','nearest');
                outvar = F(double(outlon),double(outlat));
            else 
                for t = 1:x(end)
                    if xi == 3
                        tmpvar = squeeze(invar(:,:,t));
                    elseif xi == 4
                        tmpvar = squeeze(invar(:,:,t));
                    end
                    % - build interpolant and interpolate
                    F           = scatteredInterpolant(double(inlon(idx)),double(inlat(idx)),...
                                                     double(tmpvar(idx)),'linear','nearest');
                    tmpout{t} = F(double(outlon),double(outlat));
                end
                outvar = cat(3,tmpout{:});
            end
        end % end static method grid_to_grid

        %--------------------------------------------------------------------------------
        function [lambda2,xi,ST,SN]=okubo_weiss(u,v,pm,pn); 
            % --------------------------------------------------------------------
            % Compute the OkuboWeiss parameter
            %
            % Usage:
            %   [lambda,x,ST,SN] = okubo_weiss(u,v,pm,pn);
            %
            % Inputs:
            %   u = ROMS u velocity
            %   v = ROMS v velocity
            %   pm = ROMS grid 'pm', or curvilinear coordinate metrix in XI
            %   pn = ROMS grid 'pn', or curvilinear coordinate metrix in ETA
            %
            % Outputs:
            %    lambda2 = Okubo-Weiss (lambda^2) in s^-2
            %   xi      = Relative vorticity in s^-1
            %   ST      = Shear strain in s^-1
            %   SN      = Normal strain in s^-1
            % --------------------------------------------------------------------
            % Get grid dimensions        
            [Mp,Lp]=size(pm);
            L=Lp-1; M=Mp-1;
            Lm=L-1; Mm=M-1;
            % Get u/v dimensions
            nt = size(u,3);
            % Initialize output arrays
            xi      = zeros(Mp,Lp,nt);
            ST      = xi;
            SN      = xi;
            lambda2 = xi;
            % Go through each time entry
            for i = 1:nt
                mn_p = zeros(M,L);
                uom  = zeros(M,Lp);
                von  = zeros(Mp,L);
                uom=2*u(1:M,:,i)./(pm(1:M,:)+pm(2:Mp,:));
                uon=2*u(1:M,:,i)./(pn(1:M,:)+pn(2:Mp,:));
                von=2*v(:,1:L,i)./(pn(:,1:L)+pn(:,2:Lp));
                vom=2*v(:,1:L,i)./(pm(:,1:L)+pm(:,2:Lp));
                mn=pm.*pn;
                mn_p=(mn(1:M,1:L)+mn(1:M,2:Lp)+...
                      mn(2:Mp,2:Lp)+mn(2:Mp,1:L))/4;
                % relative vorticity
                xi(:,:,i) = mn.*romsObj.psi2rho(von(2:Mp,:)-von(1:M,:)-uom(:,2:Lp)+uom(:,1:L));
                % Sigma_T
                ST(:,:,i) = mn.*romsObj.psi2rho(von(2:Mp,:)-von(1:M,:)+uom(:,2:Lp)-uom(:,1:L));
                % Sigma_N
                SN(2:end-1,2:end-1,i) = mn(2:end-1,2:end-1).*(uon(2:end,2:end-1)...
                                      -uon(2:end,1:end-2)...
                                      -vom(2:end-1,2:end)...
                                      +vom(1:end-2,2:end));
                % Lambda^2
                lambda2(:,:,i) = SN(:,:,i).^2 + ST(:,:,i).^2 - xi(:,:,i).^2;
            end
        end % end static method okubo_weiss

        %--------------------------------------------------------------------------------
        function [var_rho] = psi2rho(var_psi)
            % --------------------------------------------------------------------
            % Converts a field at psi points to the rho points
            %
            % Usage:
            % - [var_rho] = psi2rho(var_psi)
            % --------------------------------------------------------------------
            % Convert
            [M,L]=size(var_psi);
            Mp=M+1;
            Lp=L+1;
            Mm=M-1;
            Lm=L-1;
            var_rho=zeros(Mp,Lp);
            var_rho(2:M,2:L)=0.25*(var_psi(1:Mm,1:Lm)+var_psi(1:Mm,2:L)+...
                                   var_psi(2:M,1:Lm)+var_psi(2:M,2:L));
            var_rho(1,:)=var_rho(2,:);
            var_rho(Mp,:)=var_rho(M,:);
            var_rho(:,1)=var_rho(:,2);
            var_rho(:,Lp)=var_rho(:,L);
        end % end static method psi2rho

        %--------------------------------------------------------------------------------
        function [var_rho] = u2rho(var_u)
            % --------------------------------------------------------------------
            % Transfer a field at u points to the rho points
            %
            % Usage:
            % - [var_rho] = u2rho(var_u)
            % --------------------------------------------------------------------    
            % Convert
            dimensions = size(var_u);
            M  = dimensions(1);
            Lp = dimensions(2);
            Mp = M+1;
            Mm = M-1;
            if length(dimensions)==2
                var_rho = zeros(Mp,Lp);
                var_rho(2:M,:) = 0.5.*(var_u(1:Mm,:)+var_u(2:M,:));
                var_rho(1,:) = var_rho(2,:);
                var_rho(Mp,:) = var_rho(M,:);
            elseif length(dimensions)==3
                var_rho = zeros(Mp,Lp,dimensions(3));
                var_rho(2:M,:,:) = 0.5.*(var_u(1:Mm,:,:)+var_u(2:M,:,:));
                var_rho(1,:,:) = var_rho(2,:,:);
                var_rho(Mp,:,:) = var_rho(M,:,:);
            elseif length(dimensions)==4
                var_rho = zeros(Mp,Lp,dimensions(3),dimensions(4));
                var_rho(2:M,:,:,:) = 0.5.*(var_u(1:Mm,:,:,:)+var_u(2:M,:,:,:));
                var_rho(1,:,:,:) = var_rho(2,:,:,:);
                var_rho(Mp,:,:,:) = var_rho(M,:,:,:);
            end
        end % end static method u2rho

        %--------------------------------------------------------------------------------
        function [var_rho] = v2rho(var_v)
            % --------------------------------------------------------------------
            % Transfer a field at v points to the rho points
            %
            % Usage:
            % - [var_rho] = v2rho(var_v)
            % --------------------------------------------------------------------    
            % Convert
            dimensions = size(var_v);
            Mp = dimensions(1);
            L  = dimensions(2);
            Lp = L+1;
            Lm = L-1;
            if length(dimensions)==2
                var_rho = zeros(Mp,Lp);
                var_rho(:,2:L) = 0.5.*(var_v(:,1:Lm)+var_v(:,2:L));
                var_rho(:,1) = var_rho(:,2);
                var_rho(:,Lp) = var_rho(:,L);
            elseif length(dimensions)==3
                var_rho = zeros(Mp,Lp,dimensions(3));
                var_rho(:,2:L,:) = 0.5.*(var_v(:,1:Lm,:)+var_v(:,2:L,:));
                var_rho(:,1,:) = var_rho(:,2,:);
                var_rho(:,Lp,:) = var_rho(:,L,:);
            elseif length(dimensions)==4
                var_rho = zeros(Mp,Lp,dimensions(3),dimensions(4));
                var_rho(:,2:L,:,:) = 0.5.*(var_v(:,1:Lm,:,:)+var_v(:,2:L,:,:));
                var_rho(:,1,:,:) = var_rho(:,2,:,:);
                var_rho(:,Lp,:,:) = var_rho(:,L,:,:);
            end
        end % end static method v2rho

        %--------------------------------------------------------------------------------
        function [matrixOut] = spatial_filter(matrixIn,Nr,Nc)
            % --------------------------------------------------------------------
            % Smooths 2D array data, ignores NaNs.
            %    
            % This function smooths the data in matrixIn using a mean filter over a
            % rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
            % element "i" by the mean of the rectange centered on "i".  Any NaN
            % elements are ignored in the averaging.  If element "i" is a NaN, then it
            % will be preserved as NaN in the output.  At the edges of the matrix,
            % where you cannot build a full rectangle, as much of the rectangle that
            % fits on your matrix is used (similar to the default on Matlab's builtin
            % function "smooth").
            %
            % Usage:
            % - matrixOut = spatial_filter(matrixIn,Nr,Nc)
            % 
            % Inputs:
            % - matrixIn: original matrix
            % - Nr      : number of points used to smooth rows
            % - Nc      : number of points to smooth columns.  If not specified, Nc = Nr.
            % 
            % Outputs:
            % - matrixOut: smoothed version of original matrix

            % Initial error statements and definitions
            if nargin < 2, error('Not enough input arguments!'), end

            N(1) = Nr;
            if nargin < 3, N(2) = N(1); else N(2) = Nc; end

            if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
            if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

            % Building matrices that will compute running sums.  The left-matrix, eL,
            % smooths along the rows.  The right-matrix, eR, smooths along the
            % columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
            % (2*Nc+1) rectangle centered on element "i".
            [row,col] = size(matrixIn);
            eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
            eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

            % Setting all "NaN" elements of "matrixIn" to zero so that these will not
            % affect the summation.  (If this isn't done, any sum that includes a NaN
            % will also become NaN.)
            A = isnan(matrixIn);
            matrixIn(A) = 0;

            % For each element, we have to count how many non-NaN elements went into
            % the sums.  This is so we can divide by that number to get a mean.  We use
            % the same matrices to do this (ie, "eL" and "eR").
            nrmlize = eL*(~A)*eR;
            nrmlize(A) = NaN;

            % Actually taking the mean.
            matrixOut = eL*matrixIn*eR;
            matrixOut = matrixOut./nrmlize;
        end % end static method spatial_filter
    
        %--------------------------------------------------------------------------------
        function [o2_sat] = o2_sat(pT,s)
            % --------------------------------------------------------------------
            % Saturated O2 in mmol/m3
            %
            % Usage:
            %    [o2_sat] = o2_sat(pT,s)
            %
            % Inputs:
            %    pT = potential temperature
            %    s  = salinity;
            % --------------------------------------------------------------------
            % Constants
            a_0 = 2.00907;
            a_1 = 3.22014;
            a_2 = 4.05010;
            a_3 = 4.94457;
            a_4 = -2.56847E-1;
            a_5 = 3.88767;
            b_0 = -6.24523E-3;
            b_1 = -7.37614E-3;
            b_2 = -1.03410E-2;
            b_3 = -8.17083E-3;
            c_0 = -4.88682E-7;
            TS = log(((273.16+25.0)-pT)./(273.16+pT));
            O2SAT = exp(a_0+TS.*(a_1+TS.*(a_2+TS.*(a_3+TS.*(a_4+TS.*a_5)))) + ...
                    s.*((b_0+TS.*(b_1+TS.*(b_2+TS.*b_3))) + s.*c_0));
            o2_sat = O2SAT.*44.6596;
        end % end static method o2_sat

        %--------------------------------------------------------------------------------
        function [n2o_sat] = n2o_sat(pT,s)
            % --------------------------------------------------------------------
            % Saturated N2O in mmol/m3
            %
            % Usage:
            %    [n2o_sat] = n2o_sat(pT,s)
            %
            % Inputs:
            %    pT = potential temperature
            %    s  = salinity;
            % --------------------------------------------------------------------
            % Constants
            a_1 = -165.8802;
            a_2 = 222.8743;
            a_3 = 92.0792;
            a_4 = -1.48425;
            b_1 = -0.056235;
            b_2 = 0.031619;
            b_3 = -0.0048472;
            TS = 273.16 + pT;
            N2OSAT =  exp(a_1+a_2.*(100.0./TS)+a_3.*log(TS./100.0)+ ...
                          a_4.*(TS./100.0).^2+s.*(b_1+b_2.*(TS./100.0)+b_3.*(TS./100.0).^2));
            n2o_sat = N2OSAT.*1000.*1000.*(300.0*1e-9);
        end % end static method n2o_sat

        %--------------------------------------------------------------------------------
        function params = parse_pv_pairs(params,pv_pairs) 
            % --------------------------------------------------------------------
            npv = length(pv_pairs);
            n = npv/2;

            if n~=floor(n)
              error 'Property/value pairs must come in PAIRS.'
            end
            if n<=0
              % just return the defaults
              return
            end

            if ~isstruct(params)
              error 'No structure for defaults was supplied'
            end

            % there was at least one pv pair. process any supplied
            propnames = fieldnames(params);
            lpropnames = lower(propnames);
            for i=1:n
              p_i = lower(pv_pairs{2*i-1});
              v_i = pv_pairs{2*i};

              ind = strmatch(p_i,lpropnames,'exact');
              if isempty(ind)
                ind = find(strncmp(p_i,lpropnames,length(p_i)));
                if isempty(ind)
                  error(['No matching property found for: ',pv_pairs{2*i-1}])
                elseif length(ind)>1
                  error(['Ambiguous property name: ',pv_pairs{2*i-1}])
                end
              end
              p_i = propnames{ind};

              % override the corresponding default in params
              params = setfield(params,p_i,v_i);
            end
        end % end static method parse_pv_pairs

        %--------------------------------------------------------------------------------
        function fig = piofigs(pick,ymult) 
            % --------------------------------------------------------------------
            % pick options:
            % sfig = 90mm
            % mfig = 140mm
            % lfig = 190mm
            %
            % ymult = number to multiply xwidth by
            % Return errors
            if ymult <= 0
                disp('ymult must be greater than 0');
                return
            end
            % Grab dim
            if strmatch(pick,'sfig')
                xwidth = 9;
                str    = 'sfig: 9cm';
            elseif strmatch(pick,'mfig');
                xwidth = 14;
                str    = 'mfig: 14cm';
            elseif strmatch(pick,'lfig');
                xwidth = 19;
                str    = 'lfig: 19cm';
            end
            % Settings (90cm wide)
            disp(['    Initializing ',str])
            romsOpt;
            fig             = figure('visible',figVisible);
            fig.PaperUnits  = 'Centimeters';
            fig.Units       = 'Centimeters';
            fig.Position(3) = xwidth;
            fig.Position(4) = xwidth*ymult;
            fig.Color       = [1 1 1];
            if strcmp(figVisible,'off')
                disp('    Done (invisible)');
            elseif strcmp(figVisble','on');
                disp('    Done');
            end
        end % end static method piofigs

        %--------------------------------------------------------------------------------
        function pltjpg(pltnum) 
            % --------------------------------------------------------------------
            % Create a temporary figure in tmpfigsPath
            % Usage:
            % - pltjpg          <-- makes (and overrides) tmpfigs1.jpg
            % - pltjpg(pltnum)  <-- makes tmpfigs($pltnum).jpg
            % - pltjpg(0)       <-- remove all tmpfigs 
            % --------------------------------------------------------------------
            % Turn off warnings
            warning off
            % Get figpath
            romsOpt;
            A.figpath = tmpfigsPath;
            A.figname = 'tmpfig';
            if nargin < 1
                % Get number of tmpfigs
                temp   = dir(A.figpath);
                fnames = {temp.name};
                isdir  = {temp.isdir};
                fnames = fnames(cell2mat(isdir)==0);
                numfig = length(fnames);
                clear temp fnames isdir
                % Set figpath
                pltnum  = numfig + 1;
            end
            % Clean directory?
            if pltnum == 0
                if strcmp(figsFormat,'jpg');
                    cmd = ['rm ',tmpfigsPath,'/*.jpg'];
                    system(cmd);
                    return
                elseif strcmp(figsFormat,'png');
                    cmd = ['rm ',tmpfigsPath,'/*.png'];
                    system(cmd);
                    return
                elseif strcmp(figsFormat,'pdf');
                    cmd = ['rm ',tmpfigsPath,'/*.pdf'];
                    system(cmd);
                    return
                end
            end
            % Save figure as tmpfig
            figsPath = [A.figpath,A.figname,num2str(pltnum)];
            fig     = get(gcf);
            disp(['    Saving ',A.figname,num2str(pltnum)]);
            export_fig(figsFormat,figsPath,figsQuality);
        end % end static method pltjpg
    end % end static methods declarations
end % end classdef
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------------

%--------------------------------------------------------------------------------
function [constants] = getConstants
    % ----------------------
    % Save constants
    constants.mmolNs_TgNy  = [14*86400*365.25]/(1e15); % mmolN/s to TgN/yr
    constants.mmolN2s_TgNy = [28*86400*365.25]/(1e15); % mmolN2/s to TgN/yr
    constants.s_d          = 86400;                    % seconds to days
    constants.s_yr         = 86400.*365.25;            % seconds to years
    constants.mmol_umol    = 1000;                     % mmol to umol
    constants.umol_nmol    = 1000;                     % umol to nmol
    constants.mmol_nmol    = 1000*1000;                % mmol to nmol
    constants.nmol_umol    = 1/1000;                   % nmol to umol
    constants.umol_mmol    = 1/1000;                   % umol to mmol
    constants.nmol_mmol    = (1/1000)*(1/1000);        % nmol to mmol
end % end method getConstants 

%--------------------------------------------------------------------------------
function [paths] = getDiagPaths;    
    % -------------------
    % Initialize diagnostic plots: set metadata, save paths, get axis limits 
    % temperature
    romsOpt
    paths.diag.temp.file  = {[valiPath,'temp/woa18_temp_clim.nc'],...
                             [valiPath,'temp/GLODAPv2.2016b.temperature.nc']};
    paths.diag.temp.type  = {'nc','nc'};
    paths.diag.temp.var   = {'t_an','temperature'};
    paths.diag.temp.zvar  = {'depth','Depth'};
    paths.diag.temp.dim   = {'xyzt','xyz'};
    paths.diag.temp.lon   = {'lon','lon'};
    paths.diag.temp.lat   = {'lat','lat'};
    paths.diag.temp.name  = {'WOA-18','GLODAPv2'};
    paths.diag.temp.units = {'$^oC$','$^oC$'};

    % salinity
    paths.diag.salt.file  = {[valiPath,'salt/woa18_salt_clim.nc'],...
                             [valiPath,'salt/GLODAPv2.2016b.salinity.nc']};
    paths.diag.salt.type  = {'nc','nc'};
    paths.diag.salt.var   = {'s_an','salinity'};
    paths.diag.salt.zvar  = {'depth','Depth'};
    paths.diag.salt.dim   = {'xyz','xyz'};
    paths.diag.salt.lon   = {'lon','lon'};
    paths.diag.salt.lat   = {'lat','lat'};
    paths.diag.salt.name  = {'WOA-18','GLODAPv2'};
    paths.diag.salt.units = {'PSU','PSU'};

    % density
    paths.diag.rho.file  = {[valiPath,'rho/woa18_rho_clim.nc']};
    paths.diag.rho.type  = {'nc'};
    paths.diag.rho.var   = {'I_an'};
    paths.diag.rho.zvar  = {'depth'};
    paths.diag.rho.dim   = {'xyzt'};
    paths.diag.rho.lon   = {'lon'};
    paths.diag.rho.lat   = {'lat'};
    paths.diag.rho.name  = {'WOA-18 Density'};
    paths.diag.rho.units = {'kg m$^{-3}$'};

    % sea surface height
    paths.diag.SSH.file  = {[valiPath,'SSH/monthly_AVISO.mat']};
    paths.diag.SSH.type  = {'mat'};
    paths.diag.SSH.var   = {'adt_month_av'};
    paths.diag.SSH.zvar  = {[]};
    paths.diag.SSH.dim   = {'xyt'};
    paths.diag.SSH.lon   = {'lon_aviso'};
    paths.diag.SSH.lat   = {'lat_aviso'};
    paths.diag.SSH.name  = {'AVISO+ Absolute Dynamic Topography'};
    paths.diag.SSH.units = {'m'};

    % zonal wind stress
    paths.diag.sustr.file   = {[valiPath,'wind_stress/monthly_WSTRESS.mat']};
    paths.diag.sustr.type   = {'mat'};
    paths.diag.sustr.var    = {'u'};
    paths.diag.sustr.zvar   = {[]};
    paths.diag.sustr.dim    = {'xyt'};
    paths.diag.sustr.lon    = {'lon'};
    paths.diag.sustr.lat    = {'lat'};
    paths.diag.sustr.name   = {'SCOW-2010 Zonal Wind Stress'};
    paths.diag.sustr.units  = {'m$^{2}$ s$^{-2}$'};
    paths.diag.sustr.factor = {[1e-3]};

    % meridional wind stress
    paths.diag.svstr.file   = {[valiPath,'wind_stress/monthly_WSTRESS.mat']};
    paths.diag.svstr.type   = {'mat'};
    paths.diag.svstr.var    = {'v'};
    paths.diag.svstr.zvar   = {[]};
    paths.diag.svstr.dim    = {'xyt'};
    paths.diag.svstr.lon    = {'lon'};
    paths.diag.svstr.lat    = {'lat'};
    paths.diag.svstr.name   = {'SCOW-2010 Meridional Wind Stress'};
    paths.diag.svstr.units  = {'m$^{2}$ s$^{-2}$'};
    paths.diag.svstr.factor = {[1e-3]};
    
    % wind stress
    paths.diag.ws.file   = {[valiPath,'wind_stress/monthly_WSTRESS.mat']};
    paths.diag.ws.type   = {'mat'};
    paths.diag.ws.var    = {'ws'};
    paths.diag.ws.zvar   = {[]};
    paths.diag.ws.dim    = {'xyt'};
    paths.diag.ws.lon    = {'lon'};
    paths.diag.ws.lat    = {'lat'};
    paths.diag.ws.name   = {'SCOW-2010 Wind Stress'};
    paths.diag.ws.units  = {'m$^{2}$ s$^{-2}$'};
    paths.diag.ws.factor = {[1e-3]};

    % wind stress curl
    paths.diag.wsc.file  = {[valiPath,'wind_stress/monthly_WSTRESS.mat']};
    paths.diag.wsc.type  = {'mat'};
    paths.diag.wsc.var   = {'wsc'};
    paths.diag.wsc.zvar  = {[]};
    paths.diag.wsc.dim   = {'xyt'};
    paths.diag.wsc.lon   = {'lon'};
    paths.diag.wsc.lat   = {'lat'};
    paths.diag.wsc.name  = {'SCOW-2010 Wind Stress Curl'};
    paths.diag.wsc.units = {'N m$^{-2}$'};

    % u velocity
    paths.diag.u.file  = {[valiPath,'uv/GODAS_uv.mat']};
    paths.diag.u.type  = {'mat'};
    paths.diag.u.var   = {'U'};
    paths.diag.u.zvar  = {'dep'};
    paths.diag.u.dim   = {'xyzt'};
    paths.diag.u.lon   = {'lon'};
    paths.diag.u.lat   = {'lat'};
    paths.diag.u.name  = {'GODAS U-Velocity'};
    paths.diag.u.units = {'m s$^{-1}'}; 

    % v velocity
    paths.diag.v.file  = {[valiPath,'uv/GODAS_uv.mat']};
    paths.diag.v.type  = {'mat'};
    paths.diag.v.var   = {'V'};
    paths.diag.v.zvar  = {'dep'};
    paths.diag.v.dim   = {'xyzt'};
    paths.diag.v.lon   = {'lon'};
    paths.diag.v.lat   = {'lat'};
    paths.diag.v.name  = {'GODAS V-Velocity'};
    paths.diag.v.units = {'m s$^{-1}'}; 

    % Eddy kinetic energy
    paths.diag.EKE.file  = {[valiPath,'EKE/aviso_1993_2022_EKE_monthly_climatology.mat']};
    paths.diag.EKE.type  = {'mat'};
    paths.diag.EKE.var   = {'eke'};
    paths.diag.EKE.zvar  = {[]}; 
    paths.diag.EKE.dim   = {'xyt'};
    paths.diag.EKE.lon   = {'lon'};
    paths.diag.EKE.lat   = {'lat'};
    paths.diag.EKE.name  = {'AVISO+ Ssalto/Duacs EKE climatology (1993-2023)'};
    paths.diag.EKE.units = {'cm$^2$ s$^{-2}$'};

    % oxygen
    paths.diag.O2.file  = {[valiPath,'O2/woa18_o2_clim.nc'],...
                           [valiPath,'O2/GLODAPv2.2016b.oxygen.nc']};
    paths.diag.O2.type  = {'nc','nc'};
    paths.diag.O2.var   = {'o_an','oxygen'};
    paths.diag.O2.zvar  = {'depth','Depth'};
    paths.diag.O2.dim   = {'xyzt','xyz'};
    paths.diag.O2.lon   = {'lon','lon'};
    paths.diag.O2.lat   = {'lat','lat'};
    paths.diag.O2.name  = {'WOA-18','GLODAPv2'};
    paths.diag.O2.units = {'mmol O$_2$ m$^{-3}$','mmol O$_2$ m$^{-3}$'};

    % nitrate 
    paths.diag.NO3.file  = {[valiPath,'NO3/woa18_no3_clim.nc'],...
                            [valiPath,'NO3/GLODAPv2.2016b.NO3.nc']};
    paths.diag.NO3.type  = {'nc','nc'};
    paths.diag.NO3.var   = {'n_an','NO3'};
    paths.diag.NO3.zvar  = {'depth','Depth'};
    paths.diag.NO3.dim   = {'xyzt','xyz'};
    paths.diag.NO3.lon   = {'lon','lon'};
    paths.diag.NO3.lat   = {'lat','lat'};
    paths.diag.NO3.name  = {'WOA-18','GLODAPv2'};
    paths.diag.NO3.units = {'mmol N m$^{-3}$','mmol N m$^{-3}$'};
    
    % Copy for NOx (NO3 is typically reported via (NO3+NO2))
    paths.diag.NOX = paths.diag.NO3;

    % phosphate
    paths.diag.PO4.file  = {[valiPath,'PO4/woa18_po4_clim.nc'],...
                            [valiPath,'PO4/GLODAPv2.2016b.PO4.nc']};
    paths.diag.PO4.type  = {'nc','nc'};
    paths.diag.PO4.var   = {'p_an','PO4'};
    paths.diag.PO4.zvar  = {'depth','Depth'};
    paths.diag.PO4.dim   = {'xyzt','xyz'};
    paths.diag.PO4.lon   = {'lon','lon'};
    paths.diag.PO4.lat   = {'lat','lat'};
    paths.diag.PO4.name  = {'WOA-18','GLODAPv2'};
    paths.diag.PO4.units = {'mmol P m$^{-3}$','mmol P m$^{-3}$'};

    % N* (only 1.00 available)
    paths.diag.nstar.file  = {[valiPath,'nstar/woa18_nstar_clim.nc']};
    paths.diag.nstar.type  = {'nc'};
    paths.diag.nstar.var   = {'nstar'};
    paths.diag.nstar.zvar  = {'depth'};
    paths.diag.nstar.dim   = {'xyzt'};
    paths.diag.nstar.lon   = {'lon'};
    paths.diag.nstar.lat   = {'lat'};
    paths.diag.nstar.name  = {'WOA-18'};
    paths.diag.nstar.units = {'mmol N m$^{-3}$'};
                 
    % silicate 
    paths.diag.SiO3.file  = {[valiPath,'SiO3/woa18_sio3_clim.nc'],...
                             [valiPath,'SiO3/GLODAPv2.2016b.silicate.nc']};
    paths.diag.SiO3.type  = {'nc','nc'};
    paths.diag.SiO3.var   = {'i_an','silicate'};
    paths.diag.SiO3.zvar  = {'depth','Depth'};
    paths.diag.SiO3.dim   = {'xyzt','xyz'};
    paths.diag.SiO3.lon   = {'lon','lon'};
    paths.diag.SiO3.lat   = {'lat','lat'};
    paths.diag.SiO3.name  = {'WOA-18','GLODAPv2'};
    paths.diag.SiO3.units = {'mmol Si m$^{-3}$','mmol Si m$^{-3}$'};

    % N2O 
    paths.diag.N2O.file   = {[valiPath,'N2O/n2o_NN_format.mat']};
    paths.diag.N2O.type   = {'mat'};
    paths.diag.N2O.var    = {'n2o'};
    paths.diag.N2O.zvar   = {'depth'};
    paths.diag.N2O.dim    = {'xyzt'};
    paths.diag.N2O.lon    = {'lon'};
    paths.diag.N2O.lat    = {'lat'};
    paths.diag.N2O.name   = {'Machine Learning Estimate (Clements et al.)'};
    paths.diag.N2O.units  = {'mmol N$_2$O m$^{-3}$'};
    paths.diag.N2O.factor = {0.001}; %nmol to mmol

    % NO2 
    paths.diag.NO2.file   = {[valiPath,'NO2/no2predRF_05-27-2020.nc']};
    paths.diag.NO2.type   = {'nc'};
    paths.diag.NO2.var    = {'no2'};
    paths.diag.NO2.zvar   = {'depth'};
    paths.diag.NO2.dim    = {'xyzt'};
    paths.diag.NO2.lon    = {'lon'};
    paths.diag.NO2.lat    = {'lat'};
    paths.diag.NO2.name   = {'Machine Learning Estimate (Clements et al.)'};
    paths.diag.NO2.units  = {'mmol N m$^{-3}$'};

    % Chl 
    paths.diag.SFC_CHL.file  = {[valiPath,'SFC_CHL/chl_clim_0p25.nc']};
    paths.diag.SFC_CHL.type  = {'nc'};
    paths.diag.SFC_CHL.dim   = {'xyt'};
    paths.diag.SFC_CHL.var   = {'chl'};
    paths.diag.SFC_CHL.zvar  = {[]};
    paths.diag.SFC_CHL.lon   = {'lon'};
    paths.diag.SFC_CHL.lat   = {'lat'};
    paths.diag.SFC_CHL.name  = {'MODIS-Aqua Chlorophyll'};
    paths.diag.SFC_CHL.units = {'mg C m$^{-3}$'};
        
    % NPP products (VGPM, CBPM, CAFE)    
    paths.diag.NPP.file  = {[valiPath,'NPP/std-VGPMnpp_MODIS_clim2002-2018.nc'],...
                            [valiPath,'NPP/CbPM2npp_MODIS_clim2002-2018.nc'],...
                            [valiPath,'NPP/nppclim_CAFE_MODIS_9km.nc']}; 
    paths.diag.NPP.type  = {'nc','nc','nc'};
    paths.diag.NPP.var   = {'npp','npp','npp'};
    paths.diag.NPP.zvar  = {[],[],[]};
    paths.diag.NPP.dim   = {'xyt','xyt','xyt'};
    paths.diag.NPP.lon   = {'lon','lon','lon'};
    paths.diag.NPP.lat   = {'lat','lat','lat'};
    paths.diag.NPP.name  = {'NPP-VGPM','NPP-CBPM','NPP-CAFE'};
    paths.diag.NPP.units = {'mg C m$^{-2}$ d$^{-1}$','mg C m$^{-2}$ d$^{-1}$','mg C m$^{-2}$ d$^{-1}$'};

    % MLD products (WOCE/NODC/ARGO or ARGO only);
    paths.diag.MLD.file  = {[valiPath,'MLD/Argo_mixedlayers_monthlyclim_12112019.nc'],...
                            [valiPath,'MLD/mld_DR003_c1m_reg2.0.nc']};
    paths.diag.MLD.type  = {'nc','nc'};
    paths.diag.MLD.var   = {'mld_da_mean','mld'};
    paths.diag.MLD.zvar  = {[],[]};
    paths.diag.MLD.dim   = {'txy','xyt'};
    paths.diag.MLD.lon   = {'lon','lon'};
    paths.diag.MLD.lat   = {'lat','lat'};
    paths.diag.MLD.name  = {'Argo Mixed Layer Depth','IFREMER Mixed Layer Depth'}; 
    paths.diag.MLD.units = {'m','m'};

    % POC_FLUX_IN
    paths.diag.POC_FLUX_IN.file   = {[valiPath,'POC_FLUX_IN/clements_100m_flux.mat'],...
                                     [valiPath,'POC_FLUX_IN/Euphotic_Export_2023.nc'],...
                                     [valiPath,'POC_FLUX_IN/TotEZ_Annual_Monthly_July2020.mat'],...
                                     [valiPath,'POC_FLUX_IN/biopump_model_output.nc'],...
                                     [valiPath,'POC_FLUX_IN/Dunne_1deg_DM.mat']};
    paths.diag.POC_FLUX_IN.type   = {'mat','nc','mat','nc','mat'};
    paths.diag.POC_FLUX_IN.var    = {'flux_mean','Flux','avmonTotEZ','POCflux','DunnePOCflux'};
    paths.diag.POC_FLUX_IN.zvar   = {[],[],[],'DEPTH',[]};
    paths.diag.POC_FLUX_IN.dim    = {'xyt','xyt','yxt','yxz','yxt'};
    paths.diag.POC_FLUX_IN.lon    = {'lon','lon','lon','LON','lon'};
    paths.diag.POC_FLUX_IN.lat    = {'lat','lat','lat','LAT','lat'};
    paths.diag.POC_FLUX_IN.name   = {'Clements et al. (2023) 100m POC Flux',...
                                     'Clements et al. (2023) Euphotic POC Flux',...
                                     'Siegel et al. (2014) Euphotic POC Flux',...
                                     'Nowicki et al. (2022) Sinking POC flux',...
                                     'Dunne et al. (2007) Euphotic POC flux' };
    paths.diag.POC_FLUX_IN.units  = {'mmol C m$^{-2}$ s$^{-1}$',...
                                     'mmol C m$^{-2}$ s$^{-1}$',...
                                     'mmol C m$^{-2}$ s$^{-1}$',...
                                     'mmol C m$^{-2}$ s$^{-1}$',...
                                     'mmol C m$^{-2}$ s$^{-1}$'};
    paths.diag.POC_FLUX_IN.factor = {1/(12.01*86400),... % mgC/m2/d to mmolC/m2/s
                                     1/(12.01*86400),... % mgC/m2/d to mmolC/m2/s
                                     1/(12.01*86400),... % mgC/m2/d to mmolC/m2/s 
                                     1/(86400*365.25),...% per year to per s
                                     1/(12.01*86400)};   % mgC/m2/d to mmolC/m2/s


    % Fe data from Tagliabue
    paths.diag.Fe.file   = {[valiPath,'Fe/Monthly_dFe_V2.nc']};
    paths.diag.Fe.type   = {'nc'};
    paths.diag.Fe.var    = {'dFe_RF'};
    paths.diag.Fe.zvar   = {'Depth'};
    paths.diag.Fe.dim    = {'xyzt'};
    paths.diag.Fe.lon    = {'Longitude'};
    paths.diag.Fe.lat    = {'Latitude'};
    paths.diag.Fe.name   = {'Huang et al. (2022)'};
    paths.diag.Fe.units  = {'mmol Fe m$^{-3}$'};
    paths.diag.Fe.factor = {(1/1000)}; % nmol/L to mmol/m3

    % FG_N2O from Simon Yang
    paths.diag.FG_N2O.file   = {[valiPath,'N2O/yang_n2o_flux.mat']};
    paths.diag.FG_N2O.type   = {'mat'};
    paths.diag.FG_N2O.var    = {'dat'};
    paths.diag.FG_N2O.zvar   = {[]};
    paths.diag.FG_N2O.dim    = {'xyt'};
    paths.diag.FG_N2O.lon    = {'lon'};
    paths.diag.FG_N2O.lat    = {'lat'};
    paths.diag.FG_N2O.name   = {'Yang et al. (2020)'};
    paths.diag.FG_N2O.units  = {'mmol N$_2$O m$^{-2}$ s$^{-1}$'};
    paths.diag.FG_N2O.factor = {[(1000*0.5)./(86400*365.25*14)]}; % g N m-2 yr-1 to mmol N2O m-2 s-1 

    % Total Alkalinity via GLODAPv2
    paths.diag.Alk.file   = {[valiPath,'Alk/GLODAPv2.2016b.TAlk.nc']};
    paths.diag.Alk.type   = {'nc'}; 
    paths.diag.Alk.var    = {'TAlk'}; 
    paths.diag.Alk.zvar   = {'Depth'};
    paths.diag.Alk.dim    = {'xyz'};
    paths.diag.Alk.lon    = {'lon'};
    paths.diag.Alk.lat    = {'lat'};
    paths.diag.Alk.name   = {'GLODAPv2'};
    paths.diag.Alk.units  = {'mmol m$^{-3}$'};  
    
    % DIC via GLODAPv2
    paths.diag.DIC.file   = {[valiPath,'DIC/GLODAPv2.2016b.TCO2.nc']};
    paths.diag.DIC.type   = {'nc'}; 
    paths.diag.DIC.var    = {'TCO2'}; 
    paths.diag.DIC.zvar   = {'Depth'};
    paths.diag.DIC.dim    = {'xyz'};
    paths.diag.DIC.lon    = {'lon'};
    paths.diag.DIC.lat    = {'lat'};
    paths.diag.DIC.name   = {'GLODAPv2'};
    paths.diag.DIC.units  = {'mmol m$^{-3}$'};  

    % OMZ thickness
    paths.diag.OMZ.file  = {[valiPath,'OMZ/WOA18_OMZ_thickness.mat'],...
                            [valiPath,'OMZ/Bianchi2012_OMZ_thickness.mat']};
    paths.diag.OMZ.type  = {'mat','mat'};
    paths.diag.OMZ.var   = {'OMZ','OMZ'};
    paths.diag.OMZ.zvar  = {'omzthresh','omzthresh'}; % trick into thinking thresholds are depths
    paths.diag.OMZ.dim   = {'xyz','xyz'};
    paths.diag.OMZ.lon   = {'lon','lon'};
    paths.diag.OMZ.lat   = {'lat','lat'};
    paths.diag.OMZ.name  = {'WOA-18','Bianchi et al. (2012)'};
    paths.diag.OMZ.units = {'m','m'};
end % end method getDiagPaths
