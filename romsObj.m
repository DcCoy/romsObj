%------------------------------------------------------------------------------------------------
classdef romsObj
%-------------------
% Matlab class used to load ROMS data and make various plots,
% comparison figures against validation data, and many
% other actions.
%
% To begin, simply initialize romsObj
% - obj = initROMS(romsObj,sim)
%
% To view available routines
% - methods(obj)
%
% For help
% - help obj.(method)
%------------------
    %----------------------------------------------------------------------------------------
    properties
        info                      % struct containing simulation and variable information (via romsInfo)
        grid                      % struct containing grid data and dimensions (via loadGrid)
        slice                     % struct containing slice coordinates (via sliceROMS/sliceDiag)
        profile                   % struct containing profile coordinates (via getProfile)
        budget                    % struct containing BGC tracer budget results (via getBudg)
        data                      % struct containing ROMS data for comparisons (via loadData, computeVar, sliceROMS)
        diag                      % struct containing validation data for comparions (via loadDiag, sliceDiag)
        paths = getDiagPaths;     % struct containing paths to data, directories, diagnostic data, etc (via initROMS)
        params = getParams;       % struct containing N-cycle parameters
        constants = getConstants; % struct containing constants and conversions
    end % end properties
    %----------------------------------------------------------------------------------------
    methods
        %--------------------------------------------------------------------------------
        function obj = initROMS(obj,simName,varargin)
            % -------------------
            % Initialization method: gathers paths and coordinate variables 
            %
            % Usage: 
            % - obj = initROMS(obj,sim)
            % 
            % Inputs:
            % - simName: ROMS simulation (peru_chile_0p1, peru_chile_0p05, or pacmed_0p25 only)
            %
            % Optional Inputs:
            % - runName:  Simulation name (i.e. spinup, VKV4, specific runs)
            % - domain:   [xi_start xi_end eta_start eta_end] for domain 
            % - coast:    [#] to remove data within #km of coast
            %
            % Example:
            % - obj = initROMS(obj,'peru_chile_0p1','runName','VKV4_tune2')         <-- if obj defined
            % - obj = initROMS(romsObj,'peru_chile_0p1','runName','VKV4_tune2')  <-- if obj undefined
            % -------------------

            %  Begin
            disp(' ');
            disp('---------------------------------');
            disp('---------------------------------');
            disp('            romsObj              ');
            disp('---------------------------------');
            disp('---------------------------------');
            disp(' ');                   
            disp('Initializing simulation paths');

            % Suppress warning messages and addpath to Danny's scripts
            warning off
            addpath /data/project1/demccoy/matlab_scripts/
            addpath /data/project1/demccoy/ROMS/tools/Roms_tools_MF/Preprocessing_tools/

            % Options    
            % Defaults
            default_settings = 0;
            if nargin < 2
                % Load default settings
                default_settings = 1;
                resolution = 10;
                simName    = 'peru_chile_0p1';    
                A.runName  = 'dccoy_VKV4_tune7_fixAx';
                A.gridName = [simName,'_grd.nc'];
                A.domain   = [];
                A.coast    = [];
                A.depth    = [];
                A.bathy    = [];
                A.header   = [];
            end

            % Other run defaults
            if ~default_settings
                if strcmp(simName,'peru_chile_0p1')
                    resolution = 10;
                    A.runName  = 'dccoy_VKV4_tune7_fixAx';
                    A.gridName = [simName,'_grd.nc'];
                    A.domain   = []; %[31 341 301 461];
                    A.depth    = []; %[-750 inf];
                    A.coast    = []; %[20];
                    A.bathy    = [];
                    A.header   = [];
                elseif strcmp(simName,'peru_chile_0p05');
                    resolution = 5;
                    A.runName  = 'dccoy_VKV4_tune7_fixAx';
                    A.gridName = [simName,'_grd.nc']; 
                    A.domain   = [];
                    A.coast    = [];
                    A.depth    = [];
                    A.bathy    = [];
                    A.header   = [];
                elseif strcmp(simName,'pacmed_0p25');
                    resolution = 25;
                    A.runName  = 'pdamien_VKV4_tune9';
                    A.gridName = [simName,'_grd_corrected.nc']; 
                    A.domain   = [];
                    A.coast    = [];
                    A.depth    = [];
                    A.bathy    = [];
                    A.header   = []; % monthly often output as pacmed_avg_*nc
                elseif strcmp(simName,'USSW1');
                    resolution = 1;
                    A.runName  = 'dccoy_files';
                    A.gridName = ['usw1_grd.nc'];
                    A.domain   = [];
                    A.coast    = [];
                    A.depth    = [];
                    A.bathy    = [];
                    A.header   = ['ussw1_'];
                elseif strcmp(simName,'pacmed_0p12');
                    resolution = 12;
                    A.runName  = 'pdamien_VKV4_tune9';
                    A.gridName = [simName,'_grd.nc']; 
                    A.domain   = [];
                    A.coast    = [];
                    A.depth    = [];
                    A.bathy    = [];
                    A.header   = []; % monthly often output as pacmed_avg_*nc            
                end
            end
            % Set default data to all
            A = parse_pv_pairs(A,varargin);

            % Clear any loaded ROMS data
            obj = clearROMS(obj);

            % Set file paths
            obj.paths.simPath   = ['/data/project2/model_output/',simName,'/'];
            obj.paths.runPath   = [obj.paths.simPath,simName,'_',A.runName,'/'];
            obj.paths.config    = ['/data/project2/demccoy/ROMS_configs/',simName,'/'];
            obj.paths.grid      = [obj.paths.config,'grid/',A.gridName];
            obj.paths.phy_avg   = [obj.paths.runPath,'phy/avg/'];            
            obj.paths.phy_his   = [obj.paths.runPath,'phy/his/'];            
            obj.paths.bgc_avg   = [obj.paths.runPath,'bgc/avg/'];
            obj.paths.bgc_his   = [obj.paths.runPath,'bgc/his/'];
            obj.paths.dia_avg   = [obj.paths.runPath,'dia/avg/'];
            obj.paths.dia_his   = [obj.paths.runPath,'dia/his/'];
            obj.paths.flux_avg  = [obj.paths.runPath,'flux/avg/'];
            obj.paths.flux_his  = [obj.paths.runPath,'flux/his/'];

            % initiate directories if they dont exist
            mkdir([obj.paths.runPath,'Figures/']);
            mkdir([obj.paths.runPath,'Figures/Diagnostic']);
            mkdir([obj.paths.runPath,'Figures/Budget']);
            mkdir([obj.paths.runPath,'phy']);
            mkdir([obj.paths.runPath,'bgc']);
            mkdir([obj.paths.runPath,'dia']);
            mkdir([obj.paths.runPath,'phy/avg/']);
            mkdir([obj.paths.runPath,'bgc/avg/']);
            mkdir([obj.paths.runPath,'dia/avg/']);
            mkdir([obj.paths.runPath,'flux/avg/']);
            mkdir([obj.paths.runPath,'phy/his/']);
            mkdir([obj.paths.runPath,'bgc/his/']);
            mkdir([obj.paths.runPath,'dia/his/']);
            mkdir([obj.paths.runPath,'flux/his/']);

            % grab plot paths
            obj.paths.plots.diag    = [obj.paths.runPath,'Figures/Diagnostic/'];
            obj.paths.plots.budget  = [obj.paths.runPath,'Figures/Budget/'];
            obj.paths.plots.tmpfigs = ['/data/project1/demccoy/tmpfigs/'];

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
            obj.info.resolution = resolution;

            % Get info for files
            obj = romsInfo(obj);
            % Load grid
            obj = loadGrid(obj);
        end % end method initROMS

        %--------------------------------------------------------------------------------
        function [obj] = romsInfo(obj)
            % ----------------------
            % Obtains info from sample roms files in run directory
            %
            % Usage:
            % - [obj] = romsInfo(obj) 
            % ----------------------
            disp('Loading file meta data (ncinfo) into obj.info');
            rmpath /data/project1/demccoy/ROMS/ROMS_tools/nc_tools/

            % List file_types
            file_types  = {'phy','bgc','dia','flux'};
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
                        disp(['No ',[file_types{i},'_',file_dirs{j}],' files found...skipping']);
                        obj.info.([file_types{i},'_',file_dirs{j}]) = [];
                        continue
                    else
                        good_dir = 1;
                    end

                    % Get sample file (for ncinfo)
                    file_path = obj.paths.([file_types{i},'_',file_dirs{j}]);

                    % Grab ncinfo for each file
                    for k = 1:length(tmpfiles) 
                        % Call ncinfo
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
                        
                        % Save into struct
                        obj.info.([file_types{i},'_',file_dirs{j}])(k) = tmp;
                    end
                end
            end

            % Kill object if no files are found
            if good_dir == 0
                disp('ERROR: No ROMS files found, check initROMS inputs or fix output directories');
                FAIL
            end
        end % end method romsInfo

        %--------------------------------------------------------------------------------
        function [obj] = loadGrid(obj)
            % ----------------------
            % Loads 2D grid information into obj.grid
            %
            %
            % Usage:
            % - obj = loadGrid(obj,varargin)
            %
            % Optional inputs (varargin):
            % - full = if (1), loads all grid info (area, etc).
            %           defaults to (0) and only loads basic grid fields 
            %
            % Example:
            % - obj = loadGrid(obj);
            % ----------------------

            % ----------------------
            disp('Loading grid data into obj.grid');
            rmpath('/data/project1/demccoy/ROMS/ROMS_tools/nc_tools/');
            addpath('/usr/local/MATLAB/R2019b/toolbox/matlab/imagesci/');

            % Requires romsInfo
            try; obj.info;
            catch; obj = romsInfo(obj); 
            end

            % Grab z_avg_deps (WOA-18 depths, useful for diagnostics/validation)
            tmp = load('/data/project1/demccoy/ROMS/validation/wcoord.mat');
            obj.grid.z_avg_dep = tmp.wcoord0p25.depth;

            % Extract dimensions from file
            file_types = {'phy_avg','phy_his','bgc_avg','bgc_his','dia_avg','dia_his'};
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

            % Load u/v fields
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

            % Make coord data all single
            obj.grid = romsObj.struct2double(obj.grid);

            % (OPTIONAL)
            if obj.grid.region.coast_lim>(-inf)
                % Calculate distance-from-coast (m)
                if isfield(obj.grid,'coastdist') == 0 
                    disp('Refining mask via Dist2Coast');
                    obj = Dist2Coast(obj);
                end
                % Update mask
                obj.grid.mask_rho(obj.grid.coastdist < obj.grid.region.coast_lim*1000) = NaN;
            end
            if ~isempty(obj.grid.region.bathy_lim)
                % Update mask to exclude shallow points
                disp('Refining mask to exclude shallow points');
                obj.grid.mask_rho(obj.grid.h <= obj.grid.region.bathy_lim) = NaN;
            end
        end % end method loadGrid

        %--------------------------------------------------------------------------------
        function [obj] = loadDepth(obj,file,varargin) 
            % ----------------------
            % Load Z-grid info (z_r,z_w), which are output-dependent.
            % Also updates 3D mask for raw data.
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
            disp('Grabbing z_r, z_w, Hz');

            % Call loadGrid
            try; obj.grid.h;
            catch; obj = loadGrid(obj);
            end
            
            % Process optional inputs
            A.type = 'avg';
            if isempty(obj.info.phy_avg)
                A.type = 'his';
                if isempty(obj.info.phy_his)
                    disp('No physical output found!');
                    kill
                end
            end
            A.full = 0;
            A.mean = 0;
            A = parse_pv_pairs(A,varargin);
            file_type = ['phy_',A.type];

            % Grab files
            files = {obj.info.(file_type).Filename};

            % Load vertical coordinates (and zlevs4 info)
            atts = {obj.info.(file_type)(1).Attributes.Name};
            atts_to_read = {'theta_s','theta_b','hc'};
            for i = 1:length(atts_to_read);
                ind = find(strcmp(atts_to_read{i},atts)==1);
                obj.grid.(atts_to_read{i}) = obj.info.(file_type)(file(1)).Attributes(ind).Value;
            end

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
                if length(files)>1
                    fprintf([num2str(ff),'..']);
                    if ff == length(files)
                        fprintf(['\n']);
                    end
                end
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
                    z_r = zlevs4(tmp.h,tmp.zeta(:,:,i),obj.grid.theta_s,obj.grid.theta_b,obj.grid.hc,obj.grid.s_rho,'r','new2012');
                    z_w = zlevs4(tmp.h,tmp.zeta(:,:,i),obj.grid.theta_s,obj.grid.theta_b,obj.grid.hc,obj.grid.s_rho,'w','new2012');
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
            if A.mean
                obj.grid.(A.type).Hz  = (out.Hz ./ ff);
                obj.grid.(A.type).z_r = (out.z_r ./ ff);
                if A.full
                    obj.grid.(A.type).z_u = (out.z_u ./ ff);
                    obj.grid.(A.type).z_v = (out.z_v ./ ff);
                    obj.grid.(A.type).z_w = (out.z_w ./ ff);
                end
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

            % Apply depth limits to mask (optional)
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

            % Defaults for optional  arguments
            A.type   = ['avg']; % input file type
            A.mean   = [0];
            A = parse_pv_pairs(A,varargin);

            % Restrict file types to 'avg' or 'his'
            file_types = {'phy_avg','phy_his',...
                          'bgc_avg','bgc_his',...
                          'dia_avg','dia_his',...
                          'flux_avg','flux_his'};
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
                        disp(' ');
                        disp('vars must be a cell array');
                        disp(' ');
                        return
                    end
                end
            end

            % Load each variable into structure
            for i = 1:length(vars)
                disp(['Loading ',vars{i},' into struct']);
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
                    % Call computeVar
                    try; obj = computeVar(obj,vars(i),file,'type',A.type); return
                    catch; disp([vars{i},' is not a variable, and isnt available via computeVar']);
                    end
                end

                % Check for dimensions, grab index
                idx      = find(strcmp(vars{i},{obj.info.(file_types{file_dir})(file(1)).Variables.Name})==1);
                tmp_size = squeeze(obj.info.(file_types{file_dir})(file(1)).Variables(idx).Size); 
                if length(tmp_size)==4
                    ind = obj.grid.rho4D;    
                    if strcmp(vars{i},'u');
                        ind = obj.grid.u4D;
                    elseif strcmp(vars{i},'v');
                        ind = obj.grid.v4D;
                    end
                elseif length(tmp_size)==3
                    ind = obj.grid.rho3D;
                    if strcmp(vars{i},'u');
                        ind = obj.grid.u3D;
                    elseif strcmp(vars{i},'v');
                        ind = obj.grid.v3D;
                    end
                elseif length(tmp_size)==2
                    ind = obj.grid.rho2D;
                    if strcmp(vars{i},'u');
                        ind = obj.grid.u2D;
                    elseif strcmp(vars{i},'v');
                        ind = obj.grid.v2D;
                    end
                else
                    ind = [1 inf];
                end

                % Grab variable info from output file
                obj.data.(A.type).(vars{i}).name  = obj.info.(file_types{file_dir})(file(1)).Variables(idx).Attributes(1).Value;
                obj.data.(A.type).(vars{i}).units = obj.info.(file_types{file_dir})(file(1)).Variables(idx).Attributes(2).Value;
                for ff = 1:length(file)
                    if length(file)>1
                        fprintf([num2str(ff),'..']);
                        if ff == length(file)
                            fprintf(['\n']);
                        end
                    end
                    % Load data
                    tmp.data = ncread([obj.info.(file_types{file_dir})(file(ff)).Filename],vars{i},[ind(1,:)],[ind(2,:)]);
                    tmp.dims = {obj.info.(file_types{file_dir})(file(1)).Variables(idx).Dimensions.Name}; 
                    % Check for averaging switch
                    if A.mean
                        if ff == 1
                            out.data = [];
                        end
                        if length(ind)>=2
                            % Apply mask
                            if strcmp(vars{i},'u');
                                out.data = [out.data + (tmp.data .* obj.grid.mask_u)];
                            elseif strcmp(vars{i},'v');
                                out.data = [out.data + (tmp.data .* obj.grid.mask_v)];
                            else
                                out.data = [out.data + (tmp.data .* obj.grid.mask_rho)];
                            end
                        else
                            out.data = [out.data + tmp.data];
                        end
                    else    
                        if length(ind)>=2
                            % Apply mask
                            if strcmp(vars{i},'u');
                                out.data{ff} = tmp.data .* obj.grid.mask_u;
                            elseif strcmp(vars{i},'v');
                                out.data{ff} = tmp.data .* obj.grid.mask_v;
                            else
                                out.data{ff} = tmp.data .* obj.grid.mask_rho;
                            end
                        else
                            out.data{ff} = tmp.data;
                        end
                    end
                end
                if A.mean
                    obj.data.(A.type).(vars{i}).data = (out.data ./ ff);
                else
                    obj.data.(A.type).(vars{i}).data = cat(length(ind),out.data{:});
                end
                if strcmp(vars{i},'rho');
                    % Remove 'anomaly'
                    disp('NOTE: Correcting density anomaly --> density');
                    obj.data.(A.type).rho.data = obj.data.(A.type).rho.data + 27.5;
                    obj.data.(A.type).rho.name = 'averaged density';
                end
                obj.data.(A.type).(vars{i}).dims = tmp.dims;
                clear out

                % Replace Unit Strings with Latex version
                strings_to_replace = {'kilogram meter-3','kg m$^{-3}$';...
                                      'meter second-1'  ,'m s$^{-1}$';...
                                      'meter2 second-1' ,'m$^{2}$ s$^{-1}$';...
                                      'meter'           ,'m';...
                                      'Celsius'         ,'$^{o}$C';...
                                      'mMol m-3'        ,'mmol m$^{-3}$';...
                                      'mMol P m-3'      ,'mmol P m$^{-3}$';...
                                      'mMol C m-3'      ,'mmol C m$^{-3}$';...
                                      'mMol N m-3'      ,'mmol N m$^{-3}$';...
                                      'mMol Fe m-3'     ,'mmol Fe m$^{-3}$';...
                                      'mMol O2 m-3'     ,'mmol O$_2$ m$^{-3}$';...
                                      'mMol N2O m-3'    ,'mmol N$_2$O m$^{-3}$';...
                                      'mMol N2 m-3'     ,'mmol N$_2$ m$^{-3}$';...
                                      'mMol CaCO3 m-3'  ,'mmol CaCO$_3$';...
                                      'W m-2'           ,'W m$^{-2}$';...
                                      'mmol N/m3/s'     ,'mmol N m$^{-3}$ s$^{-1}$';...
                                      'mmol N2O/m3/s'   ,'mmol N$_2$O m$^{-3}$ s$^{-1}$';...
                                      'mmol N2/m3/s'    ,'mmol N$_2$ m$^{-3}$ s$^{-1}$';...
                                      'mmol/m3/s'       ,'mmol m$^{-3}$ s$^{-1}$';...
                                      'mmol/m2/s'       ,'mmol m$^{-2}$ s$^{-1}$'};

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
        function obj = v_interp(obj,vars,z_r,deps)
            % -------------------
            % Method to vertically interpolate loaded data 
            %
            % Usage:
            % - obj = v_interp(obj,vars,file,varargin)
            %
            % Inputs:
            % - vars   = variable(s) to load, as a cell array
			% - z_r    = depth or rho points, via load depth
			% - deps   = 1xnz array of depths to interpolate to
            %
            % Examples:
            % - obj = loadData(obj,{'temp','salt'},1);
			% - obj = loadDepth(obj,1);
			% - obj = v_interp(obj,{'temp','salt'},obj.grid.avg.z_r,[5:10:1000]);
            % -------------------		

			% Get size of rho depths (z_r)
			if length(size(z_r))==4
				[Lp,Mp,N,T] = size(z_r);
			elseif length(size(z_r))==3
				[Lp,Mp,N] = size(z_r);
			end

			% Find the grid position of the nearest vertical levels
			a=z<depth;
			levs=squeeze(sum(a,1));
			levs(levs==N)=N-1;
			warning off
			mask=levs./levs;
			warning on
			[imat,jmat]=meshgrid((1:Lp),(1:Mp));
			pos=N*Mp*(imat-1)+N*(jmat-1)+levs;
			pos(isnan(mask))=1;
			%
			% Do the interpolation
			%
			z1=z(pos+1);
			z2=z(pos);
			v1=var(pos+1);
			v2=var(pos);
			vnew=mask.*(((v1-v2)*depth+v2.*z1-v1.*z2)./(z1-z2));

		end % end method vinterp

        %--------------------------------------------------------------------------------
        function obj = computeVar(obj,vars,file,varargin)
            % ------------------
            % Computes additional fields like spiciness, potential density etc.
            % 
            % List of available variables:
            % - pres     (pressure)
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
            %
            % Example:
            % - obj = computeVar(obj,{'MLD'},1,'type','avg');
            % ------------------

            % process inputs
            A.type    = 'avg';
            A.ip      = [];
            A.dep     = [];
            A.thresh  = 20; %uM O2 for OMZ thickness
            A.zlim    = inf;
            A         = parse_pv_pairs(A,varargin);

            % dims for vertical computations
            nx = obj.grid.nx;
            ny = obj.grid.ny;
            nz = obj.grid.nz;

            % process requests
            for i = 1:length(vars)
                % pressure calc
                if strcmp(vars{i},'pres') % & ~isfield(obj.data,'pres');
                    disp('Calculating sea pressure');
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
                % bvf calc
                elseif strcmp(vars{i},'bvf') | strcmp(vars{i},'pv');
                    disp('Calculating bvf (Brunt Vaisala)');
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
                    obj = loadData(obj,{'temp','salt'},file,'type',A.type);
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
                    obj = loadData(obj,{'NO3','NO2','PO4'},file,'type',A.type);
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
                    obj     = loadData(obj,{'TOT_PROD'},file,'type',A.type);
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
                    obj = loadData(obj,{'zeta'},file,'type',A.type);
                    obj = loadDiag(obj,{'SSH'},0);
                    slacorr = nanmedian(obj.diag.SSH.slice(:)) - nanmedian(obj.data.(A.type).zeta.data(:));
                    disp(' '); disp(['Adding correction of ',num2str(slacorr),'m to ROMS SSH']);
                    obj.data.(A.type).SSH.data  = obj.data.(A.type).zeta.data + slacorr;
                    obj.data.(A.type).SSH.name  = 'averaged sea-surface height';
                    obj.data.(A.type).SSH.units = 'm'; 
                    obj.data.(A.type).SSH.dims  = obj.data.(A.type).zeta.dims; 
                % Wind Stress or Wind Stress Curl
                elseif strcmp(vars{i},'WS') | strcmp(vars{i},'ws') | strcmp(vars{i},'WSC') | strcmp(vars{i},'wsc');
                    obj    = loadData(obj,{'sustr','svstr'},file,'type',A.type);
                    tmpu   = obj.data.(A.type).sustr.data;
                    tmpv   = obj.data.(A.type).svstr.data;
                    tmpang = obj.grid.angle;
                    [tmpws,tmpwsc] = romsObj.WindStress(tmpu,tmpv,obj.grid.lon_rho,obj.grid.lat_rho,tmpang);
                    obj.data.(A.type).ws.data   = tmpws;
                    obj.data.(A.type).wsc.data  = tmpwsc;
                    obj.data.(A.type).ws.data   = obj.data.(A.type).ws.data .* obj.grid.mask_rho;
                    obj.data.(A.type).wsc.data  = obj.data.(A.type).wsc.data .* obj.grid.mask_rho;
                    obj.data.(A.type).ws.name   = 'averaged wind-stress';
                    obj.data.(A.type).wsc.name  = 'averaged wind-stress curl';
                    obj.data.(A.type).ws.units  = 'N m$^{-2}$';
                    obj.data.(A.type).wsc.units = 'N m$^{-2}$';
                    obj.data.(A.type).ws.dims   = obj.data.(A.type).sustr.dims;
                    obj.data.(A.type).wsc.dims  = obj.data.(A.type).sustr.dims;
                % Mixed-layer depth (MLD)
                elseif strcmp(vars{i},'MLD') | strcmp(vars{i},'mld');
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
                    if isempty(A.ip) & isempty(A.dep)
                        disp('Supply a depth input via ''ip'' or ''dep''');
                        return
                    elseif ~isempty(A.ip);
                        obj = ipslice(obj,{'u','v'},A.ip,file);
                    elseif ~isempty(A.dep);
                        if ~ismember(A.dep,obj.grid.z_avg_dep);
                            disp('Choose a depth from obj.grid.z_avg_dep only');
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
                % Jden_N2O
                elseif strcmp(vars{i},'Jden_N2O');
                    obj = loadData(obj,{'DENITRIF2','N2OSODEN_CONS'},file,'type',A.type);
                    obj.data.(A.type).Jden_N2O.data  = (obj.data.(A.type).DENITRIF2.data ./ 2) - obj.data.(A.type).N2OSODEN_CONS.data; 
                    obj.data.(A.type).Jden_N2O.name  = 'Net $N_2O$ production from denitrification';
                    obj.data.(A.type).Jden_N2O.units = 'mmol N$_2$O m$^{-3}$ s$^{-1}$';
                    obj.data.(A.type).Jden_N2O.dims  = obj.data.(A.type).DENITRIF2.dims;
                % Jnit_N2O
                elseif strcmp(vars{i},'Jnit_N2O');
                    obj = loadData(obj,{'N2OAMMOX'},file,'type',A.type);
                    obj.data.(A.type).Jnit_N2O.data  = obj.data.(A.type).N2OAMMOX.data - obj.data.(A.type).N2OAO1_CONS.data;    
                    obj.data.(A.type).Jnit_N2O.name  = 'Net $N_2O$ production from nitrification';    
                    obj.data.(A.type).Jnit_N2O.units = 'mmol N$_2$O m$^{-3}$ s$^{-1}$';
                    obj.data.(A.type).Jnit_N2O.dims  = obj.data.(A.type).N2OAMMOX.dims;
                % J_N2O
                elseif strcmp(vars{i},'J_N2O');
                    obj = loadData(obj,{'N2OAMMOX','DENITRIF2','DENITRIF3'},file,'type',A.type);
                    obj.data.(A.type).J_N2O.data  = (obj.data.(A.type).DENITRIF2.data ./ 2) + obj.data.(A.type).N2OAMMOX.data - ...
                                                     obj.data.(A.type).DENITRIF3.data;
                    obj.data.(A.type).J_N2O.name  = 'N$_2$O sources-minus-sinks';
                    obj.data.(A.type).J_N2O.units = 'mmol N$_2$O m$^{-3}$ s$^{-1}$';
                    obj.data.(A.type).J_N2O.dims  = obj.data.(A.type).DENITRIF2.dims;
                % Apparent Oxygen Utiliziation (AOU)
                elseif strcmp(vars{i},'AOU');
                    obj = loadData(obj,{'temp','salt','O2'},file,'type',A.type);
                    o2_sat = romsObj.o2_sat(obj.data.(A.type).temp.data,obj.data.(A.type).salt.data);
                    obj.data.(A.type).AOU.data  = o2_sat - obj.data.(A.type).O2.data;
                    obj.data.(A.type).AOU.name  = 'averaged AOU';
                    obj.data.(A.type).AOU.units = 'mmol O$_2 m$^{-3}$';
                    obj.data.(A.type).AOU.dims  = obj.data.(A.type).O2.dims;
                % Sat/Delta N2O
                elseif strcmp(vars{i},'DeltaN2O');
                    obj = loadData(obj,{'temp','salt','N2O'},file,'type',A.type);
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
                    obj = loadDepth(obj,file);
                    obj = loadData(obj,{'O2'},file,'type',A.type);
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
                    obj = loadData(obj,{'NH4','NO2'},file,'type',A.type);
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
                    obj = loadData(obj,{'NO3','NO2'},file,'type',A.type);
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
                % N2O_BRY (N2O_ATM + N2O_BOU)
                elseif strcmp(vars{i},'N2O_BRY');
                    obj = loadData(obj,{'N2O_SIDEN','N2O_ATM'},file,'type',A.type);
                    tmpbou = obj.data.(A.type).N2O_SIDEN.data;
                    tmpatm = obj.data.(A.type).N2O_ATM.data;
                    tmpbou(tmpbou<0) = 0;
                    tmpatm(tmpatm<0) = 0;
                    obj.data.(A.type).N2O_BRY.data  = tmpbou + tmpatm;
                    obj.data.(A.type).N2O_BRY.name  = 'averaged Nitrous oxide from boundaries';
                    obj.data.(A.type).N2O_BRY.units = obj.data.(A.type).N2O_ATM.units;
                    obj.data.(A.type).N2O_BRY.dims  = obj.data.(A.type).N2O_ATM.dims;
                % NO2AMMOX (AMMOX - 2*N2OAMMOX)
                elseif strcmp(vars{i},'NO2AMMOX');
                    obj = loadData(obj,{'AMMOX','N2OAMMOX'},file,'type',A.type);
                    tmpammox    = obj.data.(A.type).AMMOX.data;
                    tmpn2oammox = obj.data.(A.type).N2OAMMOX.data;
                    tmpammox(tmpammox<0) = 0;
                    tmpn2oammox(tmpn2oammox<0) = 0;
                    obj.data.(A.type).NO2AMMOX.data  = tmpammox - 2.*(tmpn2oammox);
                    obj.data.(A.type).NO2AMMOX.name  = 'averaged NH$^{+}_4$ oxidation to NO$^{-}_2$';
                    obj.data.(A.type).NO2AMMOX.units = obj.data.(A.type).AMMOX.units;
                    obj.data.(A.type).NO2AMMOX.dims  = obj.data.(A.type).AMMOX.dims;
                else
                    disp([vars{i},' is already, or cant be, calculated']);
                end
            end
        end % end method computeVar

        %--------------------------------------------------------------------------------
        function obj = sliceROMS(obj,vars,choice,deg,file,varargin);
            % -------------------
            % Takes 2D depth slice of ROMS along a given latitude or longitude
            % 
            % Usage:
            % - obj = sliceROMS(obj,vars,choice,deg,file,varargin);
            % 
            % Inputs:
            % - vars   = ROMS variable(s) to slice, as a cell array
            % - choice = 'lon','lat','xi','eta' 
            %              (lat/lon slices along a given lat/lon degree)
            %              (xi/eta slices along a given xi or eta index, use gridView(obj))
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

            % Clear slice struct
            obj.slice = [];
    
            % Check inputs
            if nargin<5
                disp('Incorrect number of inputs');
                help sliceROMS
                return
            end    

            % Grab user inputs
            A.type = ['avg'];
            A.zdep = [];

            % If calling lat/lon, force zslice
            if strcmp(choice,'lat') | strcmp(choice,'lon');
                A.zdep = obj.grid.z_avg_dep;    
            end
            A = parse_pv_pairs(A,varargin);

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
                try
                    obj = loadData(obj,vars(i),file,'type',A.type);
                catch
                    obj = computeVar(obj,vars(i),file,'type',A.type);
                end
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
                            disp('Slice outside ROMS domain (check latitude)')
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
                                outdat(i,:) = interp1(-z_raw,dat_raw,A.zdep);
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
                                if isempty(find(isnan(z_raw)==0)) | isempty(find(isnan(dat_raw)==0))
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
                            disp('Slice outside ROMS domain (check longitude)')
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
                                outdat(i,:) = interp1(-z_raw,dat_raw,A.zdep);
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
                                outdat(i,:,j) = interp1(-z_raw,dat_raw,A.zdep);
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
            disp('Grabbing zsliced variables');
            
            % Process optional inputs
            A.type = 'avg';
            A.clear = 'on';
            A = parse_pv_pairs(A,varargin);

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
        function obj = zslice_data(obj,vars,zdep); 
            % -----------------------
            % Shortcut to zslice loaded data
            %
            % Usage:
            % - obj = zslice_data(obj,vars,zdep);
            %
            % Inputs:
            % - obj  = romsObject
            % - vars = variable to slice 
            % - zdep = depth to slice
            %
            % Steps:
            %    Load z_r (or z_u, z_v) into obj.grid.avg.z_r;
            %    Load data into obj.grid.avg.(vars).data      
            % -----------------------

            % Grab nt
            dims = size(obj.data.avg.(vars{1}).data);
            if length(dims) == 3
                nt = 1;
            elseif length(dims)>4
                nt = dims(end);
            end

            % Interpolate variable to zdep level
            for i = 1:length(vars)
                for t = 1:nt
                    vnew   = nan(obj.grid.nx,obj.grid.ny,length(zdep));
                    tmpvar = obj.data.avg.(vars{i}).data(:,:,:,t);
                    if strcmp(vars{i},'u');
                        tmpz = -obj.grid.avg.z_u(:,:,:,t);
                        vnew = vnew(1:end-1,:,:);
                    elseif strcmp(vars{i},'v');
                        tmpz = -obj.grid.avg.z_v(:,:,:,t);
                        vnew = vnew(:,1:end-1,:);
                    else
                        tmpz = -obj.grid.avg.z_r(:,:,:,t);
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
                    obj.data.avg.(vars{i}).slice(:,:,:,t)  = squeeze(vnew);
                end
            end
        end % end method zslice_data

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
            disp('Grabbing ipsliced variables');

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

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
        function obj = getBudg(obj,vars,file,varargin)
            % --------------------
            % Main method to perform budget analysis on a ROMS tracer (vars).
            %
            % Usage:
            % - obj = getBudg(obj,vars,file,varargin)
            %
            % Inputs:
            % - vars = cell array of budgets to close 
            % - file = (#s) load specific time-averages 
            %        = 0 load all files and average (for monthly)
            %
            % Optional Inputs:
            % - hisfile  = alternate file#s for history files (in case output is limited) 
            % - physfile = alternate file#s for flux files (in case output is limited) 
            % - srcsnk   = (1) to call sources_sinks (default 0)
            % - int      = (1) to call intBudg (default 0)
            % - clean    = (0) to keep all data loaded (default 1)
            % - fix      = (1) to adjust x/y advective fluxes based on net
            %               For some simulations, the coast has large budget imbalances
            %               due to issues with x (main issue) and y advective fluxes
            %
            % Example:
            % - obj = getBudg(obj,'N2O',1)
            %
            %   This will instruct budget methods on which
            %   tracer to perform budget on. Here, tracking N2O.
            %
            % Available tracer budgets:
            % NO3
            % NO2
            % NH4
            % N2O
            % N2
            % O2
            % --------------------

            % Toggles
            A.int      = [0];
            A.srcsnk   = [0];
            A.clean    = [1];
            A.hisfile  = [];
            A.physfile = [];
            A.fix      = [0];
            A          = parse_pv_pairs(A,varargin);

            % Check inputs
            if isempty(vars)
                disp('vars must be defined, see romsObj.getBudg')
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

            % Constants (these may change!!! check ROMS code)
            param.DONrefract = 0.0115;   % Fraction of DON to refractory pool
            param.Q          = 0.137;    % N/C ratio

            % Get budget
            for i = 1:length(vars)
                disp(['Computing budget for ',vars{i}]);
                if strcmp(vars{i},'NO3')
                    obj.budget.(vars{i}).info.tits   = {'NO$^{-}_3$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'NITROX','DENITRIF1','SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE'};
                    obj.budget.(vars{i}).info.rtits  = {'NO$^{-}_2$ oxidation','NO$^{-}_3$ reduction',...
                                                        'NO$^{-}_3$ uptake via small phyto','NO$^{-}_3$ uptake via diatoms',...
                                                        'NO$^{-}_3$ uptake via diazotrophs'};
                    obj.budget.(vars{i}).info.smseq  = [(1) (-1) (-1) (-1) (-1)];
                    obj.budget.(vars{i}).info.fluxes = {'SED_DENITRIF'};
                    obj.budget.(vars{i}).info.lvls   = {'sed'};
                    obj.budget.(vars{i}).info.ftits  = {'Sediment denitrification'};
                elseif strcmp(vars{i},'NO2')
                    obj.budget.(vars{i}).info.tits   = {'NO$^{-}_2$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'AMMOX','N2OAMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2',...
                                                        'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
                    obj.budget.(vars{i}).info.rtits  = {'NH$^{+}_4$ oxidation','NH$^{+}_4$ oxidation to N$_2$O',...
                                                        'NO$^{-}_2$ oxidation','Anammox','NO$^{-}_3$ reduction',...
                                                        'NO$^{-}_2$ reduction',...
                                                        'NO$^{-}_2$ uptake via small phyto','NO$^{-}_2$ uptake via diatoms',...
                                                        'NO$^{-}_2$ uptake via diazotrophs'};
                    obj.budget.(vars{i}).info.smseq  = [(1) (-2) (-1) (-1) (1) (-1) (-1) (-1) (-1)];    
                    obj.budget.(vars{i}).info.fluxes = {[]};
                    obj.budget.(vars{i}).info.lvls   = {[]};
                    obj.budget.(vars{i}).info.ftits  = {[]};
                elseif strcmp(vars{i},'NH4')
                    obj.budget.(vars{i}).info.tits   = {'NH$^{+}_4$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'SP_NH4_UPTAKE','DIAT_NH4_UPTAKE','DIAZ_NH4_UPTAKE','AMMOX','ANAMMOX',...
                                                        'DON_REMIN','DONr_REMIN','ZOO_LOSS_DIC','SP_LOSS_DIC','DIAT_LOSS_DIC',...
                                                        'DIAZ_LOSS_DIC','SP_GRAZE_DIC','DIAT_GRAZE_DIC','DIAZ_GRAZE_DIC',...
                                                        'POC_REMIN'}; 
                    obj.budget.(vars{i}).info.rtits  = {'NH$^{+}_4$ uptake via small phyto','NH$^{+}_4$ uptake via diatoms',...
                                                        'NH$^{+}_4$ uptake via diazotrophs','NH$^{+}_4$ oxidation',...
                                                        'Anammox','DON remineralization','Refractory DON remineralization',...
                                                        'Zooplankton mortality','Small phyto mortality',...
                                                        'Diatom mortality','Diazotroph mortality','Small phyto grazing loss',...
                                                        'Diatom grazing loss','Diatotroph grazing loss','POC remineralization'};
                    %obj.budget.(vars{i}).info.smseq  = [(-1) (-1) (-1) (-1) (-1) (1) (1) (param.Q) (param.Q) (param.Q) (param.Q) ...
                    %                                    (param.Q) (param.Q) (param.Q) (param.Q.*(1-param.DONrefract))];
                    obj.budget.(vars{i}).info.smseq  = [(-1) (-1) (-1) (-1) (-1) (1) (1) (param.Q) (param.Q) (param.Q) (param.Q + (param.Q./2.23)) ...
                                                        (param.Q) (param.Q) (param.Q + (param.Q./2.23)) (param.Q.*(1-param.DONrefract))];
                    obj.budget.(vars{i}).info.fluxes = {[]};
                    obj.budget.(vars{i}).info.lvls   = {[]};    
                    obj.budget.(vars{i}).info.ftits  = {[]};
                elseif strcmp(vars{i},'N2')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'DENITRIF3','ANAMMOX'};
                    obj.budget.(vars{i}).info.rtits  = {'N$_2$O reduction','Anammox'};
                    obj.budget.(vars{i}).info.smseq  = [(1) (1)]; 
                    %obj.budget.(vars{i}).info.fluxes = {'FG_N2','N2_SED'};
                    %obj.budget.(vars{i}).info.lvls   = {'sfc',  'sed'};
                    %obj.budget.(vars{i}).info.ftits  = {'$\Phi$ N$_2$','N$_2$ from sediment'};
                    obj.budget.(vars{i}).info.fluxes = {[]};
                    obj.budget.(vars{i}).info.lvls   = {[]};
                    obj.budget.(vars{i}).info.ftits  = {[]};
                elseif strcmp(vars{i},'N2O')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$O'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
                    obj.budget.(vars{i}).info.rtits  = {'NO$^{-}_2$ reduction','NH$^{+}_4$ oxidation to N$_2$O','N$_2$O reduction'};
                    obj.budget.(vars{i}).info.smseq  = [(0.5) (1) (-1)]; 
                    obj.budget.(vars{i}).info.fluxes = {'FG_N2O'};
                    obj.budget.(vars{i}).info.lvls   = {   'sfc'};
                    obj.budget.(vars{i}).info.ftits  = {'$\Phi$ N$_2$O'};
                elseif strcmp(vars{i},'N2O_SODEN')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$O$_{den}$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'DENITRIF2','N2OSODEN_CONS'};
                    obj.budget.(vars{i}).info.smseq  = [(0.5) (-1)]; 
                    obj.budget.(vars{i}).info.fluxes = {'FG_N2O_SODEN'};
                    obj.budget.(vars{i}).info.lvls   = {   'sfc'};
                elseif strcmp(vars{i},'N2O_AO1')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$O$_{nit}$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'N2OAMMOX','N2OAO1_CONS'};
                    obj.budget.(vars{i}).info.smseq  = [(1) (-1)]; 
                    obj.budget.(vars{i}).info.fluxes = {'FG_N2O_AO1'};
                    obj.budget.(vars{i}).info.lvls   = {   'sfc'};        
                elseif strcmp(vars{i},'N2O_ATM')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$O$_{atm}$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'N2OATM_CONS'};
                    obj.budget.(vars{i}).info.smseq  = [(-1)]; 
                    obj.budget.(vars{i}).info.fluxes = {'FG_N2O_ATM'};
                    obj.budget.(vars{i}).info.lvls   = {   'sfc'};
                elseif strcmp(vars{i},'N2O_SIDEN')
                    obj.budget.(vars{i}).info.tits   = {'N$_2$O_{bou}$'};
                    obj.budget.(vars{i}).info.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
                    obj.budget.(vars{i}).info.rates  = {'N2OSIDEN_CONS'};
                    obj.budget.(vars{i}).info.smseq  = [(-1)]; 
                    obj.budget.(vars{i}).info.fluxes = {'FG_N2O_SIDEN'};
                    obj.budget.(vars{i}).info.lvls   = {   'sfc'};
                elseif strcmp(vars{i},'O2');
                    obj.budget.(vars{i}).info.tits   = {'O$_2$'};
                    obj.budget.(vars{i}).info.units  = {'mmol O$_2$ m$^{-3}$ s$^{-1}$','mmol O$_2$ m$^{-2}$ s$^{-1}$','mmol O$_2$ s$^{-1}$'};
                    obj.budget.(vars{i}).info.rates  = {'O2_PRODUCTION','O2_CONSUMPTION'};
                    obj.budget.(vars{i}).info.rtits  = {'O$_2$ production','O$_2$ consumption'};
                    obj.budget.(vars{i}).info.smseq  = [(1) (-1)];
                    %obj.budget.(vars{i}).info.fluxes = ['FG_O2']; 
                    %obj.budget.(vars{i}).info.lvls   = {'sfc'};
                    %obj.budget.(vars{i}).info.ftits  = {'$\Phi$ O$_2$'};
                    obj.budget.(vars{i}).info.fluxes = {[]}; 
                    obj.budget.(vars{i}).info.lvls   = {[]};
                    obj.budget.(vars{i}).info.ftits  = {[]};
                end
            
                % Load all 3D output terms
                terms = [vars{i},obj.budget.(vars{i}).info.rates,obj.budget.(vars{i}).info.fluxes];
                terms = [terms(find(~cellfun(@isempty,terms)))];
                obj   = loadDepth(obj,file,'full',1);
                obj   = loadData(obj,terms,file,'type','avg');

                % Integrate 3D variables, rates vertically
                terms = [vars{i},obj.budget.(vars{i}).info.rates];
                terms = [terms(find(~cellfun(@isempty,terms)))];
                obj   = intVar(obj,terms);

                % Load and process 2D fluxes
                obj = getFluxes(obj,vars(i),obj.budget.(vars{i}).info.fluxes,obj.budget.(vars{i}).info.lvls);

                % Get dCdt, advection, sms, net terms
                obj = computeDcDt(obj,vars(i),hisfile);
                obj = computeXYZflux(obj,vars(i),physfile);
                obj = computeSMS(obj,vars(i),obj.budget.(vars{i}).info.rates,...
                    obj.budget.(vars{i}).info.smseq,obj.budget.(vars{i}).info.rtits);
                obj = computeNet(obj,vars(i),'fix',A.fix);

                % OPTIONAL SWITCHES
                % Integrate vertically (and horizontally)
                if A.int == 1
                    obj = intBudg(obj,vars(i),'fix',A.fix); 
                end
                % Split sources and sinks
                if A.srcsnk == 1
                    obj = sources_sinks(obj,vars(i));
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
        function obj = getFluxes(obj,vars,fluxes,lvls,varargin);
            % -------------------
            % Grab 2D fluxes for budget, convert to 3D 
            % Called in getBudg
            % Fluxes in mmol/m2/s
            %
            % Usage:
            % - obj = getFluxes(obj,vars,fluxes,lvls)
            %
            % Inputs:
            % - vars = BGC budget that you are closing (i.e. 'NO2')
            % - fluxes  = 2D fluxes (air-sea, sediment, etc)
            % - lvls    = levels to apply 2D flux ('sfc','sed', as a cell array)
            %
            % Optional inputs:
            % - type    = file type (avg, his, rst)
            %
            % Example:
            % - obj = getFluxes(obj,'N2O',{'FG_N2O'},{'sfc'})
            % -------------------
            disp('---------------------------------');
            disp('Get 2D fluxes, apply to 3D grid')

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);
            
            % Kill if no fluxes
            if isempty(fluxes{1});
                obj.budget.(vars{1}).fg = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
                obj.budget.(vars{1}).sed = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
                return
            end

            % Convert 2D flux to 3D based on lvls
            for i = 1:length(fluxes)
            
                % Initialize matrices-to-fill
                obj.budget.(vars{1}).fg = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
                obj.budget.(vars{1}).sed = zeros(size(obj.grid.(A.type).mask_rho3d)).*obj.grid.(A.type).mask_rho3d;
    
                % Apply 2D mask to 2D data
                obj.data.(A.type).(fluxes{i}).data = obj.data.(A.type).(fluxes{i}).data .* obj.grid.mask_rho;
                
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
            disp('---------------------------------');
            disp('Computing dC/dt')

            % Grab time info from average file
            fieldnames = {obj.info.bgc_avg(1).Attributes.Name};
            fieldvalue = {obj.info.bgc_avg(1).Attributes.Value};
            tmpdt      = fieldvalue(strcmp(fieldnames,'dt'));
            tmpdt      = double(tmpdt{1});
            nwrt       = fieldvalue(strcmp(fieldnames,'nwrt'));
            if nwrt{1} == 999999 % fill value, occasionally these files are created separately
                nwrt = fieldvalue(strcmp(fieldnames,'navg'));
                if nwrt{1} == 999999
                    disp('Check output frequency, then update computeDcDt');
                    kill
                end
            end
            nwrt = double(nwrt{1});

            % Compute dt according to output frequency
            dt = tmpdt*nwrt;

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

                % Mask padding
                obj.budget.(vars{v}).dcdt(end,:,:,:) = nan;
                obj.budget.(vars{v}).dcdt(:,end,:,:) = nan;
            end 
        end % end method computeDcDt

        %--------------------------------------------------------------------------------
        function obj = computeXYZflux(obj,vars,file)
            % -------------------
            % Compute HorXAdvFlux, HorYAdvFlux, and top/bottom advection
            % Also get diffusion, if it is available
            % Called in getBudg
            % End result is mmol/m3/s
            %
            % Usage:
            % - obj = computeXYZflux(obj,vars)
            %
            % Inputs:
            % - vars = variable to get flux terms for, as a cell array
            %
            % Example:
            % - obj = computeXYZflux(obj,{'NO2'});
            % -------------------
            disp('---------------------------------');
            disp('Get advective/diffusion terms')
        
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
            disp('---------------------------------');
            disp('Computing sources-minus-sinks (SMS)');

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

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
        end % end method getSMS

        %--------------------------------------------------------------------------------
        function obj = computeNet(obj,vars,varargin)
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
            disp('---------------------------------');
            disp('Computing budget remainder (net)');

            % Process inputs
            A.fix      = [0];
            A          = parse_pv_pairs(A,varargin);

            % Calculate remainder (net)
            obj.budget.(vars{1}).net = obj.budget.(vars{1}).dcdt - (obj.budget.(vars{1}).adx + ...
                                       obj.budget.(vars{1}).ady  +  obj.budget.(vars{1}).adz + ...
                                       obj.budget.(vars{1}).dfz  +  obj.budget.(vars{1}).sms + ...
                                       obj.budget.(vars{1}).fg   +  obj.budget.(vars{1}).sed);

            % Apply fixes?
            if A.fix == 1 & strcmp(vars{1},'N2')
                % Set 'net' at surface to FG_N2
                disp('Set N2 net at surface to FG_N2');
                obj.budget.N2.fg(:,:,end,:)  = obj.budget.N2.net(:,:,end,:);
                obj.budget.N2.net(:,:,end,:) = zeros(size(obj.budget.N2.net(:,:,end,:))) .* obj.grid.mask_rho; 
            end
            if A.fix == 1
                % Fix x and y advection (near coast, advection mismatches budget) 
				disp('Adjusting x/y advection to balance budget');
				tmpnet = abs(obj.budget.(vars{1}).net);
				tmpadx = abs(obj.budget.(vars{1}).adx);
				tmpady = abs(obj.budget.(vars{1}).ady);
				tmpdiffx = tmpnet - tmpadx;
				tmpdiffy = tmpnet - tmpady;
				advfix = tmpdiffx - tmpdiffy;

				% If > 0, fix y-advection
				obj.budget.(vars{1}).ady(advfix>0) = obj.budget.(vars{1}).ady(advfix>0) + obj.budget.(vars{1}).net(advfix>0);
				% If < 0, fix x-advection
				obj.budget.(vars{1}).adx(advfix<0) = obj.budget.(vars{1}).adx(advfix<0) + obj.budget.(vars{1}).net(advfix<0);

				% Recalculate horizontal and total advection
				obj.budget.(vars{1}).adxy = obj.budget.(vars{1}).adx + obj.budget.(vars{1}).ady;
				obj.budget.(vars{1}).adv  = obj.budget.(vars{1}).adx + obj.budget.(vars{1}).ady + obj.budget.(vars{1}).adz;

                % Re-Calculate remainder (net)
				% Should be 0!!!!
                obj.budget.(vars{1}).net = obj.budget.(vars{1}).dcdt - (obj.budget.(vars{1}).adx + ...
                                           obj.budget.(vars{1}).ady  +  obj.budget.(vars{1}).adz + ...
                                           obj.budget.(vars{1}).dfz  +  obj.budget.(vars{1}).sms + ...
                                           obj.budget.(vars{1}).fg   +  obj.budget.(vars{1}).sed);
            end
        end % end method computeNet
    
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
            disp('---------------------------------');
            disp('Integrating 3D variables');

            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);

            % Load depths?
            try; obj.grid.(A.type).Hz;
            catch; disp('call loadDepth first'); return
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
        function obj = intFlux(obj,vars,varargin)
            % ------------------
            % Totally integrate 2D variable(s) 
            %
            % Usage:
            % - obj = intFlux(obj,vars,varargin)
            %
            % Inputs:
            % - vars = 2D variables to integrate, as a cell array
            %
            % Optional:
            % - type = file type ('avg','his',etc)
            %
            % Example:
            % - obj = loadData(obj,{'FG_N2O'},file);
            % - obj = intVar(obj,{'FG_N2O'});
            % ------------------
            disp('---------------------------------');
            disp('Integrating 2D variables');
        
            % Process optional inputs
            A.type = 'avg';
            A = parse_pv_pairs(A,varargin);
    
            % Go through each 3D rate, integrate vertically (and totally)
            for i = 1:length(vars)        
                obj.data.(A.type).(vars{i}).int = obj.data.(A.type).(vars{i}).data;
                for t = 1:size(obj.data.(A.type).(vars{i}).data,3);
                    obj.data.(A.type).(vars{i}).tot(t) = nansum(obj.data.(A.type).(vars{i}).data(:,:,t) .*obj.grid.area_rho,'all');
                end
            end
        end % end method intVar

        %--------------------------------------------------------------------------------
        function obj = intBudg(obj,vars,varargin)
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
			% Optional Inputs:
			% - fix = (1) to apply corrections to x/y advective terms near coasts 
            %
            % Example:
            % - obj = intBudg(obj,{'NO2'},'fix',1); <-- to apply fixes
            % ------------------

			% Optional inputs
            A.fix      = [0];
            A          = parse_pv_pairs(A,varargin);

            disp('---------------------------------');
            disp('Computing vertical integration...');
            disp('...and getting grid area fluxes');
    
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

			% Apply correction to total advection
			if A.fix == 1
				% Correcting x/y advection near coasts
				disp('Correcting x/y errors near coast');
				% Calculate budget remainder
				mynet = obj.budget.(vars{1}).intdcdt - ...
					 (obj.budget.(vars{1}).intadv + obj.budget.(vars{1}).intfg + ...
					  obj.budget.(vars{1}).intsms + obj.budget.(vars{1}).intsed + obj.budget.(vars{1}).intdfz);
				% Add fixes to horizontal advection
				obj.budget.(vars{1}).intadxy = obj.budget.(vars{1}).intadxy + mynet;
				% Update total advection
				obj.budget.(vars{1}).intadv  = obj.budget.(vars{1}).intadxy + obj.budget.(vars{1}).intadz;
				% Save budget remainder
				obj.budget.(vars{1}).intmynet = mynet;
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

            % process inputs
            A.type     = 'avg';
            A          = parse_pv_pairs(A,varargin);    

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
                disp(['...grabbing ',vars{v},' profile(s)']);
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
        function [obj] = sources_sinks(obj,vars)
            % ------------------------------------------------
            % Function to extract sources and sinks
            % Rates are in the correct units (see obj.budget.(vars).info)
            %
            % Usage: 
            % - [obj] = sources_sinks(obj,vars) 
            % 
            % Inputs:
            % - obj:    romsObj object
            % - vars:   budget variable
            %
            % Example:
            % - [obj] = sources_sinks(obj,{'NH4'});
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
        end % end method sources_and_sinks
    
        %--------------------------------------------------------------------------------
        function plotIntBudg(obj,vars,varargin);
            % ------------------
            % Plot the vertcially integrated budget terms
            %
            % Usage:
            % - plotIntBudg(obj,vars,varargin)
            %
            % Inputs (varargin):
            % - time = time record to plot (if length > 1, plot the average)
            % - prc  = percentile to limit colorbar
            %       ...if prc == 5, then caxis = [5th %, 95th %]
            %
            % Example:
            % - plotIntBudg(obj,{'N2O'},'time',10,'prc',5)
            % ------------------

            % defaults for optional  arguments
            A.time      = [];
            A.prc       = [2];
            A.lonbounds = [floor(min(obj.grid.lon_rho(:))) ceil(max(obj.grid.lon_rho(:)))];
            A.latbounds = [floor(min(obj.grid.lat_rho(:))) ceil(max(obj.grid.lat_rho(:)))];
            A.font      = 8;
            A           = parse_pv_pairs(A,varargin); % parse method arguments to A
            
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
            terms = {'dcdt','adv','adxy','adx','ady','adz','dfz','sms','fg','sed','net','mynet'};
            tits  = {'d$C$/dt','$T$','$T_{h}$','$T_{x}$','$T_{y}$','T$_{z}$','$D$','$J$','$\Phi$','Sed','Net','MyNet'};
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
                    tmplims(1:2) = prclims(tmp.(terms{j}),'prc',A.prc,'bal',1);
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
                        'ticks',1,...
                        'fontsize',A.font);
                    caxis([lims(1) lims(end)]);
                    hold on
                    m_plot(obj.grid.lon_rho(1,:),obj.grid.lat_rho(1,:),'--k','linewidth',2);
                    m_plot(obj.grid.lon_rho(:,1),obj.grid.lat_rho(:,1),'--k','linewidth',2);
                    m_plot(obj.grid.lon_rho(end,:),obj.grid.lat_rho(end,:),'--k','linewidth',2);
                    m_plot(obj.grid.lon_rho(:,end),obj.grid.lat_rho(:,end),'--k','linewidth',2);
                    title([obj.budget.(vars{i}).info.tits{1},': ',tits{j}],'Interpreter', 'Latex','FontSize',A.font+2);
                    ylabel(cb,obj.budget.(vars{1}).info.units{2},'Interpreter','Latex');
                    fname = [vars{i},'_int',terms{j}];
                    if exist([obj.paths.plots.budget,fname,'.png']) == 2
                        cmd = ['rm ',obj.paths.plots.budget,fname,'.png'];
                        system(cmd);
                    end
                    export_fig('-png',[obj.paths.plots.budget,fname],'-m5');
                    close(fig)    
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
            A           = parse_pv_pairs(A,varargin); % parse method arguments to A

            % Check for integrated budget
            try; obj.budget.(vars{1}).totdcdt;
            catch; obj = intBudg(obj,vars);
            end

            % Get terms to plot
            terms = {'dcdt','adv','dfz','sms','fg','sed','net'};
            tits  = {'d$C$/dt','$T$','$D$','$J$','$\Phi$','Sed','Net'};
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

				% Grab 'mynet'
				mynet = obj.budget.(vars{i}).totdcdt - (...
					obj.budget.(vars{i}).totadv + ...
					obj.budget.(vars{i}).totsms + ...
					obj.budget.(vars{i}).totfg  + ...
					obj.budget.(vars{i}).totsed + ...
					obj.budget.(vars{i}).totdfz);
				if length(A.time)==1
					tmp.mynet.mean = mynet(A.time);
					tmp.mynet.std  = 0;
				else
					tmp.mynet.mean = nanmean(mynet);
					tmp.mynet.std  = nanstd(mynet);
				end
                terms{end+1} = 'mynet';
                tits{end+1}  = 'MyNet';

                % Plot total terms
                fig = piofigs('mfig',1);
                for j = 1:length(terms)
                    y(j)   = [tmp.(terms{j}).mean];
                    ye(j)  = [tmp.(terms{j}).std];
					if ismember(vars{i},{'NO3','NO2','NH4'})
						y(j) = y(j) .* obj.constants.mmolNs_TgNy;
						ye(j) = ye(j) .* obj.constants.mmolNs_TgNy;
					elseif ismember(vars{i},{'N2O','N2'}); 
						y(j) = y(j) .* obj.constants.mmolN2s_TgNy;
						ye(j) = ye(j) .* obj.constants.mmolN2s_TgNy;
					end
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
				if ismember(vars{i},{'NO3','NO2','NH4','N2','N2O'})
					ylabel('Tg N yr$^{-1}$','Interpreter','Latex');
				else
					ylabel(obj.budget.(vars{i}).info.units{3},'Interpreter','Latex');
				end
                fname = [vars{i},'_budget'];
                if exist([obj.paths.plots.budget,fname,'.png']) == 2
                    cmd = ['rm ',obj.paths.plots.budget,fname,'.png'];
                    system(cmd);
                end
                % Save figure
                export_fig('-png',[obj.paths.plots.budget,fname],'-m5');
            end
        end % end method plotTotBudg

        %--------------------------------------------------------------------------------
        function dispVars(obj,varargin)
            % ----------------------
            % Lists output variables in command window
            % Assumes all files are similar to the first file
            %
            % Usage:
            % dispVars(obj,varargin)    
            %
            % Optional inputs (varargin):
            % - type = 'avg','his','rst','bgc','phy','phys_flux' 
            %
            % Example:
            % - dispVars(obj,'type','his');
            % ----------------------
        
            % Default arguments    
            A.type = 'avg';
            A      = parse_pv_pairs(A,varargin);
    
            % Call ncdisp    
            ncdisp([obj.paths.(A.type),obj.info.(A.type).files{1}]);    
        end % end method dispVars

        %--------------------------------------------------------------------------------
        function [fig,ax] = gridView(obj,varargin)
            % ----------------------
            % Plots the lat/lon indices of a grid file
            % Useful for identifying transect, or lat/lon indices
            %
            % Usage:
            % - [fig,ax] = gridView(obj,varargin)
            %
            % Inputs (varargin):
            % - full  = 1 (default), view entire grid. Use 0 for regional grid
            % - dx    = plot lon lines separated by dx (default = 20)
            % - dy    = plot lat lines separated by dy (default = 20)
            % - ticks = 0 (no lon/lat labels), 1 (yes), 2 (fancy box)
            % - font  = tick font size (default = 12) 
            % - save  = 1 (print), 0 (no print)
            %
            % Example:
            % - [fig,ax] = gridView(obj,'dx',20,'dy',20)
            % ----------------------

            % Grab inputs (varargin)
            A.full  = 1;
            A.dx    = [20];
            A.dy    = [20];
            A.ticks = [1];
            A.font  = [10];
            A.save  = [0];
            A       = parse_pv_pairs(A,varargin);

            % Get grid
            if A.full == 1
                % Get original lon/lat
                tmp.lon_rho = double(ncread(obj.paths.grid,'lon_rho'));
                tmp.lat_rho = double(ncread(obj.paths.grid,'lat_rho'));
            else
                tmp.lon_rho = obj.grid.lon_rho;
                tmp.lat_rho = obj.grid.lat_rho;
            end

            % Plot lon/lat lines
            fig(1) = piofigs('mfig',1.5);
            ax(1)  = map_plot(fig(1),tmp.lon_rho,tmp.lat_rho,'ticks',A.ticks,'font',A.font);    
            [a,b]  = size(tmp.lon_rho);
            for i = 1:A.dx:a
                m_plot(tmp.lon_rho(i,:),tmp.lat_rho(i,:),'r');
                m_text(tmp.lon_rho(i,end),tmp.lat_rho(i,end),num2str(i),'fontsize',8);
                for j = 1:A.dy:b
                    hold on
                    m_plot(tmp.lon_rho(:,j),tmp.lat_rho(:,j),'b');
                    m_text(tmp.lon_rho(end,j),tmp.lat_rho(end,j),num2str(j),'fontsize',8);
                end
            end
            hold on
            m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'k','linewidth',1);
            fname = [obj.info.simName,'_grid'];
            if A.save == 1
                export_fig('-jpg',[obj.paths.runPath,fname]);
            end
            if nargout < 1
                pltjpg(1);
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
            A.ticks = [1];
            A.font  = [10];
            A.save  = [0];
            A       = parse_pv_pairs(A,varargin);

            % Get original lon/lat
            tmp.lon_rho = ncread(obj.paths.grid,'lon_rho');
            tmp.lat_rho = ncread(obj.paths.grid,'lat_rho');
    
            % Map of region
            if ~isempty(obj.grid.lon_rho)
                % Generate whole map
                fig(1)  = piofigs('mfig',1.5);
                set(0,'CurrentFigure',fig(1));
                [ax(1)] = map_plot(fig(1),obj.grid.lon_rho,obj.grid.lat_rho,'ticks',A.ticks,'font',A.font);
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
                    export_fig('-jpg',[obj.paths.runPath,fname]);
                else
                    pp = 1;
                    pltjpg(1);
                end
                if nargout < 1 & pp ~= 1;
                    pltjpg(1);
                end
            end
        end % end method regionView

        %--------------------------------------------------------------------------------
        function obj = Dist2Coast(obj)
            % --------------------
            % Call this function to calculate each grid cell's distance-to-coast
            % --------------------
            disp('Calculating distance from coast');

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
            nch(1) = 1;
            nch(2) = 3;
            nch(3) = 4;
            disp(' Getting distance from coast');
            disp(' Starting the looper')
            for it = 1:1
                ncx = nch(it);
                ncy = nch(it);
                for ic = 1:ncx
                    for jc = 1:ncy
                        [ic jc];
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
                        else
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
            % - meta:        option to include path to name and units (e.g., obj.data.avg.temp)
            % - lonbounds:  x-boundaries (defaults to whole domain)
            % - latbounds:  y-boundaries (defaults to whole domain)
            % - ticks:      2 = fancy, 1 = on, 0 = off
            % - background: background color (default 'LightGray');
            % - coastcolor: coast color (default 'DimGray');
            % - fontsize:   default 10
            % - figtype:    default 'mfig' (140mm wide), can use 'sfig' (90mm) or 'lfig' (190mm)
            % - figdim:     default 1 (same height as width, its a multiplier)
            % - prc:        percentage to limit colorbar axes
            % - bal:        balance colorbar around 0 (1 == yes, 0 == no)
            % - levels:     hard-coded levels to plot (can't be used with A.prc)
            % - cmap:       colormap(default = thermal for no balance, balance for balance)
            % - caxis       colorbar limits 
            % - log         log-scale (1), use with caxis to set limits
            % -----------------------
            
            % User-inputs
            A.meta       = [];
            A.lonbounds  = [];
            A.latbounds  = [];
            A.lonticks   = [];
            A.latticks   = [];
            A.ticks      = 0;
            A.background = rgb('LightGray');
            A.coastcolor = rgb('DimGray');
            A.fontsize   = 10;
            A.figtype    = 'mfig';
            A.figdim     = 1;
            A.prc        = 2;
            A.bal        = 0;
            A.levels     = [];
            A.cmap       = [];
            A.caxis      = [];
            A.log        = 0;
            A = parse_pv_pairs(A,varargin);

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

            % Initiate figure
            fig = piofigs(A.figtype,A.figdim);

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
			if isempty(A.cmap)
				if min(clims(:))<0 & max(clims(:)) > 0
					A.cmap = cmocean('balance',30);
				else
					A.cmap = cmocean('thermal',30);
				end
			end

            % Make map
            set(0,'CurrentFigure',fig);
            m_proj('mercator','lat',latbounds,'lon',lonbounds); drawnow
            hold on
            if A.ticks == 0
                m_grid('box','on','linestyle','none','xtick',0,'ytick',0,...
                       'xticklabels',[],'yticklabels',[],'backgroundcolor',rgb('LightGray')); drawnow
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
            m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
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
            
            % Print figure
            if nargout<1
                pltjpg(1);
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
            % - cmap:       colormap(default = thermal for no balance, balance for balance)
            % - units:      string containing units for colorbar
            % ----------------------

            % User-inputs
            A.lonbounds  = [];
            A.latbounds  = [];
            A.ticks      = 0;
            A.background = rgb('LightGray');
            A.coastcolor = rgb('DimGray');
            A.fontsize   = 10;
            A.figtype    = 'mfig';
            A.figdim     = 1;
            A.prc        = 0.1;
            A.bal        = 0;
            A.levels     = [];
            A.difflevels = [];
            A.cmap       = 'thermal';
            A.units      = [];
            A = parse_pv_pairs(A,varargin);

            % Balance override
            if A.bal == 1
                A.cmap = cmocean('balance');
            end        
            
            % Get universal levels
            if isempty(A.levels)
                all_dat = [dat1(:) dat2(:)];
                A.levels = romsMaster.prclims(all_dat,'prc',A.prc,'bal',A.bal);
                A.levels = linspace(A.levels(1),A.levels(2),20);
            end

            % Make figs(1) and figs(2)
            [figs(1),cbs(1)] = mapPlot(obj,dat1,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',cmocean(A.cmap,length(A.levels)-1));
            [figs(2),cbs(2)] = mapPlot(obj,dat2,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.levels,'cmap',cmocean(A.cmap,length(A.levels)-1));

            % Get differences
            diff_dat  = dat1 - dat2;
            if isempty(A.difflevels)
                A.difflevels = romsMaster.prclims(diff_dat,'prc',A.prc,'bal',1); 
                A.difflevels = linspace(A.difflevels(1),A.difflevels(2),20);
            end
            
            % Make figs(3)
            [figs(3),cbs(3)] = mapPlot(obj,diff_dat,...
                'lonbounds',A.lonbounds,'latbounds',A.lonbounds,'ticks',A.ticks,...
                'background',A.background,'coastcolor',A.coastcolor,'fontsize',A.fontsize,...
                'figtype',A.figtype,'figdim',A.figdim,'levels',A.difflevels,'cmap',cmocean('balance',length(A.difflevels)-1));

            % Add units?
            if ~isempty(A.units)
                ylabel(cbs(1),A.units,'Interpreter','Latex','FontSize',A.fontsize);
                ylabel(cbs(2),A.units,'Interpreter','Latex','FontSize',A.fontsize);
                ylabel(cbs(3),A.units,'Interpreter','Latex','FontSize',A.fontsize);
            end

            % Auto-print if no output
            if nargout == 0
                set(0,'CurrentFigure',figs(1));
                pltjpg(1);
                
                set(0,'CurrentFigure',figs(2));
                pltjpg(2);

                set(0,'CurrentFigure',figs(3));
                pltjpg(3);
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
            A.lonbounds  = [];
            A.latbounds  = [];
            A.lonticks   = [];
            A.latticks   = [];
            A.ticks      = 0;
            A.ytick      = 1;
            A.xtick      = 1;
            A.background = rgb('LightGray');
            A.coastcolor = rgb('DimGray');
            A.fontsize   = 10;
            A.prc        = 2;
            A.bal        = 0;
            A.levels     = [];
            A.cmap       = cmocean('thermal');
            A.caxis      = [];
            A.log        = 0;
            A = parse_pv_pairs(A,varargin);

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
                       'xticklabels',[],'yticklabels',[],'backgroundcolor',rgb('LightGray')); drawnow
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
            m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
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
                disp('No point using this without [fig] output');
                return
            end

            % User-inputs
            A.lonbounds  = [obj.grid.minlon_rho obj.grid.maxlon_rho];
            A.latbounds  = [obj.grid.minlat_rho obj.grid.maxlat_rho];
            A.latticks   = [];
            A.lonticks   = [];
            A.ticks      = 0;
            A.box        = 'on';
            A.background = rgb('LightGray');
            A.coastcolor = rgb('DimGray');
            A.fontsize   = 7;
            A.figtype    = 'mfig';
            A.figdim     = 1;
            A.poly       = 1;
            A = parse_pv_pairs(A,varargin);

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
            if isempty(A.lonticks) & ismember(A.ticks,[1 2]);
                if abs(diff(lonbounds)) > 50
                    dx = 10;
                else
                    dx = 5;
                end
                lonticks  = (lonbounds(1):floor(range(lonbounds)/dx):lonbounds(2));
            elseif ismember(A.ticks,[1 2]);
                lonticks = A.lonticks;
            end
            if isempty(A.latticks) & ismember(A.ticks,[1 2])
                if abs(diff(latbounds)) > 50
                    dy = 10;
                else
                    dy = 5;
                end
                latticks  = (latbounds(1):floor(range(latbounds)/dy):latbounds(2));
            elseif ismember(A.ticks,[1 2]);
                latticks = A.latticks;
            end

            % Make map
            fig = piofigs(A.figtype,A.figdim);
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
            m_coast('patch',rgb('DimGray'),'edgecolor','k'); drawnow
    
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
            A.xlims      = [];
            A.zlims      = [];
            A.cmap       = [];
            A.fontsize   = 10;
            A.figtype    = 'mfig';
            A.figdim     = 0.33;
            A.prc        = 0.5;
            A.bal        = 0;
            A.levels     = [];
            A.background = rgb('DimGray');
            A = parse_pv_pairs(A,varargin);

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
            fig = piofigs(A.figtype,A.figdim);
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
                pltjpg(1);
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

            % Get optional inputs
            A.zlim = inf;
            A = parse_pv_pairs(A,varargin);

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
                disp(' '); disp(['Slicing validation ',vars{ff},' data...']);disp(' ');
                % get path and coords for current variable
                if ~isfield(obj.paths.diag,(vars{ff}));
                    disp([vars{ff},' is not a diagnostic variable']);
                    disp(['Filling with NaN']);
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
                        % record progress
                        fprintf([num2str(rcrd),'...']);
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
                    fprintf(['\n']);
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
            A.xlims      = [];
            A.zlims      = [];
            A.cmap       = [];
            A.fontsize   = 10;
            A.figtype    = 'mfig';
            A.figdim     = 0.33;
            A.prc        = 2;
            A.bal        = 0;
            A.levels     = [];
            A.difflevels = [];
            A = parse_pv_pairs(A,varargin);

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

            % Optional arguments
            A.outer = [3]; % interpolant padding (degrees)
            A=parse_pv_pairs(A,varargin);

            % Check inputs
            diagfields = fieldnames(obj.paths.diag); 
            for i = 1:length(vars)
                if ~strcmp(vars{i},diagfields) & ~strcmp(vars{i},upper(diagfields));
                    disp(' ');
                    disp([vars{i},' is not a diagnostic variable']);
                    disp(['Filling with NaN']);
                    obj.diag.(vars{i}).slice = nan(obj.grid.nx,obj.grid.ny,length(depth),12);
                    obj.diag.(vars{i}).name  = ' ';
                    obj.diag.(vars{i}).units = ' ';
                    obj.diag.(vars{i}).depth = depth;
                    disp(' ');
                    skip(i) = 1;
                else
                    skip(i) = 0;
                end
            end
            if min(depth) < 0 | max(depth) > max(obj.grid.z_avg_dep)
                disp(' ');
                disp('Check depth input');
                disp(' '); return
            end
            for i = 1:length(depth)
                diffd = abs(depth(i) - obj.grid.z_avg_dep);
                ind   = find(diffd == min(diffd));
                depth(i) = obj.grid.z_avg_dep(ind);
                disp(['Closest depth = ',num2str(depth(i))])
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
                fprintf(['\n Processing ', vars{i}]);

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
                    fprintf('\n month:');
                    skip = 0;
                    for k = 1:12
                        tmpdata = []; tmpdepth = [];
                        fprintf([num2str(k),'...']);

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
                                skip = 1;
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
                                skip = 1;
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
                                disp(['No z-data for ',vars{i},', skipping']);
                                continue
                            end
                            % Interpolate
                            F = scatteredInterpolant(double(tmp.lon(idx)),double(tmp.lat(idx)),...
                                                     double(data(idx)),'linear','nearest');
                            tmpout{z,k} = F(double(tmp.outlon),double(tmp.outlat));
                            tmpout{z,k}(isnan(obj.grid.mask_rho)) = nan;
                        end % end z-loop
                    end % end k-loop
                    fprintf('\n');
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
            disp('Clearing loaded data');            

            % Clear fields
            obj.slice = [];
            obj.profile = [];
            obj.diag = [];
            obj.data = [];
            obj.budget = [];

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
        function o2_corr = woao2corr(o2)
            % ------------------
            % -  applies Bianchi 2012 correction to O2
            % ------------------
                
            o2_corr=1.009*o2-2.523;
            o2_corr(o2_corr<0)=0;
        end % end static method woao2corr

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
        function loglevs=dfloglevs(levs,logminparam)
            % ------------------
            % - converts positive levels into symmetrical centered around 0 levels
            % - in log space. For use in difference plots with log scales  
            % ------------------
                
            levs(abs(levs)<logminparam)=sign(levs(abs(levs)<logminparam)).*logminparam;
            loglevs=levs;
            loglevs(levs>0)=+abs(log10(logminparam) - (log10(levs(levs>0))));
            loglevs(levs<0)=-abs(log10(logminparam) - log10(abs(levs(levs<0))));
        end % end static method dfloglevs

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
        function [mnmx,levs] = oom_levs(mnmx,varargin);
            % ----------------
            % - Gets order of magnitude estimate and sets difference levels for axis labels
            % - Essentially, used to automatically set axes such that there are 10 levels
            % - and the levels are nice rounded values rather than 16.666666 etc
            % ----------------
        
            A.diff = 0; % default to not centered around 0
            A = parse_pv_pairs(A,varargin);
    
            % Adjust min and max levels
            di     = diff(mnmx);
            di_oom = floor(log10(di));
            if di_oom >= 0
                mnmx = round(mnmx,-di_oom+1);
            else
                mnmx = round(mnmx,-di_oom);
            end 
            if A.diff == 1        
                mnmx = [-max(abs(mnmx)) max(abs(mnmx))];
            end

            % Get levels
            dabs = [diff(mnmx)/9];
            oom  = floor(log10(dabs));
            doom = [dabs/10^(oom)];
            doom = round(doom);
            dabs = doom*10^(oom);
            levs = [mnmx(1):dabs:mnmx(2)];
            if mnmx(2) > levs(end)
                levs = [levs levs(end)+dabs];
            end
            if A.diff == 0
                mnmx = [levs(1) levs(end)];
            elseif A.diff == 1
                tmplevs = [0:dabs:mnmx(2)];
                if mnmx(2) > tmplevs(end)
                    tmplevs = [tmplevs tmplevs(end)+dabs];
                end
                levs = unique([-tmplevs tmplevs]);
                mnmx = [levs(1) levs(end)];
            end    
        end % end static method oom_levs

        %--------------------------------------------------------------------------------
        function [ws,wsc] = WindStress(u,v,lon,lat,ang)
            % ------------------
            % - calculates wind stress and wind stress curl from zonal/meridional wind stress 
            % ------------------
            
            % - calculate wind stress
            ws = sqrt(u.*u + v.*v);

            % - check that lon/lat are gridded
            [a,b] = size(lon);
            if a == 1 | b == 1
                [lon,lat] = meshgrid(lon,lat);
            end
        
            % - check for yearly or monthly file
            [a,b,c] = size(u);

            % - track progress
            d        = a*b*c;
            progress = [d/10:d/10:d];
            cnt      = 0;
            pcnt     = 10;
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
                        cnt = cnt + 1;
                        if ismember(cnt,progress)
                            disp([num2str(pcnt),'% complete']);
                            pcnt = pcnt + 10;
                        end
                    end
                end
            end
            end % end static method WindStress

        %--------------------------------------------------------------------------------
        function [lon,lat] = glodap_coord
            % -------------------
            % - Function to grab allowable coordinates for GLODAPv2 comparisons
            % -------------------

            % Load GLODAPv2.2020 data
            load('/data/project1/demccoy/ROMS/validation/GLODAPv2/glodap_v2_2020_format.mat');
            tmplon  = [glodap.lon]; tmplon(tmplon<0) = tmplon(tmplon<0)+360;
            tmplat  = [glodap.lat];
            tmptime = [glodap.time];

            % Plot map of locations
            figure
            set(gcf,'Position',[451          75        1376         946],'color','w');
            plot(tmplon,tmplat,'.k'); hold on; plot_coast; xlim([0 360]); ylim([-90 90]);
            set(gca,'XTick',[0:10:360],'YTick',[-90:10:90]);
            grid on

            % Get user input of lon/lat
            disp('Choose location for 5x5 grid box');
            for i = 1:10
                if i > 1
                    disp('Try again');
                    delete(r);
                end

                % Query lon/lat box
                x = []; y = []; idx = [];
                            [x] = input('---------------------------\nLon/lat estimate ([lon,lat]):\n---------------------------\n>> ');
                [y] = x(2); 
                [x] = x(1);

                % Find number of profiles within 5x5 box
                idx   = find(x-2.5 <= tmplon & tmplon <= x+2.5 & y-2.5 <= tmplat & tmplat <= y+2.5);
                r     = rectangle('Position',[x-2.5 y-2.5 5 5],'FaceColor','b');
                
                % Query keep or reject
                disp([num2str(length(idx)),' profiles found']);
                            q = input('---------------------------\nOK? (1 == YES):\n---------------------------\n>> ');
                if q == 1 & ~isempty(idx)
                    break
                else
                    continue
                end
            end

            % Save results
            lon = x;
            lat = y;

        end % end static method glodap_coord

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
        function [salt_abs] = getSA(salt,pres,lon,lat)
            % -------------------
            % Function to call gsw_SA_from_SP, but in parallel
            % 
            % Usage:
            % - [salt_abs] = getSA(salt,pres,lon,lat);
            %
            % Inputs:
            % - salt, pres, lon, lat = ROMS variables...must be the same size!
            % 
            % -------------------

            % Create filler
            salt_abs = NaN(size(salt));
            dims     = size(salt);
            if length(dims)>3
                nt = dims(4);
            else
                nt = 1;
            end
            nz = dims(3);
            if nt == 1
                i = 1;
                for z = 1:nz
                    tmp_salt = squeeze(salt(:,:,z,i));
                    tmp_pres = squeeze(pres(:,:,z,i));
                    tmp_lon  = squeeze(lon(:,:,z,i));
                    tmp_lat  = squeeze(lat(:,:,z,i));
                    tmp_SA   = gsw_SA_from_SP(tmp_salt(:),tmp_pres(:),tmp_lon(:),tmp_lat(:));
                    salt_abs(:,:,z,i) = reshape(tmp_SA,dims(1),dims(2));
                end
            else
                tic
                % Start parpool
                delete(gcp('nocreate'));
                parpool(6);
                parfor i = 1:nt
                    for z = 1:nz;
                        tmp_salt = squeeze(salt(:,:,z,i));
                        tmp_pres = squeeze(pres(:,:,z,i));
                        tmp_lon  = squeeze(lon(:,:,z,i));
                        tmp_lat  = squeeze(lat(:,:,z,i));
                        tmp_SA   = gsw_SA_from_SP(tmp_salt(:),tmp_pres(:),tmp_lon(:),tmp_lat(:));
                        salt_abs(:,:,z,i) = reshape(tmp_SA,dims(1),dims(2));
                    end
                end
                toc
                delete(gcp('nocreate'));
            end
        end % end static method getSA

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
            A     = parse_pv_pairs(A,varargin);

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
                    fprintf([num2str(t),'...']);
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
            % Transfert a field at psi points to the rho points
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
            % Transfer a field at v points to the rho points
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
        function [out,grd] = scatter_density(xdata,ydata,xbounds,ybounds)
            % --------------------------------------------------------------------
            % Probability heatmap given xbounds and ybounds, based on
            % independent X variable and dependent Y variable
            % Usage: 
            % - [out,grd] = scatter_density(xdata,ydata,xbounds,ybounds)  
            % 
            % Inputs:
            % - xdata, ydata     = independent, dependent variables (matrices ok)
            % - xbounds, ybounds = increasing xbin/ybin limits
            %
            % Outputs:
            % - out = count, probability, mean, and median results for each xbin
            % - grd = meshgrid for plotting
            % --------------------------------------------------------------------

            % Initialize matrices
            try
                dat = xdata(:)+ydata(:);
            catch
                disp('debug scatter_density');
                keyboard
            end
            xdata = xdata(~isnan(dat)); 
            ydata = ydata(~isnan(dat));
            xdata = xdata(:);
            ydata = ydata(:);
            xbounds = xbounds(:); ybounds = ybounds(:);
            out.count  = nan(length(ybounds)-1,length(xbounds)-1);
            out.prob   = out.count;
            out.mean   = nan(1,length(xbounds)-1);
            out.median = out.mean;

            % Get output grid
            grd.xgrid = 0.5*(xbounds(1:end-1)+xbounds(2:end));
            grd.ygrid = 0.5*(ybounds(1:end-1)+ybounds(2:end));
            [grd.X grd.Y]  = meshgrid(grd.xgrid,grd.ygrid);

            % Get 2D histogram
            for i = 1:length(xbounds)-1
                xi = xbounds(i);
                xf = xbounds(i+1);
                indx = find(xi<=xdata & xdata<xf);
                if isempty(indx)
                    out.count(:,i) = zeros(length(ybounds)-1,1);
                    out.prob(:,i)  = zeros(length(ybounds)-1,1);;
                    out.mean(1,i)  = NaN;
                    out.median(1,i) = NaN;
                    continue
                else
                    tmpdat = ydata(indx);
                    for j = 1:length(ybounds)-1
                        yi = ybounds(j);
                        yf = ybounds(j+1);
                        if j == 1
                            indy = find(tmpdat<yf);
                        elseif j == length(ybounds)-1
                            indy = find(tmpdat>=yi);
                        else
                            indy = find(yi<=tmpdat & tmpdat <yf); 
                        end
                        if isempty(indy)
                            out.count(j,i) = 0;
                            out.prob(j,i)  = 0;
                        else
                            out.count(j,i) = length(indy);
                            out.prob(j,i)  = 100*(length(indy)./length(tmpdat));
                        end
                    end
                    out.mean(1,i)   = nanmean(tmpdat);
                    out.median(1,i) = nanmedian(tmpdat); 
                end
            end
        end % end static method scatter_density

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
    end % end static methods declarations
end % end classdef
%------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
% Local functions
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
function [params] = getParams 
    % ----------------------
    % Get N-cycle parameters
    filedir = ['/data/project1/demccoy/iNitrOMZ/analysis/run/bgc_params.mat'];
    tmp = load(filedir);
    params = tmp.params;
end % end method getParams

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
    paths.diag.temp.file  = {'/data/project3/data/woa18/temperature/0p25/temp_woa18_clim.nc',...
                             '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.temperature.nc'};
    paths.diag.temp.type  = {'nc','nc'};
    paths.diag.temp.var   = {'temp','temperature'};
    paths.diag.temp.zvar  = {'depth','Depth'};
    paths.diag.temp.dim   = {'xyzt','xyz'};
    paths.diag.temp.lon   = {'lon','lon'};
    paths.diag.temp.lat   = {'lat','lat'};
    paths.diag.temp.name  = {'WOA-18','GLODAPv2'};
    paths.diag.temp.units = {'$^oC$','$^oC$'};

    % salinity
    paths.diag.salt.file  = {'/data/project3/data/woa18/salinity/0p25/salt_woa18_clim.nc',...
                             '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.salinity.nc'};
    paths.diag.salt.type  = {'nc','nc'};
    paths.diag.salt.var   = {'salt','salinity'};
    paths.diag.salt.zvar  = {'depth','Depth'};
    paths.diag.salt.dim   = {'xyzt','xyz'};
    paths.diag.salt.lon   = {'lon','lon'};
    paths.diag.salt.lat   = {'lat','lat'};
    paths.diag.salt.name  = {'WOA-18','GLODAPv2'};
    paths.diag.salt.units = {'PSU','PSU'};

    % density
    paths.diag.rho.file  = {'/data/project3/data/woa18/density/0p25/sigma_woa18_clim.nc'};
    paths.diag.rho.type  = {'nc'};
    paths.diag.rho.var   = {'sigma'};
    paths.diag.rho.zvar  = {'depth'};
    paths.diag.rho.dim   = {'xyzt'};
    paths.diag.rho.lon   = {'lon'};
    paths.diag.rho.lat   = {'lat'};
    paths.diag.rho.name  = {'WOA-18 Density'};
    paths.diag.rho.units = {'kg m$^{-3}$'};

    % sea surface height
    paths.diag.SSH.file  = {'/data/project1/demccoy/ROMS/validation/AVISO/monthly_AVISO.mat'};
    paths.diag.SSH.type  = {'mat'};
    paths.diag.SSH.var   = {'adt_month_av'};
    paths.diag.SSH.zvar  = {[]};
    paths.diag.SSH.dim   = {'xyt'};
    paths.diag.SSH.lon   = {'lon_aviso'};
    paths.diag.SSH.lat   = {'lat_aviso'};
    paths.diag.SSH.name  = {'AVISO+ Absolute Dynamic Topography'};
    paths.diag.SSH.units = {'m'};

    % zonal wind stress
    paths.diag.zws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
    paths.diag.zws.type  = {'mat'};
    paths.diag.zws.var   = {'u'};
    paths.diag.zws.zvar  = {[]};
    paths.diag.zws.dim   = {'xyt'};
    paths.diag.zws.lon   = {'lon'};
    paths.diag.zws.lat   = {'lat'};
    paths.diag.zws.name  = {'SCOW-2010 Zonal Wind Stress'};
    paths.diag.zws.units = {'N m$^{-2}$'};

    % meridional wind stress
    paths.diag.mws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
    paths.diag.mws.type  = {'mat'};
    paths.diag.mws.var   = {'v'};
    paths.diag.mws.zvar  = {[]};
    paths.diag.mws.dim   = {'xyt'};
    paths.diag.mws.lon   = {'lon'};
    paths.diag.mws.lat   = {'lat'};
    paths.diag.mws.name  = {'SCOW-2010 Meridional Wind Stress'};
    paths.diag.mws.units = {'N m$^{-2}$'};
    
    % wind stress
    paths.diag.ws.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
    paths.diag.ws.type  = {'mat'};
    paths.diag.ws.var   = {'ws'};
    paths.diag.ws.zvar  = {[]};
    paths.diag.ws.dim   = {'xyt'};
    paths.diag.ws.lon   = {'lon'};
    paths.diag.ws.lat   = {'lat'};
    paths.diag.ws.name  = {'SCOW-2010 Wind Stress'};
    paths.diag.ws.units = {'N m$^{-2}$'};

    % wind stress curl
    paths.diag.wsc.file  = {'/data/project1/demccoy/ROMS/validation/wind_stress/monthly_WSTRESS.mat'};
    paths.diag.wsc.type  = {'mat'};
    paths.diag.wsc.var   = {'wsc'};
    paths.diag.wsc.zvar  = {[]};
    paths.diag.wsc.dim   = {'xyt'};
    paths.diag.wsc.lon   = {'lon'};
    paths.diag.wsc.lat   = {'lat'};
    paths.diag.wsc.name  = {'SCOW-2010 Wind Stress Curl'};
    paths.diag.wsc.units = {'N m$^{-2}$'};

    % u velocity
    paths.diag.u.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
    paths.diag.u.type  = {'mat'};
    paths.diag.u.var   = {'U'};
    paths.diag.u.zvar  = {'dep'};
    paths.diag.u.dim   = {'xyzt'};
    paths.diag.u.lon   = {'lon'};
    paths.diag.u.lat   = {'lat'};
    paths.diag.u.name  = {'GODAS U-Velocity'};
    paths.diag.u.units = {'m s$^{-1}'}; 

    % v velocity
    paths.diag.v.file  = {'/data/project1/demccoy/ROMS/validation/uv/GODAS_uv.mat'};
    paths.diag.v.type  = {'mat'};
    paths.diag.v.var   = {'V'};
    paths.diag.v.zvar  = {'dep'};
    paths.diag.v.dim   = {'xyzt'};
    paths.diag.v.lon   = {'lon'};
    paths.diag.v.lat   = {'lat'};
    paths.diag.v.name  = {'GODAS V-Velocity'};
    paths.diag.v.units = {'m s$^{-1}'}; 

    % Eddy kinetic energy
    paths.diag.EKE.file  = {'/data/project1/data/AVISO/EKE_clim/aviso_1993_2022_EKE_monthly_climatology.mat'};
    paths.diag.EKE.type  = {'mat'};
    paths.diag.EKE.var   = {'eke'};
    paths.diag.EKE.zvar  = {[]}; 
    paths.diag.EKE.dim   = {'xyt'};
    paths.diag.EKE.lon   = {'lon'};
    paths.diag.EKE.lat   = {'lat'};
    paths.diag.EKE.name  = {'AVISO+ Ssalto/Duacs EKE climatology (1993-2023)'};
    paths.diag.EKE.units = {'cm$^2$ s$^{-2}$'};

    % oxygen
    paths.diag.O2.file  = {'/data/project3/data/woa18/oxygen/1p0/o2_woa18_clim.nc',...
                           '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.oxygen.nc'};
    paths.diag.O2.type  = {'nc','nc'};
    paths.diag.O2.var   = {'o2','oxygen'};
    paths.diag.O2.zvar  = {'depth','Depth'};
    paths.diag.O2.dim   = {'xyzt','xyz'};
    paths.diag.O2.lon   = {'lon','lon'};
    paths.diag.O2.lat   = {'lat','lat'};
    paths.diag.O2.name  = {'WOA-18','GLODAPv2'};
    paths.diag.O2.units = {'mmol O$_2$ m$^{-3}$','mmol O$_2$ m$^{-3}$'};

    % nitrate 
    paths.diag.NO3.file  = {'/data/project3/data/woa18/nitrate/1p0/no3_woa18_clim.nc',...
                            '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.NO3.nc'};
    paths.diag.NO3.type  = {'nc','nc'};
    paths.diag.NO3.var   = {'no3','NO3'};
    paths.diag.NO3.zvar  = {'depth','Depth'};
    paths.diag.NO3.dim   = {'xyzt','xyz'};
    paths.diag.NO3.lon   = {'lon','lon'};
    paths.diag.NO3.lat   = {'lat','lat'};
    paths.diag.NO3.name  = {'WOA-18','GLODAPv2'};
    paths.diag.NO3.units = {'mmol N m$^{-3}$','mmol N m$^{-3}$'};
    
    % Copy for NOx (NO3 is typically reported via (NO3+NO2))
    paths.diag.NOX = paths.diag.NO3;

    % phosphate
    paths.diag.PO4.file  = {'/data/project3/data/woa18/phosphate/1p0/po4_woa18_clim.nc',...
                            '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.PO4.nc'};
    paths.diag.PO4.type  = {'nc','nc'};
    paths.diag.PO4.var   = {'po4','PO4'};
    paths.diag.PO4.zvar  = {'depth','Depth'};
    paths.diag.PO4.dim   = {'xyzt','xyz'};
    paths.diag.PO4.lon   = {'lon','lon'};
    paths.diag.PO4.lat   = {'lat','lat'};
    paths.diag.PO4.name  = {'WOA-18','GLODAPv2'};
    paths.diag.PO4.units = {'mmol P m$^{-3}$','mmol P m$^{-3}$'};

    % N* (only 1.00 available)
    paths.diag.nstar.file  = {'/data/project3/data/woa18/nstar/1p0/nstar_woa18_clim.nc'};
    paths.diag.nstar.type  = {'nc'};
    paths.diag.nstar.var   = {'nstar'};
    paths.diag.nstar.zvar  = {'depth'};
    paths.diag.nstar.dim   = {'xyzt'};
    paths.diag.nstar.lon   = {'lon'};
    paths.diag.nstar.lat   = {'lat'};
    paths.diag.nstar.name  = {'WOA-18'};
    paths.diag.nstar.units = {'mmol N m$^{-3}$'};
                 
    % silicate 
    paths.diag.SiO3.file  = {'/data/project3/data/woa18/silicate/1p0/si_woa18_clim.nc',...
                            '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.silicate.nc'};
    paths.diag.SiO3.type  = {'nc','nc'};
    paths.diag.SiO3.var   = {'si','silicate'};
    paths.diag.SiO3.zvar  = {'depth','Depth'};
    paths.diag.SiO3.dim   = {'xyzt','xyz'};
    paths.diag.SiO3.lon   = {'lon','lon'};
    paths.diag.SiO3.lat   = {'lat','lat'};
    paths.diag.SiO3.name  = {'WOA-18','GLODAPv2'};
    paths.diag.SiO3.units = {'mmol Si m$^{-3}$','mmol Si m$^{-3}$'};

    % N2O 
    paths.diag.N2O.file   = {'/data/project1/demccoy/ROMS/validation/n2o/n2o_NN_format.mat'};
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
    paths.diag.NO2.file   = {'/data/project2/yangsi/analysis/no2Interior/processed/no2predRF_05-27-2020.nc'};
    paths.diag.NO2.type   = {'nc'};
    paths.diag.NO2.var    = {'no2'};
    paths.diag.NO2.zvar   = {'depth'};
    paths.diag.NO2.dim    = {'xyzt'};
    paths.diag.NO2.lon    = {'lon'};
    paths.diag.NO2.lat    = {'lat'};
    paths.diag.NO2.name   = {'Machine Learning Estimate (Clements et al.)'};
    paths.diag.NO2.units  = {'mmol N m$^{-3}$'};

    % Chl 
    paths.diag.SFC_CHL.file  = {'/data/project1/data/MODIS-Aqua/backup/res0p25/chl_clim_0p25.nc'};
    paths.diag.SFC_CHL.type  = {'nc'};
    paths.diag.SFC_CHL.dim   = {'xyt'};
    paths.diag.SFC_CHL.var   = {'chl'};
    paths.diag.SFC_CHL.zvar  = {[]};
    paths.diag.SFC_CHL.lon   = {'lon'};
    paths.diag.SFC_CHL.lat   = {'lat'};
    paths.diag.SFC_CHL.name  = {'MODIS-Aqua Chlorophyll'};
    paths.diag.SFC_CHL.units = {'mg C m$^{-3}$'};
        
    % NPP products (VGPM, CBPM, CAFE)    
    paths.diag.NPP.file  = {'/data/project2/yangsi/analysis/NPPcode/std-VGPM/std-VGPMnpp_MODIS_clim2002-2018.nc',...
                            '/data/project2/yangsi/analysis/NPPcode/CbPM2/CbPM2npp_MODIS_clim2002-2018.nc',...
                            '/data/project1/data/MODIS-Aqua/CAFE_NPP/climatology/nppclim_CAFE_MODIS_9km.nc'}; 
    paths.diag.NPP.type  = {'nc','nc','nc'};
    paths.diag.NPP.var   = {'npp','npp','npp'};
    paths.diag.NPP.zvar  = {[],[],[]};
    paths.diag.NPP.dim   = {'xyt','xyt','xyt'};
    paths.diag.NPP.lon   = {'lon','lon','lon'};
    paths.diag.NPP.lat   = {'lat','lat','lat'};
    paths.diag.NPP.name  = {'NPP-VGPM','NPP-CBPM','NPP-CAFE'};
    paths.diag.NPP.units = {'mg C m$^{-2}$ d$^{-1}$','mg C m$^{-2}$ d$^{-1}$','mg C m$^{-2}$ d$^{-1}$'};

    % MLD products (WOCE/NODC/ARGO or ARGO only);
    paths.diag.MLD.file  = {'/data/project1/demccoy/ROMS/validation/MLD/Argo_mixedlayers_monthlyclim_12112019.nc',...
                            '/data/project1/demccoy/ROMS/validation/MLD/mld_DR003_c1m_reg2.0.nc'};
    paths.diag.MLD.type  = {'nc','nc'};
    paths.diag.MLD.var   = {'mld_da_mean','mld'};
    paths.diag.MLD.zvar  = {[],[]};
    paths.diag.MLD.dim   = {'txy','xyt'};
    paths.diag.MLD.lon   = {'lon','lon'};
    paths.diag.MLD.lat   = {'lat','lat'};
    paths.diag.MLD.name  = {'Argo Mixed Layer Depth','IFREMER Mixed Layer Depth'}; 
    paths.diag.MLD.units = {'m','m'};

    % POC_FLUX_IN
    paths.diag.POC_FLUX_IN.file   = {'/data/project1/demccoy/ROMS/validation/POC_FLUX_IN/clements_100m_flux.mat',...
                                     '/data/project1/data/particle_flux/Clements_2023/Euphotic_Export_2023.nc',...
                                     '/data/project1/data/particle_flux/Siegel_2014/TotEZ_Annual_Monthly_July2020.mat',...
                                     '/data/project1/data/particle_flux/Nowicki_2022/biopump_model_output.nc',...
                                     '/data/project1/data/particle_flux/Dunne_2005/Dunne_1deg_DM.mat'};
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
    paths.diag.Fe.file   = {'/data/project1/data/Iron_data_AP/Tagliabue_data/Monthly_dFe_V2.nc'};
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
    paths.diag.FG_N2O.file   = {'/data/project1/demccoy/ROMS/validation/n2o/yang_n2o_flux.mat'};
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
    paths.diag.Alk.file   = {'/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TAlk.nc'};
    paths.diag.Alk.type   = {'nc'}; 
    paths.diag.Alk.var    = {'TAlk'}; 
    paths.diag.Alk.zvar   = {'Depth'};
    paths.diag.Alk.dim    = {'xyz'};
    paths.diag.Alk.lon    = {'lon'};
    paths.diag.Alk.lat    = {'lat'};
    paths.diag.Alk.name   = {'GLODAPv2'};
    paths.diag.Alk.units  = {'mmol m$^{-3}$'};  
    
    % DIC via GLODAPv2
    paths.diag.DIC.file   = {'/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc'};
    paths.diag.DIC.type   = {'nc'}; 
    paths.diag.DIC.var    = {'TCO2'}; 
    paths.diag.DIC.zvar   = {'Depth'};
    paths.diag.DIC.dim    = {'xyz'};
    paths.diag.DIC.lon    = {'lon'};
    paths.diag.DIC.lat    = {'lat'};
    paths.diag.DIC.name   = {'GLODAPv2'};
    paths.diag.DIC.units  = {'mmol m$^{-3}$'};  

    % OMZ thickness
    paths.diag.OMZ.file  = {'/data/project1/demccoy/ROMS/validation/OMZ/WOA18_OMZ_thickness.mat',...
                            '/data/project1/demccoy/ROMS/validation/OMZ/Bianchi2012_OMZ_thickness.mat'};
    paths.diag.OMZ.type  = {'mat','mat'};
    paths.diag.OMZ.var   = {'OMZ','OMZ'};
    paths.diag.OMZ.zvar  = {'omzthresh','omzthresh'}; % trick into thinking thresholds are depths
    paths.diag.OMZ.dim   = {'xyz','xyz'};
    paths.diag.OMZ.lon   = {'lon','lon'};
    paths.diag.OMZ.lat   = {'lat','lat'};
    paths.diag.OMZ.name  = {'WOA-18','Bianchi et al. (2012)'};
    paths.diag.OMZ.units = {'m','m'};
end % end method getDiagPaths
