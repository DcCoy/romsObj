# romsObj
romsObj is Matlab class developed at UCLA to load and process ROMS data and make various plots, comparison figures against validation data, and many other actions. It requires MATLAB and command line NCO tools (ncread, ncinfo, ncdump, etc). Contact Daniel McCoy (dmccoy@carnegiescience.edu) for assistance.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Using romsObj](#using-romsObj)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

## Updates
* 07/17/2024 -- First commit of romsObj, romsDiag, romsComp, and romsOpt 
* 07/29/2024 -- Included diagnostic LaTeX file and generation shell script 

## Getting started
#### Update romsOpt.m for your system and settings
     - Set paths for your system / directory structure
     - Set default figure options
     - Update list of LaTeX-friendly labels for ROMS variables (used for output PDFs)
     - (OPTIONAL) Provide variable-specific info for biogeochemical budget routine

#### Create a directory in 'simPath' that corresponds to your simulation configuration (simName)
    - in MATLAB:
    >> simName = 'pacmed_0p25';
    >> mkdir([simPath,simName])

    - Example:
		- /data/project2/model_output/pacmed_0p25/

#### Create a subdirectory for your grid file, then copy (or create a symbolic link) your grid file there
    - in MATLAB:
    >> mkdir([simPath,simName,'/grid'])
    >> cd([simPath,simName,'/grid'])
    >> cmd=['ln -s /path/to/grid/pacmed_0p25_grid_file.nc .']; system(cmd)

    - Example:
		- /data/project2/model_output/pacmed_0p25/grid
		- /data/project2/model_output/pacmed_0p25/grid/pacmed_0p25_grid_file.nc

#### Create a subdirectory for your specific run (runName)
    - in MATLAB:
    >> runName = 'test';
    >> mkdir([simPath,simName,runName]) 

    - Example:
		- /data/project2/model_output/pacmed_0p25/test/

#### Create the following subdirectories in your new directory 
    - bgc corresponds to biogeochemical output (O2, NO3, etc.)
    - phy corresponds to physical output (temp, salt, etc.)
    - dia corresponds to diagnostic output (O2_CONSUMPTION, etc.)
    - flux corresponds to tracer advective flux output
    - flx corresponds to surface flux output (shflux, etc.)

    - in MATLAB:
    >> cd([[simPath,simName,runName])
    >> mkdir bgc   bgc/avg   bgc/his
    >> mkdir phy   phy/avg   phy/his
    >> mkdir dia   dia/avg   dia/his
    >> mkdir flux flux/avg  flux/his
    >> mkdir flx   flx/avg   flx/his

    - Example (here, for bgc output...should be identical for other types):
		- /data/project2/model_output/pacmed_0p25/test/bgc
		- /data/project2/model_output/pacmed_0p25/test/bgc/avg
		- /data/project2/model_output/pacmed_0p25/test/bgc/his


#### Copy or create symbolic links to your output in these directories
    - averaged outputs live in their respective 'avg' directories
    - history snapshots live in their respective 'his' directories 
    - NOTE: if ROMS data is NOT in separate output files (as in older ROMS versions), you only need to create 1
      link in the 'phy/avg' and/or 'phy/his' directory

    - Example (year 2021 annually-averaged files)
		- /data/project2/model_output/pacmed_0p25/test/bgc/avg/pacmed_bgc_avg.20210101113000.nc
		- /data/project2/model_output/pacmed_0p25/test/phy/avg/pacmed_avg.20210101120000.nc
		- /data/project2/model_output/pacmed_0p25/test/dia/avg/pacmed_bgc_dia_avg.20210101120000.nc
		- /data/project2/model_output/pacmed_0p25/test/flx/avg/pacmed_flx_avg.20200201000000.nc

## Using romsObj
#### Create an empty object
    - in MATLAB:
    >> obj = romsObj;

#### Initialize the object with your configuration and runName
    - in MATLAB:
	>> simName = 'pacmed_0p25';
	>> runName = 'test';
    >> obj = initROMS(obj,'pacmed_0p25','test')

#### Object structure
    -  obj.info      = struct containing simulation and variable information (via romsInfo)
    -  obj.grid      = struct containing grid data and dimensions (via loadGrid)
    -  obj.slice     = struct containing slice coordinates (via sliceROMS/sliceDiag)
    -  obj.profile   = struct containing profile coordinates (via getProfile)
    -  obj.budget    = struct containing biogeochemical budget output (via getBudg) 
    -  obj.data      = struct containing ROMS output (via loadData, computeVar, sliceROMS)
    -  obj.diag      = struct containing validation data for comparions (via loadDiag, sliceDiag)
    -  obj.paths     = struct containing paths to data, directories, diagnostic data, etc (via initROMS)
    -  obj.constants = struct containing useful constants and conversions

## Code structure 
#### Initialization methods for class romsObj 
    -  initROMS  = Initialization method: gathers paths and coordinate variables 
    -  romsInfo  = Obtains info from sample roms files in run directory (called during initROMS) 
    -  loadGrid  = Loads 2D grid information into obj.grid (called during initROMS)
    -  clearROMS = Clears data structures

#### Loading methods for class romsObj 
    -  dispVars   = Lists available output variables in command window
    -  loadDepth  = Load Z-grid info (z_r,z_w), which are output-file dependent.
    -  loadData   = Main method to load ROMS data from specific files
    -  loadDiag   = Loads and interpolates monthly 2D validation data to the ROMS grid
    -  computeVar = Computes additional fields like wind stress, Okubo-Weiss, AOU, etc. 

#### Data processing methods for class romsObj 
    -  sliceROMS  = Takes depth slice of ROMS along a given latitude, longitude, xi-index, or eta-index
    -  sliceDiag  = Takes depth transect of diagnostic data along a given latitude or longitude
    -  zslice     = Slices ROMS along a given depth
    -  ipslice    = Slices ROMS along a given isopycnal
    -  intVar     = Interpolate 3D variables both vertically ('int') and over domain volumes ('tot')
    -  getProfile = Loads profile data at the nearest lon/lat point

#### Biogeochemical budget methods for class romsObj 
    -  getBudg       = Main method, calls other routines below 
    -  computeDcDt   = Computes tracer tendencies from history ('his') snapshot output (called in getBudg) 
    -  computeAdv    = Computes tracer advective fluxes from flux ('flux') output (called in getBudg)
    -  computeSMS    = Computes tracer sources and sinks (called in getBudg, settings required in romsOpt) 
    -  computeFluxes = Computes tracer fluxes from surface and sediment interfaces (called in getBudg, settings required in romsOpt) 
    -  computeNet    = Computes budget remainder (should be close to 0, called in getBudg) 
    -  intBudg       = Integrates budget terms vertically ('int') and over domain volumes ('tot') (OPTIONAL in getBudg)
    -  sourcesSinks  = Extracts individual sources and sinks (OPTIONAL in getBudg)
    -  plotIntBudg   = Makes mapPlot(s) of vertically-integrated budget terms (dDdt, adv, sms, air-sea fluxes, etc) 
    -  plotTotBudg   = Makes bar chart of volume-integrated budget terms, useful in checking for closed budgets 

#### Plotting methods for class romsObj 
    -  gridView   = Plots the xi/eta indices of a grid file in a map figure
    -  regionView = Plots a ROMS grid map along with a regional grid (if requested during initROMS)
    -  Dist2Coast = Calculates the distance of each grid point to the nearest masked point
    -  mapPlot    = Plots 2D fields onto a map projection
    -  mapCmp     = Takes 2 2D fields, plots them, and produces a difference plot (3 figures)
    -  axPlot     = Similar to mapPlot, but attempts to place map inside a pre-defined axis
    -  quickMap   = Similar to mapPlot, but produces an empty map
    -  slicePlot  = Plots transect variables obtained from sliceROMS or sliceDiag
    -  sliceCmp   = Similar to mapCmp, but used to compare transects

#### Static methods for class romsObj
    - lon360         = Converts longitude points to    0 to 360 format
    - lon180         = Converts longitude poitns to -180 to 180 format
    - struct2double  = Converts all fields in a structure to double format 
    - struct2single  = Converts all fields in a structure to single format
    - WindStress     = Calculates wind stress and wind stress curl from zonal/meridional wind stress input
    - gc_dist        = Calculates distance between 2 points along a great circle (adapted from J. Molemaker)
    - prclims        = Script to automatically create limits of colorbar from input data
    - grid_to_grid   = Function to grid 2D data from one meshgrid to another (i.e., for validation data)
    - okubo_weiss    = Computes the Okubo-Weiss parameter, relative vorticity, and shear/normal strain components
    - psi2rho        = Converts a field at psi points to the rho points
    - u2rho          = Converts a field at u points to the rho points
    - v2rho          = Converts a field at v points to the rho points
    - spatial_filter = Smooths 2D array data, ignores NaNs.
    - o2_sat         = Calculates O2 @ saturation with respect to SST/SSS
    - n2o_sat        = Calculates N2O @ saturation with respect to SST/SSS
    - parse_pv_pairs = Used to process optional inputs (varargin)
    - piofigs        = 'Progress in Oceanography' figure dimensions. Creates figure width/height based on input
    - pltjpg         = Prints current figure to temporary figure folder (set in romsOpt.m)

#### Additional methods for class romsObj
    - getDiagPaths = Sets paths for individual validation products (set in romsOpt.m, used in loadDiag & sliceDiag)

## Support
Contact Daniel McCoy at the Carnegie Institution for Science (dmccoy@carnegiescience.edu) 

