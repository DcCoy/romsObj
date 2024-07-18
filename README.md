# romsObj
romsObj is Matlab class developed at UCLA to load and process ROMS data and make various plots, comparison figures against validation data, and many other actions. It requires MATLAB and command line NCO tools (ncread, ncinfo, ncdump, etc). Contact Daniel McCoy (dmccoy@carnegiescience.edu) for assistance.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

## Updates
* 07/17/2024 -- First commit of romsObj, romsDiag, romsComp, and romsOpt 

## Getting started
#### Update romsOpt.m to set your correct paths
	 - set paths to root, simulation, and run directories
	 - set paths to scripts directory
	 - set paths to validation products
	 - set path to temporary figure folder
	 - set default figure options

#### Create a directory in 'simPath' that corresponds to your simulation configuration (simName)
	- in MATLAB:
	>> simName = 'pacmed_0p25';
    >> mkdir([simPath,simName])

#### Create a subdirectory for your grid file, then copy (or create a symbolic link) your grid file there
	- in MATLAB:
	>> mkdir([simPath,simName,'/grid'])
	>> cd([simPath,simName,'/grid'])
	>> cmd=['ln -s /path/to/grid/pacmed_0p25_grid_file.nc .']; system(cmd)

#### Create a subdirectory for your specific run (runName)
	- in MATLAB:
	>> runName = 'test';
	>> mkdir([simPath,simName,runName]) 

#### Create the following subdirectories in your new directory 
	- in MATLAB:
	>> cd([[simPath,simName,runName])
	>> mkdir bgc   bgc/avg   bgc/his
	>> mkdir phy   phy/avg   phy/his
	>> mkdir dia   dia/avg   dia/his
	>> mkdir flux flux/avg  flux/his
	>> mkdir flx   flx/avg   flx/his

	- bgc corresponds to biogeochemical output (O2, NO3, etc.)
	- phy corresponds to physical output (temp, salt, etc.)
	- dia corresponds to diagnostic output (O2_CONSUMPTION, etc.)
	- flux corresponds to advective flux output
	- flx corresponds to surface flux output (shflux, etc.)

#### Copy or create symbolic links to your output in these directories
	- averaged outputs live in their respective 'avg' directories
	- history snapshots live in their respective 'his' directories 
	- NOTE: if ROMS data is NOT in separate output files (as in older ROMS versions), you only need to create 1
	  link in the 'phy/avg' and/or 'phy/his' directory

## Using romsObj
#### Create an empty object
	- in MATLAB:
	>> obj = romsObj;

#### Initialize the object with your configuration and runName
	- in MATLAB:
	>> obj = initROMS(obj,'pacmed_0p25','runName','test')

#### Object structure
    -  obj.info      = struct containing simulation and variable information (via romsInfo)
    -  obj.grid      = struct containing grid data and dimensions (via loadGrid)
    -  obj.slice     = struct containing slice coordinates (via sliceROMS/sliceDiag)
    -  obj.profile   = struct containing profile coordinates (via getProfile)
    -  obj.data      = struct containing ROMS output (via loadData, computeVar, sliceROMS)
    -  obj.diag      = struct containing validation data for comparions (via loadDiag, sliceDiag)
    -  obj.paths     = struct containing paths to data, directories, diagnostic data, etc (via initROMS)
    -  obj.constants = struct containing useful constants and conversions

## Code structure 
#### Methods for class romsObj
    -  initROMS   = Initialization method: gathers paths and coordinate variables 
    -  romsInfo   = Obtains info from sample roms files in run directory (called during initROMS) 
    -  loadGrid   = Loads 2D grid information into obj.grid (called during initROMS)
	-  loadDepth  = Load Z-grid info (z_r,z_w), which are output-file dependent.
	-  loadData   = Main method to load ROMS data from specific files
	-  computeVar = Computes additional fields like wind stress, Okubo-Weiss, AOU, etc. 
	-  sliceROMS  = Takes depth slice of ROMS along a given latitude, longitude, xi-index, or eta-index
	-  zslice     = Slices ROMS along a given depth
	-  ipslice    = Slices ROMS along a given isopycnal
	-  intVar     = Vertically integrates 3D variables
	-  getProfile = Loads profile data at the nearest lon/lat point
	-  dispVars   = Lists available output variables in command window
	-  gridView   = Plots the xi/eta indices of a grid file in a map figure
	-  regionView = Plots a ROMS grid map along with a regional grid (if requested during initROMS)
	-  Dist2Coast = Calculates the distance of each grid point to the nearest masked point
	-  mapPlot    = Plots 2D fields onto a map projection
	-  mapCmp     = Takes 2 2D fields, plots them, and produces a difference plot (3 figures)
	-  axPlot     = Similar to mapPlot, but attempts to place map inside a pre-defined axis
	-  quickMap   = Similar to mapPlot, but produces an empty map
	-  slicePlot  = Plots transect variables obtained from sliceROMS or sliceDiag
	-  sliceDiag  = Takes depth transect of diagnostic data along a given latitude or longitude
	-  sliceCmp   = Similar to mapCmp, but used to compare transects
	-  loadDiag   = Loads and interpolates monthly 2D validation data to the ROMS grid
	-  clearROMS  = Clears loaded data

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

