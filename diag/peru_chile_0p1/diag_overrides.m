%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permanent overrides for Peru Chile configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transects
plots(2).choice = {  'lat', 'lon', 'lon'};
plots(2).degs   = [  0       255    272 ];
plots(2).xlims =  [255 278;-40  0;-40  0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User overrides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set user-specific overrides here
%
% e.g., to only plot temperature zslices...
% plots(1).opt = [1 0 0 0 0 0 0 0 0 0 0 0 0 0]; <-- turn off plotting of all other variables
% 
% e.g., to adjust depths to zslice 
% plots(1).zdeps = [100]; <-- only zslice at 100m
%
% e.g., change the latitudes or longitudes to slice at 
% plots(2).lats = [-10 10];  <-- plot zonal sections at -10S & 10N
% plots(3).lons = [230 250]; <-- plot meridional sections at 230 & 250
