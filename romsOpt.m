% Set your romsObj options here

% PATHS
rootPath      = '/data/project1/demccoy/romsObj/';                  % where romsObj, romsDiag, romsComp, romsOpt live
simPath       = '/data/project2/model_output/';                     % where ROMS output is stored (or linked)
scriptsPath   = '/data/project1/demccoy/romsObj/scripts/';          % where additional matlab scripts live 
valiPath      = '/data/project1/demccoy/ROMS/validation/products/'; % where validation data lives
tmpfigsPath   = '/data/project1/demccoy/tmpfigs/';                  % where temporary figures will be stored 
diagPath      = '/data/project1/demccoy/romsObj/diag/';             % where diagnostic scripts, options, and plots live
compPath      = '/data/project1/demccoy/romsObj/comp/';             % where ROMS-to-ROMS comparison scripts, options, and plots live

% ADD PATHS
addpath(rootPath);
addpath(simPath);
addpath(scriptsPath);
addpath([scriptsPath,'m_map']);
addpath(valiPath);
addpath(tmpfigsPath);
addpath(diagPath);
addpath(compPath);

% DEFAULT FIGURE OPTIONS
figsFormat    = 'png';            % format of output figures (jpg, png, pdf)
figsQuality   = '-m5';            % factor to magnify onscreen figure pixel dimensions (export_fig) for higher quality plots
figVisible    = 'off';            % ('on') to make figure windows visible (not required for printing) 
background    = rgb('LightGray'); % color to be used in backgrounds of plots 
coastcolor    = rgb('DimGray');   % color to be used for NaN points (coasts, bathymetery, etc) 
fontsize      = 6;                % fontsize in figures
figtype       = 'mfig';           % 'sfig' = 9cm wide, 'mfig' = 14cm, 'lfig' = 19cm   
coast         = 'coast';          % choose from: coast, crude, low, high, intermediate, or full (map coast resolution)  
ticks         = 1;                % (1) to show lon/lat ticks on x/y axes for mapPlot figures
polygon       = 1;                % (1) to show ROMS model boundary in figures

% LATEX-FRIENDLY LABELS (LEAVE EMPTY TO DEFAULT TO NETCDF UNITS) 
strings_to_replace = {...
	'kilogram meter-3','kg m$^{-3}$';...
	'kg/m^2/s'        ,'kg m$^{-2}$ s$^{-1}$';...
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
	'W/m^2'           ,'W m$^{-2}$';...
	'mmol N/m3/s'     ,'mmol N m$^{-3}$ s$^{-1}$';...
	'mmol N2O/m3/s'   ,'mmol N$_2$O m$^{-3}$ s$^{-1}$';...
	'mmol N2/m3/s'    ,'mmol N$_2$ m$^{-3}$ s$^{-1}$';...
	'mmol/m3/s'       ,'mmol m$^{-3}$ s$^{-1}$';...
	'mmol/m2/s'       ,'mmol m$^{-2}$ s$^{-1}$';...
	'm^2/s^2'         ,'m$^{2}$ s$^{-2}$';...
	'm/s'             ,'m s$^{-1}$';...
	'mmol/m^3'        ,'mmol m$^{-3}$';...
	'mmol/m^2'        ,'mmol m$^{-2}$';...
	'mmol/m^3/s'      ,'mmol m$^{-3}$ s$^{-1}$';...
	'mmol/m^2/s'      ,'mmol m$^{-2}$ s$^{-1}$';...
	'meter2 second-1' ,'m$^{2}$ s$^{-1}$';...
	'meter second-1'  ,'m s$^{-1}$';...
	'1/y'             ,'yr$^{-1}$';...
					  };      

% BIOGEOCHEMICAL BUDGET SETTINGS (see getBudg)
%
%     budget.(var).tits   = LaTex-friendly title
%     budget.(var).units  = 3D rates units, 2D rates unit (vertical integration), 1D rates unit (volume integration)
%     budget.(var).rates  = 3D diagnostics list used in (var) tracer sources-minus-sinks
%     budget.(var).rtits  = LaTex-friendly 3D diagnostics titles
%     budget.(var).smseq  = 'Sources-minus-sinks' equation, corresponds to 3D diagnostics list
%     budget.(var).fluxes = 2D diagnostics list used in (var) tracer sources-minus-sinks 
%     budget.(var).lvls   = 'sfc' or 'sed', for surface or sediment interface fluxes
%     budget.(var).ftits  = LaTex-friendly 2D diagnostics titles
%
% NO3 (BEC-NitrOMZ)
budget.NO3.tits   = {'NO$^{-}_3$'};
budget.NO3.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
budget.NO3.rates  = {'NITROX','DENITRIF1','SP_NO3_UPTAKE','DIAT_NO3_UPTAKE','DIAZ_NO3_UPTAKE'};
budget.NO3.rtits  = {'NO$^{-}_2$ oxidation','NO$^{-}_3$ reduction',...
					 'NO$^{-}_3$ uptake via small phyto','NO$^{-}_3$ uptake via diatoms',...
					 'NO$^{-}_3$ uptake via diazotrophs'};
budget.NO3.smseq  = [(1) (-1) (-1) (-1) (-1)];
budget.NO3.fluxes = {'SED_DENITRIF'};
budget.NO3.lvls   = {'sed'};
budget.NO3.ftits  = {'Sediment denitrification'};
budget.NO3.feq    = [(1)];
% NO2 (BEC-NitrOMZ)
budget.NO2.tits   = {'NO$^{-}_2$'};
budget.NO2.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
budget.NO2.rates  = {'AMMOX','N2OAMMOX','NITROX','ANAMMOX','DENITRIF1','DENITRIF2',...
		   			 'SP_NO2_UPTAKE','DIAT_NO2_UPTAKE','DIAZ_NO2_UPTAKE'};
budget.NO2.rtits  = {'NH$^{+}_4$ oxidation','NH$^{+}_4$ oxidation to N$_2$O',...
					 'NO$^{-}_2$ oxidation','Anammox','NO$^{-}_3$ reduction',...
					 'NO$^{-}_2$ reduction',...
					 'NO$^{-}_2$ uptake via small phyto','NO$^{-}_2$ uptake via diatoms',...
					 'NO$^{-}_2$ uptake via diazotrophs'};
budget.NO2.smseq  = [(1) (-2) (-1) (-1) (1) (-1) (-1) (-1) (-1)];    
budget.NO2.fluxes = {[]};
budget.NO2.lvls   = {[]};
budget.NO2.ftits  = {[]};
budget.NO2.feq    = [];
% NH4 (BEC-NitrOMZ)
param.DONrefract  = 0.0115;% Fraction of DON to refractory pool
param.Q           = 0.137; % N/C ratio
budget.NH4.tits   = {'NH$^{+}_4$'};
budget.NH4.units  = {'mmol N m$^{-3}$ s$^{-1}$','mmol N m$^{-2}$ s$^{-1}$','mmol N s$^{-1}$'}; 
budget.NH4.rates  = {'SP_NH4_UPTAKE','DIAT_NH4_UPTAKE','DIAZ_NH4_UPTAKE','AMMOX','ANAMMOX',...
		   			 'DON_REMIN','DONr_REMIN','ZOO_LOSS_DIC','SP_LOSS_DIC','DIAT_LOSS_DIC',...
		   			 'DIAZ_LOSS_DIC','SP_GRAZE_DIC','DIAT_GRAZE_DIC','DIAZ_GRAZE_DIC',...
		   			 'POC_REMIN'}; 
budget.NH4.rtits  = {'NH$^{+}_4$ uptake via small phyto','NH$^{+}_4$ uptake via diatoms',...
		   			 'NH$^{+}_4$ uptake via diazotrophs','NH$^{+}_4$ oxidation',...
		   			 'Anammox','DON remineralization','Refractory DON remineralization',...
		   			 'Zooplankton mortality','Small phyto mortality',...
		   			 'Diatom mortality','Diazotroph mortality','Small phyto grazing loss',...
		   			 'Diatom grazing loss','Diatotroph grazing loss','POC remineralization'};
budget.NH4.smseq  = [(-1) (-1) (-1) (-1) (-1) (1) (1) (param.Q) (param.Q) (param.Q) (param.Q) ...
                     (param.Q) (param.Q) (param.Q) (param.Q.*(1-param.DONrefract))];
budget.NH4.fluxes = {[]};
budget.NH4.lvls   = {[]};    
budget.NH4.ftits  = {[]};
budget.NH4.feq    = [];
% N2 (BEC-NitrOMZ)
budget.N2.tits   = {'N$_2$'};
budget.N2.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
budget.N2.rates  = {'DENITRIF3','ANAMMOX'};
budget.N2.rtits  = {'N$_2$O reduction','Anammox'};
budget.N2.smseq  = [(1) (1)]; 
budget.N2.fluxes = {'FG_N2','SED_DENITRIF'};
budget.N2.lvls   = {'sfc','sed'};
budget.N2.ftits  = {'$\Phi$ N$_2$'};
budget.N2.feq    = [(1) (0.5)];
% N2O (BEC-NitrOMZ)
budget.N2O.tits   = {'N$_2$O'};
budget.N2O.units  = {'mmol N$_2$ m$^{-3}$ s$^{-1}$','mmol N$_2$ m$^{-2}$ s$^{-1}$','mmol N$_2$ s$^{-1}$'}; 
budget.N2O.rates  = {'DENITRIF2','N2OAMMOX','DENITRIF3'};
budget.N2O.rtits  = {'NO$^{-}_2$ reduction','NH$^{+}_4$ oxidation to N$_2$O','N$_2$O reduction'};
budget.N2O.smseq  = [(0.5) (1) (-1)]; 
budget.N2O.fluxes = {'FG_N2O'};
budget.N2O.lvls   = {'sfc'};
budget.N2O.ftits  = {'$\Phi$ N$_2$O'};
budget.N2O.feq    = [1];
% O2 (BEC-NitrOMZ)
budget.O2.tits   = {'O$_2$'};
budget.O2.units  = {'mmol O$_2$ m$^{-3}$ s$^{-1}$','mmol O$_2$ m$^{-2}$ s$^{-1}$','mmol O$_2$ s$^{-1}$'};
budget.O2.rates  = {'O2_PRODUCTION','O2_CONSUMPTION'};
budget.O2.rtits  = {'O$_2$ production','O$_2$ consumption'};
budget.O2.smseq  = [(1) (-1)];
budget.O2.fluxes = {'FG_O2'}; 
budget.O2.lvls   = {'sfc'};
budget.O2.ftits  = {'$\Phi$ O$_2$'};
budget.O2.feq    = [(1)];
% DiDA (BEC-HABS)
budget.DiDA.tits   = {'DiDA'};
budget.DiDA.units  = {'mmol m$^{-3}$ s$^{-1}$','mmol m$^{-2}$ s$^{-1}$','mmol s$^{-1}$'};
budget.DiDA.rates  = {'DiDA_PROD','DiDA_LOSS'};
budget.DiDA.rtits  = {'DiDA$_{prod}$','DiDA$_{loss}$'};
budget.DiDA.smseq  = [(1) (-1)];
budget.DiDA.fluxes = {[]};
budget.DiDA.lvls   = {[]};
budget.DiDA.ftits  = {[]};
budget.DiDA.feq    = [];
% DDA (BEC-HABS)
budget.DDA.tits   = {'DDA'};
budget.DDA.units  = {'mmol m$^{-3}$ s$^{-1}$','mmol m$^{-2}$ s$^{-1}$','mmol s$^{-1}$'};
budget.DDA.rates  = {'DDA_PROD','DDA_REMIN'};
budget.DDA.rtits  = {'DDA$_{prod}$','DDA$_{remin}$'};
budget.DDA.smseq  = [(1) (-1)];
budget.DDA.fluxes = {[]};
budget.DDA.lvls   = {[]};
budget.DDA.ftits  = {[]};
budget.DDA.feq    = [];
% zooDA (BEC-HABS)
budget.zooDA.tits   = {'zooDA'};
budget.zooDA.units  = {'mmol m$^{-3}$ s$^{-1}$','mmol m$^{-2}$ s$^{-1}$','mmol s$^{-1}$'};
budget.zooDA.rates  = {'ZooDA_PROD','ZooDA_LOSS'};
budget.zooDA.rtits  = {'ZooDA$_{prod}$','ZooDA$_{loss}$'};
budget.zooDA.smseq  = [(1) (-1)];
budget.zooDA.fluxes = {[]};
budget.zooDA.lvls   = {[]};
budget.zooDA.feq    = [];
