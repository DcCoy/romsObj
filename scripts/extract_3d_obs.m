function extract_3d_obs(obj,fname,plotsDir,rootPath)
% -------------------------------------------------
% Script to plot gridded observations vs ROMS
%
% Usage: extract_3d_obs(obj,fname)
%
% Inputs:
%	obj   = romsObj object
%   fname = filename and path to gridded 3D product
% -------------------------------------------------
warning off
% Addpath to romsOpt
addpath(rootPath); romsOpt;

% Load data
load(fname);
obs = gout3d; clear gout3d;
[grd.lon,grd.lat] = meshgrid(obs.lon,obs.lat);

% Plot options
vars    = {'o2','no3','nh4','no2','nh4vno2','no2vnh4'};
name    = {'oxygen','nitrate','ammonium','nitrite','nitrammn','ammnnitr'};
tits    = {'O$_2$','NO$^{-}_3$','NH$^{+}_4$','NO$^{-}_2$','NH$^{+}_4$ / NO$^{-}_2$','NO$^{-}_2$ / NH$^{+}_4$'};
units   = {'mmol O$_2$ m$^{-3}$','mmol NO$^{-}_3$ m$^{-3}$','mmol NH$^{+}_4$ m$^{-3}$','mmol NO$^{-}_2$ m$^{-3}$','ratio (mol/mol)','ratio (mol/mol)'};
lims    = {linspace(0,300,128),linspace(0,50,128),linspace(0,5,128),linspace(0,8,128)};
lims{3} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1]);
lims{4} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1]);
lims{5} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1 1e1:1e1:1e2 1e2:1e2:1e3]);
lims{6} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1 1e1:1e1:1e2 1e2:1e2:1e3]);
cmaps   = {cmocean('-ice',127),cmocean('tempo',127),cmocean('tempo',127),cmocean('tempo',127),cmocean('balance',127),cmocean('balance',127)}; 
deps = 1:11;
depth = obs.depth(deps); 

% Calculate ratios
obs.nh4vno2.median = obs.nh4.median ./ obs.no2.median;
obs.nh4vno2.median(obs.nh4vno2.median==Inf) = NaN;
obs.nh4vno2.median(obs.nh4vno2.median==0)   = NaN;
obs.no2vnh4.median = obs.no2.median ./ obs.nh4.median;
obs.no2vnh4.median(obs.nh4vno2.median==Inf) = NaN;
obs.no2vnh4.median(obs.no2vnh4.median==0)   = NaN;

% Vars loop
for v = 1:length(vars)
	% Go through depths
	for d = deps 
		fname  = [plotsDir,sprintf('%02d',v),'_gridded_obs_',name{v},'_',num2str(d)];
		% Get quickMap
		fig = quickMap(obj);
		hold on
		% Extract data for depth
		tmpdat = squeeze(obs.(vars{v}).median(d,:,:));
		if v > 2
			tmpdat = log10(tmpdat);	
		end
		% Get pcolor
		warning off
		m_pcolor(grd.lon,grd.lat,real(tmpdat)); shading flat;
		m_coast('patch',coastcolor,'edgecolor','k'); drawnow
		if obj.grid.nx > obj.grid.ny
			cb = colorbar('location','southoutside');
		else
			cb = colorbar('location','eastoutside');
		end
		cb.FontSize = fontsize;
		if v < 3
			caxis([lims{v}(1) lims{v}(end)]);
		else
			caxis([log10(lims{v}(1)) log10(lims{v}(end))]);
			cb.Ticks = log10(lims{v});
			str = [];
			for i = 1:length(cb.Ticks)
				if ismember(i,1:9:length(lims{v}));
					str{1,i} = num2str(lims{v}(i));
				else
					str{1,i} = ' ';
				end
			end
			cb.TickLabels = str;
		end
		set(gca,'ColorMap',cmaps{v});
		ylabel(cb,units{v},'Interpreter','Latex');
		hold on
		m_plot(obj.grid.polygon(:,1),obj.grid.polygon(:,2),'k','linewidth',1);
		title([tits{v},': ',num2str(obs.depth(d)),'m'],'Interpreter','Latex');
		set(gca,'FontSize',fontsize);
		export_fig(figsFormat,[fname],figsQuality);	
		if v == 1 & d == 1
			ax    = gca;
			axpos = ax.Position;
			cbpos = cb.Position; 
			save(['gridded_positions.mat'],'axpos','cbpos','depth');
		end
		close all
	end
end
warning on
