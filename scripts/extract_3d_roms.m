function extract_3d_roms(obj,file,fname,plotsDir,rootPath)
% Script to plot ROMS like we do validation data 
warning off
% Addpath to romsOpt
addpath(rootPath); romsOpt;

% Load ax positions
load('gridded_positions.mat');

% Plot options
vars    = {'O2','NO3','NH4','NO2','NH4vNO2','NO2vNH4'};
tits    = {'O$_2$','NO$^{-}_3$','NH$^{+}_4$','NO$^{-}_2$','NH$^{+}_4$ / NO$^{-}_2$','NO$^{-}_2$ / NH$^{+}_4$'};
name    = {'oxygen','nitrate','ammonium','nitrite','nitrammn','ammnnitr'};
units   = {'mmol O$_2$ m$^{-3}$','mmol NO$^{-}_3$ m$^{-3}$','mmol NH$^{+}_4$ m$^{-3}$','mmol NO$^{-}_2$ m$^{-3}$','ratio (mol/mol)','ratio (mol/mol)'};
lims    = {linspace(0,300,128),linspace(0,50,128),linspace(0,5,128),linspace(0,8,128)};
lims{3} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1]);
lims{4} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1]);
lims{5} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1 1e1:1e1:1e2 1e2:1e2:1e3]);
lims{6} = unique([1e-3:1e-3:1e-2 1e-2:1e-2:1e-1 1e-1:1e-1:1e0 1e0:1e0:1e1 1e1:1e1:1e2 1e2:1e2:1e3]);
cmaps   = {cmocean('-ice',127),cmocean('tempo',127),cmocean('tempo',127),cmocean('tempo',127),cmocean('balance',127),cmocean('balance',127)}; 
zdeps = load('gridded_positions.mat','depth');
zdeps = zdeps.depth;
cmd = 'rm gridded_positions.mat'; system(cmd);

% Vars loop
for v = 1:length(vars)
	obj = zslice(obj,vars(v),zdeps,file);
	for d = 1:length(zdeps) 
		fname  = [plotsDir,sprintf('%02d',v),'_gridded_roms_',name{v},'_',num2str(d)];
		% Get quickMap
		fig = quickMap(obj);
		hold on
		% Extract data for depth
		tmpdat = nanmean(squeeze(obj.data.avg.(vars{v}).slice(:,:,d,:)),3);
		tmpdat(tmpdat<0) = 0;
		if v < 3
			m_contourf(obj.grid.lon_rho,obj.grid.lat_rho,tmpdat,lims{v},'linestyle','none');
			if obj.grid.nx > obj.grid.ny
				cb = colorbar('location','southoutside');
			else
				cb = colorbar('location','eastoutside');
			end
			caxis([lims{v}(1) lims{v}(end)]);
			set(gca,'ColorMap',cmaps{v});
		else
			tmpdat(tmpdat<1e-3) = 1e-3;
			tmpdat(tmpdat>1e3)  = 1e3;
			m_contourf(obj.grid.lon_rho,obj.grid.lat_rho,log10(tmpdat),log10(lims{v}),'linestyle','none');
			if obj.grid.nx > obj.grid.ny
				cb = colorbar('location','southoutside');
			else
				cb = colorbar('location','eastoutside');
			end            
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
            set(gca,'ColorMap',cmaps{v});
		end
		cb.FontSize = fontsize;
		ylabel(cb,units{v},'Interpreter','Latex');
		title([tits{v},': ',num2str(zdeps(d)),'m'],'Interpreter','Latex');
		set(gca,'FontSize',fontsize);
		ax = gca;
		ax.Position = axpos;
		cb.Position = cbpos;
		export_fig(figsFormat,[fname],figsQuality);
		close all
	end
end
warning on
