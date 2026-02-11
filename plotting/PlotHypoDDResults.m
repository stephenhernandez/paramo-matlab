clear; close all;
cd ~/research/now/ecuador/hypodd/
regionName = "ecuador";
boundaryBox = getRegionSpatialDimensions(regionName);

phase_data_file = strcat('phase_',regionName,'.dat'); %actual PHASE _data_!
unix(['bash ~/scripts/separateCatalogs.sh ',char(phase_data_file)]);
unix('wc reloc_cat.txt');
dataDir = '~/igdata';
load(fullfile(dataDir,'ec_boundaries.mat'))
S = shaperead(fullfile(dataDir,'volcanes/volcanes_2.shp'));
load('~/igdata/soam_noec','lat_noec','lon_noec');

%%
[tReloc,eqlatReloc,eqlonReloc,eqdepthReloc,...
    eqmagReloc,idReloc] = readHypoDD(regionName);
R = readtable('orig_catalog_events_also_relocated.txt');
idOrig = R.Var1;
tOrig = datetime(R.Var2,R.Var3,R.Var4,R.Var5,R.Var6,R.Var7);
eqlatOrig = R.Var8;
eqlonOrig = R.Var9;
eqdepthOrig = R.Var10;
eqmagOrig = R.Var11;

[lia,locb]  = ismember(idReloc,idOrig);
if sum(lia) ~= length(locb(lia))
    idOrig = idOrig(locb(lia));
    tOrig = tOrig(locb(lia));
    eqlatOrig = eqlatOrig(locb(lia));
    eqlonOrig = eqlonOrig(locb(lia));
    eqdepthOrig = eqdepthOrig(locb(lia));
    eqmagOrig = eqmagOrig(locb(lia));
end

%%
fI = eqdepthReloc >= -30;
magFact = 2;

figure('units','normalized','outerposition',[0.03 0.05 0.94 0.9]);
nSubplots = 2;
tiledlayout(1,nSubplots, 'Padding', 'compact', 'TileSpacing', 'compact');
ax(1) = nexttile;
SS = scatter(eqlonOrig(fI),eqlatOrig(fI),magFact*exp(eqmagOrig(fI)),eqdepthOrig(fI),'filled');
SS.MarkerFaceAlpha = 0.5;
SS.MarkerEdgeColor = 'k';
SS.MarkerEdgeAlpha = 0.5;
axis equal; grid on;
colorbar;

set(ax(1),'ColorScale','log');
zoom on;
colormap turbo;
hold on;
plot(lonEC,latEC,'k','linewidth',3);
%plot(lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
axis(boundaryBox);
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color',[0.5 0.5 0.5],'linewidth',1);
plot([S.X],[S.Y],'k--','linewidth',1);
%geoshow(ax(1),'~/igdata/PobEc_cpv2021.shp');
clim([1 300]);

%%
ax(2) = nexttile;
SS = scatter(eqlonReloc(fI),eqlatReloc(fI),magFact*exp(eqmagReloc(fI)),eqdepthReloc(fI),'filled');
SS.MarkerFaceAlpha = 0.5;
SS.MarkerEdgeColor = 'k';
SS.MarkerEdgeAlpha = 0.5;
axis equal; grid on;
colorbar;

set(ax(2),'ColorScale','log');
zoom on;
colormap turbo;
hold on;
plot(lonEC,latEC,'k','linewidth',3);
%plot(lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
axis(boundaryBox);
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color',[0.5 0.5 0.5],'linewidth',1);
plot([S.X],[S.Y],'k--','linewidth',1);
%geoshow(ax(2),'~/igdata/PobEc_cpv2021.shp');
clim([1 300]);

%%
linkaxes(ax,'xy');