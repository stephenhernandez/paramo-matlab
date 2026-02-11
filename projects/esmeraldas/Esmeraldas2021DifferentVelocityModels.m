% esmeraldas2022 from different velocity models
clear; close all; clc;

cd ~/research/now/esmeraldas/hypodd/
models = ["esmeraldas2021.m3d.reloc";...
    "esmeraldas2021.mleon.reloc";...
    "esmeraldas2021.leon.reloc";...
    "esmeraldas2021.central.reloc";...
    "esmeraldas2021.macquet.reloc";...
    "esmeraldas2021.soto.reloc";...
    "esmeraldas2021.mfont.reloc"];

titleStrs = ["M3D";"MLEON";"LEON";"CENTRAL";"MACQUET";"SOTO";"MFONT"];
for i = 1:length(models)+1
    if i == length(models)+1
        load esmeraldas2021.orig.txt; %orig_cat.txt;
        t = datetime(esmeraldas2021_orig(:,2:7));
        eqlat = esmeraldas2021_orig(:,8);
        eqlon = esmeraldas2021_orig(:,9);
        eqdepth = esmeraldas2021_orig(:,10);
        eqmag = esmeraldas2021_orig(:,11);
        plot_fig_(eqlat,eqlon,eqdepth,eqmag,sprintf("%s, N=%d\n","IASPEI",length(t)));
        continue;
    end

    [t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD(models(i),true);
    [eqmag,sI] = sort(eqmag,'descend');
    t = t(sI);
    eqlat = eqlat(sI);
    eqlon = eqlon(sI);
    eqdepth = eqdepth(sI);
    id = id(sI);

    refVMdepth = 0;
    eqdepth = eqdepth-refVMdepth;
    plot_fig_(eqlat,eqlon,eqdepth,eqmag,sprintf("%s, N=%d",titleStrs(i),length(t)));
end

function plot_fig_(eqlat,eqlon,eqdepth,eqmag,titleStr)
[minLon,maxLon,minLat,maxLat] = deal(-80,-79.5,0.80,1.6);
colorbarRightShift = 0.8; 
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');

load ~/igdata/ecuador_slab_model.mat
latSlabI = latSlab >= minLat & latSlab <= maxLat;
slabDepth = slabDepth(latSlabI,:);
latSlab = latSlab(latSlabI);
nans = nans(latSlabI,:);

lonSlabI = lonSlab >= minLon & lonSlab <= maxLon;
slabDepth = slabDepth(:,lonSlabI);
lonSlab = lonSlab(lonSlabI);
nans = nans(:,lonSlabI);

[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
demData = double(demData);
lI = lon >= minLon & lon <= maxLon;
demData(~lI,:) = [];
lon = lon(lI);
lI = lat >= minLat & lat <= maxLat;
demData(:,~lI) = [];
lat = lat(lI);

contrast = 1;
ilumAz = -90;
wMark = 2;
elevAngle = 70;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 0.55 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1]);
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [colorbarRightShift,0.55,0.0083333333333333,0.35];

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

geoshow(ax(laxes),'~/igdata/ZonaUrbana/ZonaUrbana.shp');

plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

%% scale ruler
refShifter = (maxLat-minLat)*(0.01/0.3);
scaleRefLat0 = minLat + refShifter; 
scaleRefLon0 = minLon + refShifter;
refEllipse = referenceEllipsoid('wgs84');

[~,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);
[~,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);
[~,lonout20] = reckon(scaleRefLat0,scaleRefLon0,20*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
plot(ax(laxes),[scaleRefLon0 lonout20],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -refShifter/2;
textFontSize = 22;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout5,-textShifter+scaleRefLat2,"5","FontSize",textFontSize);

plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout10,-textShifter+scaleRefLat2,"10","FontSize",textFontSize);

plot(ax(laxes),[lonout20 lonout20],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout20,-textShifter+scaleRefLat2,"20 km","FontSize",textFontSize);

labelFontSize = 24;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);
title(ax(laxes),titleStr,'FontSize',labelFontSize)

%% bathymetry contours
laxes = laxes +1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

contour(ax(laxes),lon,lat,demData',(-(50:100:4000))','-','linewidth',0.1);
cbar12 = colorbar(ax(laxes));
cbar12.Position = [colorbarRightShift,0.55,0.0083333333333333,0.35];
axis(ax(laxes),'equal');
hold(ax(laxes),'on');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
bonemap = flipud(bone(512));
colormap(ax(laxes),bonemap);
ax(laxes).Visible = 'off';
cbar12.Visible = 'off';
hold(ax(laxes),'on');
axis(ax(laxes),'equal');

%%
zInc = 2;
zlevs = 1000*(ceil(min(min(slabDepth)/1000)):zInc:floor(max(max(slabDepth)/1000)));
for iii = 1:length(zlevs)
    Ctmp = contour(lonSlab,latSlab,slabDepth,[zlevs(iii) zlevs(iii)]);
    [~,~,xstart,ystart] = plotContourMatrix(ax(laxes),Ctmp);
    HH = text(xstart(1),ystart(1),[num2str(zlevs(iii)/1000),' km.'],'FontSize',15);
    for jj = 1:length(HH)
        HH(jj).Rotation = 63; %.Interpreter = 'none';
    end
end

disp('contouring ISC');
data = importdata('~/igdata/wsm_020.0_regularization_Corrm_exponential_030_km_sm0_000_bayesian_sol_coupling.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);

laxes = laxes +1; %new axes to have custom color scale for the isc contours
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

contour(ax(laxes),xq,yq,vq,50:10:90,'linewidth',2);
cbar2 = colorbar(ax(laxes));
axis(ax(laxes),'equal');
hold(ax(laxes),'on');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
hotmap = flipud(hot(512));
colormap(ax(laxes),hotmap);
ax(laxes).Visible = 'off';
cbar2.Position = [colorbarRightShift,0.55,0.0083333333333333,0.35];
cbar2.Limits = [50 90];
cbar2.Label.String = 'inter-seismic coupling [%]';
disp('done contouring ISC');
hold(ax(laxes),'on');
axis(ax(laxes),'equal');

%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

magFact = 0.3;
SS = bubblechart(ax(laxes),eqlon,eqlat,exp(eqmag),eqdepth);
bubblesize(magFact*[exp(min(eqmag)) min(exp(max(eqmag)))]);

SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.85;
SS.MarkerEdgeAlpha = 0.85;
SS.LineWidth = 1;
blgd = bubblelegend('Magnitude','Location','northwest');
blgd.LimitLabels = {num2str(round(min(eqmag)*100)/100),num2str(round(max(eqmag)*100)/100)};
blgd.Style = 'vertical';
blgd.NumBubbles = 2;

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

colormap(ax(laxes),'turbo');
cbar3 = colorbar(ax(laxes));
cbar3.Position = [colorbarRightShift,0.12,0.0083333333333333,0.35];
clim(ax(laxes),[0 15]);
cbar3.Label.String = 'Depth [km]';

ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

plot(ax(laxes),NaN,NaN,'.');

cbar4 = colorbar(ax(laxes));
cbar4.Visible = 'off';
cbar4.Position = [colorbarRightShift,0.01,0.0083333333333333,0.8];

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

titleSplit = split(titleStr,",");
print('-djpeg',strcat("Figure_",titleSplit(1),"_2021"));
end
