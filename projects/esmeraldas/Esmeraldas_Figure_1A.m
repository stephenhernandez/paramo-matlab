% esmeraldas maps
clear; close all; clc;

[minLon,maxLon,minLat,maxLat] = deal(-82,-77.5,-2,1.75);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
%[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
lat = ncread('~/igdata/dem/esmeraldas_1.grd','lat');
lon = ncread('~/igdata/dem/esmeraldas_1.grd','lon');
demData = ncread('~/igdata/dem/esmeraldas_1.grd','/z');
demData = double(demData);

contrast = 1/4;
ilumAz = 90;
wMark = 2.5;
elevAngle = 50;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1]);
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [0.77,0.55,0.0083333333333333,0.35];

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

% scaleruler
refShifter = (maxLat-minLat)*(0.005/0.3);
scaleRefLat0 = minLat + refShifter;
scaleRefLon0 = minLon + refShifter;
refEllipse = referenceEllipsoid('wgs84');

[latout25,lonout25] = reckon(scaleRefLat0,scaleRefLon0,25*1e3,90,refEllipse);
[latout50,lonout50] = reckon(scaleRefLat0,scaleRefLon0,50*1e3,90,refEllipse);
[latout100,lonout100] = reckon(scaleRefLat0,scaleRefLon0,100*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
scaleRefLon2 = scaleRefLon0;
plot(ax(laxes),[scaleRefLon0 lonout100],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -refShifter/2;
textFontSize = 14;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout25 lonout25],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout25,-textShifter+scaleRefLat2,"25","FontSize",textFontSize);

plot(ax(laxes),[lonout50 lonout50],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout50,-textShifter+scaleRefLat2,"50","FontSize",textFontSize);

plot(ax(laxes),[lonout100 lonout100],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout100,-textShifter+scaleRefLat2,"100 km","FontSize",textFontSize);

labelFontSize = 22;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);

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
cbar2.Position = [0.82,0.55,0.0083333333333333,0.35];
cbar2.Limits = [50 90];
cbar2.Label.String = 'inter-seismic coupling [%]';
disp('done contouring ISC');

%
laxes = laxes +1; %new axes to have custom color scale for the isc contours
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

disp('contouring coseismic slip');
data = importdata('~/igdata/slip_53_s_m.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
contour(ax(laxes),xq,yq,vq,0:1,'linewidth',4,'color',[0 0 0]);
cbar22 = colorbar(ax(laxes));
axis(ax(laxes),'equal');
hold(ax(laxes),'on');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
hotmap = flipud(hot(512));
colormap(ax(laxes),hotmap);
ax(laxes).Visible = 'off';
cbar22.Visible = 'off';
cbar22.Position = [0.82,0.55,0.0083333333333333,0.35];
disp('done contouring coseismic slip');

%
%cd ~/Dropbox/Esmeraldas2022/shapefiles/rupture-slowslip/;
cd ~/igdata/shape_files/rupture-slowslip
T = readgeotable('Swenson_Beck_maxmoment.shp');
T = geotable2table(T,["Lat","Lon"]);

for i = 1:size(T,1)
    plot(ax(laxes),cell2mat(T.Lon(i))',cell2mat(T.Lat(i))','linewidth',4,'color',[0.47 0.67 0.19]);
end

T = readgeotable('Bilek_2021_1906rup_geog.shp');
T = geotable2table(T,["Lat","Lon"]);

%% faults
% cd ~/research/now/esmeraldas/shapefiles/faults/;
% T = readgeotable('faults_polylines.shp');
% H = geoshow(T);
% lineObjs = findobj(H,'Type','Line');
% for i = 1:length(lineObjs)
%     linobj_ = lineObjs(i);
%     linobj_.Color(4) = 0.25;
% end

regionName = "esmeraldas";
[t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD('~/igdata/esmeraldas2022.hypodd.reloc',true);
[eqmag,sI] = sort(eqmag,'descend');
t = t(sI);
eqlat = eqlat(sI);
eqlon = eqlon(sI);
eqdepth = eqdepth(sI);
id = id(sI);

refVMdepth = 0;
eqdepth = eqdepth-refVMdepth;
magFact = 10;

%%
StationMarkerSize = 15;

%% temporary network
% [stlaH,stloH] = metaDataFromStationList(["HB10";"HB12";"HB15";"HB17";"HB18";"HB19";"HB56"]);
% h1 = plot(ax(laxes),stloH,stlaH,'v','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);

allStationKstnms = ["AAT1","AATC","AES1","AES2","ALIT","ALOR","ANGU",...
    "ANTG","ANTI","ANTM","ANTS","AOTA","APED","ARA2","ASAM","ASDO","ATON",...
    "BBAC","BMAS","BMOR","BNAS","BONI","BPAT","BREF","BTAM","BTER","BVC2",...
    "CAB1","CAYR","CERN","CHL1","CHL2","CHMA","CHSH","CIVI","COTA","CUIC",...
    "CUSE","CUSW","ECEN","ESM1","FLF1","GGPC","GGPT","ILLI","IMBA",...
    "JAMA","JIPI","JSCH","JUA2","LITA","LNGL","MONB","NAS2","OTAV","PAC1",...
    "PAS1","PECV","PINO","PIS1","PITA","PORT","PTGL","PUEM","PULU","REVS",...
    "SLOR","SNLR","SRAM","SUCR","TAMB","TAMH","TING","TOMA","URCU",...
    "VC1","VCES","YAHU","YANA"];

%% permanent network
[stla,stlo] = metaDataFromStationList(allStationKstnms);
h2 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);
[legh,legIcons] = legend(h2,'RENSIG','Location','NorthWest');

%% seiscomp catalog
minMag = 5;
maxDepth = 50;
[~,~,tRegional,eqlatRegional,eqlonRegional,eqdepthRegional,eqmagRegional] = readCat1(datetime(2013,01,01),datetime(2022,03,27),minMag);
eI = eqdepthRegional <= maxDepth;
tRegional(~eI) = [];
eqlatRegional(~eI) = [];
eqlonRegional(~eI) = [];
eqdepthRegional(~eI) = [];
eqmagRegional(~eI) = [];
plot(ax(laxes),eqlonRegional,eqlatRegional,'k.');

%% cross section profile
% C1lon = -80.1509364;
% C1lat = 1.0545203;
% C2lon = -79.3009142;
% C2lat = 0.6332050;
% plot(ax(laxes),[C1lon C2lon],[C1lat C2lat],'k','linewidth',2);

%%
GPS_STLA(1) = 0.78153;
GPS_STLO(1) = -80.03036;
GPS_STNM(1) = "PTGL";

GPS_STLA(2,1) = 0.86165;
GPS_STLO(2,1) = -79.85933;
GPS_STNM(2,1) = "TRAT";

GPS_STLA(3,1) = 0.97878;
GPS_STLO(3,1) = -79.62541;
GPS_STNM(3,1) = "TRTH";

% h3 = scatter(ax(laxes),GPS_STLO,GPS_STLA,StationMarkerSize.^2,"magenta",'s','filled'); %%,'markeredgecolor','k','markerfacecolor','m','MarkerSize',StationMarkerSize);
% h3.MarkerEdgeColor = 'k';
% h3.MarkerFaceAlpha = 0.85;
% [legh,legIcons] = legend([h1;h2;h3;H],'Seismic (Temp.)','Seismic (Perm.)','GPS','Faults','Location','NorthWest');
% legIcons = findobj(legIcons, '-not', 'Marker', 'none');
% legIcons = findobj(legIcons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
% legIcons(4).MarkerSize = StationMarkerSize;

linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

