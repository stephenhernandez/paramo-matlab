clear; close all; %clc;

[minLon,maxLon,minLat,maxLat] = deal(-80,-79.6,1.1,1.4);
colorbarRightShift = 0.93;
tStart = datetime(2010,01,01);
tEnd = datetime(2024,05,01);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');

% lat = ncread('~/igdata/dem/completemultigolfe2017.grd','y');
% lon = ncread('~/igdata/dem/completemultigolfe2017.grd','x');
% demData = ncread('~/igdata/dem/completemultigolfe2017.grd','/z');
% 
% lat = ncread('~/igdata/dem/EsmeraldasExpanded.grd','lat');
% lon = ncread('~/igdata/dem/EsmeraldasExpanded.grd','lon');
% demData = ncread('~/igdata/dem/EsmeraldasExpanded.grd','/z');

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
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
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

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

%% scale ruler
refShifter = (maxLat-minLat)*(0.01/0.3);
scaleRefLat0 = minLat + refShifter; 
scaleRefLon0 = minLon + refShifter;
refEllipse = referenceEllipsoid('wgs84');

[latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);
[latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);
[latout20,lonout20] = reckon(scaleRefLat0,scaleRefLon0,20*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
scaleRefLon2 = scaleRefLon0;
plot(ax(laxes),[scaleRefLon0 lonout20],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -refShifter/2;
textFontSize = 24;
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
magFact = 6;
minMag = 0;
maxDepth = 500;
[~,~,t,eqlat,eqlon,eqdepth,eqmag,ids] = readCat1(tStart,tEnd,minMag);
eI = eqdepth <= maxDepth & eqlon >= minLon & eqlon <= maxLon & eqlat >= minLat & eqlat <= maxLat;

t(~eI) = [];
eqlat(~eI) = [];
eqlon(~eI) = [];
eqdepth(~eI) = [];
eqmag(~eI) = [];
ids(~eI) = [];

[t,eI] = sort(t); %,'descend');
eqmag = eqmag(eI);
eqlat = eqlat(eI);
eqlon = eqlon(eI);
eqdepth = eqdepth(eI);
ids = ids(eI);

%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

SS = bubblechart(ax(laxes),eqlon,eqlat,magFact*exp(eqmag),eqdepth);
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
clim(ax(laxes),[0 10]);
cbar3.Label.String = 'Depth [km]';

StationMarkerSize = 20;

% contour(ax(laxes),lon,lat,demData',(-(100:100:4000))','k-','linewidth',0.2);

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

%%
cd ~/research/now/esmeraldas/hypodd/
!\rm -rf hypoDD.inp; \rm -rf ph2dt_esmeraldas.inp;
%!ln -s hypoDD_M3d_2.inp hypoDD.inp
%!ln -s hypoDD_MSERGIO.inp hypoDD.inp
%!ln -s hypoDD_SERGIO.inp hypoDD.inp
%!ln -s hypoDD_CENTRAL.inp hypoDD.inp
%!ln -s hypoDD_MACQUET.2.inp hypoDD.inp
%!ln -s hypoDD_MFONT.inp hypoDD.inp
!ln -s hypoDD_SOTO.inp hypoDD.inp
!ln -s ph2dt_sandro.inp ph2dt_esmeraldas.inp

%
regionName = "esmeraldas";
boundaryBox = [minLon; maxLon; minLat; maxLat];
depthCorrection = 0; %positive number for elevation above sea-level
diasFlag = false;
minMag = 2;
maxDepth = 50; 
depthRange = [15 0];
[E,t,eqlat,eqlon,eqdepth,eqmag,id,tOrig,eqlatOrig,eqlonOrig,eqdepthOrig,eqmagOrig,idOrig] = ...
    genHypoDDFiles(tStart,tEnd,regionName,boundaryBox,depthCorrection,diasFlag,minMag,maxDepth,depthRange);
