% esmeraldas maps
clear; close all; %clc;
E = readSCBulletin("igepn2025ceos");
%E = readSCBulletin("igepn2025dtkg");

%%
[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1(datetime(2025,01,31));

width = 0.2;
height = width;
[minLon,maxLon,minLat,maxLat] = deal(E.lon-width,E.lon+width,E.lat-height,E.lat+height);
lI = eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon;

table(evDescription(lI),evStatus(lI),t(lI),eqlat(lI),eqlon(lI),eqdepth(lI),...
    eqmag(lI),ids(lI),stderr(lI),azgap(lI),phaseTot(lI),nMLv(lI),timerr(lI),...
    eqlaterr(lI),eqlonerr(lI),eqdeptherr(lI),eqmagerr(lI),magType(lI),meanMethod(lI),...
    locMethod(lI),earthModel(lI),creationTime(lI),agencyID(lI),evMode(lI),...
    scbullmag(lI),authorID(lI),evType(lI))

%%
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
% lat = ncread('~/igdata/dem/esmeraldas_2.grd','lat');
% lon = ncread('~/igdata/dem/esmeraldas_2.grd','lon');
% demData = ncread('~/igdata/dem/esmeraldas_2.grd','/z');
demData = double(demData);

contrast = 3/4;
ilumAz = -45;
wMark = 1.75;
elevAngle = 45;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'nodecim');
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [0.86,0.55,0.0083333333333333,0.35];

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

[latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,1*1e3,90,refEllipse);
[latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,2*1e3,90,refEllipse);
[latout20,lonout20] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
scaleRefLon2 = scaleRefLon0;
plot(ax(laxes),[scaleRefLon0 lonout20],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -refShifter/2;
textFontSize = 14;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout5,-textShifter+scaleRefLat2,"1","FontSize",textFontSize);

plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout10,-textShifter+scaleRefLat2,"2","FontSize",textFontSize);

plot(ax(laxes),[lonout20 lonout20],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout20,-textShifter+scaleRefLat2,"5 km","FontSize",textFontSize);

%%
labelFontSize = 22;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);

pointColormap = 'turbo';
relocFlag = false;
if relocFlag
    [t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD('~/igdata/napo2025.hypodd.reloc',true);
    [eqmag,sI] = sort(eqmag,'descend');
    t = t(sI);
    eqlat = eqlat(sI);
    eqlon = eqlon(sI);
    eqdepth = eqdepth(sI);
    id = id(sI);
else
    eqmag = eqmag(lI);
    t = t(lI);
    eqlat = eqlat(lI);
    eqlon = eqlon(lI);
    eqdepth = eqdepth(lI);
    id = ids(lI);
    [eqmag,sI] = sort(eqmag,'descend');
    t = t(sI);
    eqlat = eqlat(sI);
    eqlon = eqlon(sI);
    eqdepth = eqdepth(sI);
    id = id(sI);
end

refVMdepth = 0;
eqdepth = eqdepth-refVMdepth;
magFact = 10;

%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

SS = bubblechart(ax(laxes),eqlon,eqlat,magFact*exp(eqmag),eqdepth); %,'filled');
SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.9;
SS.MarkerEdgeAlpha = 0.9;
SS.LineWidth = 1;
blgd = bubblelegend('Magnitude','Location','northwest');
blgd.LimitLabels = {num2str(round(min(eqmag)*10)/10),num2str(round(max(eqmag)*10)/10)};
blgd.Style = 'vertical';
blgd.NumBubbles = 2;

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

colormap(ax(laxes),'parula');
cbar3 = colorbar(ax(laxes));
cbar3.Position = [0.76,0.108,0.0083333333333333,0.82];
clim(ax(laxes),[4 20]);
cbar3.Label.String = 'Depth [km]';

%%
StationMarkerSize = 20;

allStationKstnms = ["AAT1","AATC","AES1","AES2","ALIT","ALOR","ANGU",...
    "ANTG","ANTI","ANTM","ANTS","AOTA","APED","ARA2","ASAM","ASDO","ATON",...
    "BBAC","BMAS","BMOR","BNAS","BONI","BPAT","BREF","BTAM","BTER","BVC2",...
    "CAB1","CAYR","CERN","CHL1","CHL2","CHMA","CHSH","CIVI","COTA","CUIC",...
    "CUSE","CUSW","ECEN","ESM1","FLF1","GGPC","GGPT","ILLI","IMBA",...
    "JAMA","JIPI","JSCH","JUA2","LITA","LNGL","MONB","NAS2","OTAV","PAC1",...
    "PAS1","PECV","PINO","PIS1","PITA","PORT","PTGL","PUEM","PULU","REVS",...
    "SLOR","SNLR","SRAM","SUCR","TAMB","TAMH","TING","TOMA","URCU",...
    "VC1","VCES","YAHU","YANA"];

%perm stations
% [stla,stlo] = metaDataFromStationList(allStationKstnms);
% h2 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);

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

%[legh,legIcons] = legend([h1;h2],'Seismic (Temp.)','Seismic (Perm.)','Location','NorthWest');
% legIcons = findobj(legIcons, '-not', 'Marker', 'none');
% legIcons = findobj(legIcons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
% legIcons(4).MarkerSize = StationMarkerSize;

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
cbar4.Position = [0.9,0.01,0.0083333333333333,0.8];

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);


%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);
hold(ax(laxes),'on');
axis(ax(laxes),'equal');

geoshow(ax(laxes),'~/igdata/fallasactualizado/fallas2008completas.shp','Color',[0.5 0.5 0.5],'linewidth',1);
geoshow(ax(laxes),'~/igdata/ZonaUrbana/ZonaUrbana.shp');

ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

% %igepn2022fzqp 202203270428 0.9400 -79.8100 25.0 6.0   0.961 -79.845 16.1 5.80   25 22 119   174 70 78   0.29 11 9 9 11    244 A
% mt=sdr2mt(25,22,119); mlv = 6; focalmech(mt,-79.89,1.07,0.015*sqrt(mlv));
% 
% %igepn2022fztk 202203270553 0.8027 -79.7359 17.5 5.1   0.853 -79.785 17.4 4.97   23 21 121   170 72 78   0.34 17 13 12 18  202 A
% mt=sdr2mt(23,21,121); mlv = 4.97; focalmech(mt,-79.85,1.07,0.015*sqrt(mlv));
% 
% %igepn2022fzwc 202203270714 0.8401 -79.7534 16.8 3.6   0.888 -79.792 17.4 3.68   34 28 129   171 68 71   0.11 5 2 4 5      148 A
% mt=sdr2mt(34,28,129); mlv = 3.68; focalmech(mt,-79.81,1.07,0.015*sqrt(mlv));
% 
% %igepn2020gavx 202203272017 0.8454 -79.7476 15.8 5.3   0.895 -79.797 17.3 5.32   30 23 124   173 71 76   0.22 19 14 13 19  204 A
% mt=sdr2mt(30,23,124); mlv = 5.32; focalmech(mt,-79.77,1.07,0.015*sqrt(mlv));
% 
% %igepn2022gbhe 202203280159 0.8948 -79.7720 14.0 3.4   0.897 -79.758 14.7 3.42   185 48 94   359 42 85   0.16 9 3 3 9      176 A
% mt=sdr2mt(185,48,94); mlv = 3.42; focalmech(mt,-79.73,1.07,0.015*sqrt(mlv));
% 
% %igepn2022gcdu 202203281324 0.7499 -79.7946 14.1 3.3   0.780 -79.795 15.1 3.46   31 29 115   182 63 76   0.18 14 2 4 14    103 A
% mt=sdr2mt(31,29,115); mlv = 3.46; focalmech(mt,-79.69,1.07,0.015*sqrt(mlv));
% 
% %igepn2022gmav 202204022316 0.9338 -79.8036 14.6 3.3   0.976 -79.782 19.3 3.54   357 34 138  123 68 63   0.19 7 1 1 8      193 B
% mt=sdr2mt(357,34,138); mlv = 3.54; focalmech(mt,-79.65,1.07,0.015*sqrt(mlv));
% 
% %igepn2022gmhu 202204030246 0.9408 -79.7924 13.4 3.4   0.919 -79.773 18.0 3.40   313 38 70   157 54 105  0.12 7 3 4 8      202 A
% mt=sdr2mt(313,38,70); mlv = 3.4; focalmech(mt,-79.61,1.07,0.015*sqrt(mlv));
% 
% %igepn2022ixab 202205071027 0.8080 -79.7189 11.5 3.1   0.758 -79.763 24.2 3.11   223 74 86   57 16 103   0.18 4 1 2 4      171 A
% mt=sdr2mt(223,74,86); mlv = 3.11; focalmech(mt,-79.57,1.07,0.015*sqrt(mlv));
