% esmeraldas maps
clear; close all; %clc;

%[minLon,maxLon,minLat,maxLat] = deal(-80.2,-79.2,0.6,1.1);
[minLon,maxLon,minLat,maxLat] = deal(-80.15,-79.35,0.6,1.1);
%[minLon,maxLon,minLat,maxLat] = deal(-79.9,-79.65,0.7,1.05);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);

contrast = 1/2;
ilumAz = 90;
wMark = 2;
elevAngle = 50;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

[I,dem_cmap] = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'NoDecim','legend');
imagesc(ax(laxes),I.x,I.y,I.rgb);

cbar1 = colorbar(ax(laxes));
colormap(cbar1,dem_cmap)
clim([min(min(demData)) max(max(demData))]);
cbar1.Visible = 'off';

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

%% scaleruler
scaleRefLat0 = minLat + (maxLat-minLat)*(0.005/0.3); %0.755;
scaleRefLon0 = minLon + (maxLat-minLat)*(0.005/0.3); %-79.895;
refEllipse = referenceEllipsoid('wgs84');
% [latout1,lonout1] = reckon(scaleRefLat0,scaleRefLon0,1*1e3,90,refEllipse);
% [latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);
% [latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);
% 
% scaleRefLat2 = scaleRefLat0 + 0.0035;
% scaleRefLon2 = scaleRefLon0;
% plot(ax(laxes),[scaleRefLon0 lonout5],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);
% 
% textShifter = -0.003;
% textFontSize = 12;
% plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
% text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);
% 
% plot(ax(laxes),[lonout1 lonout1],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
% text(textShifter+lonout1,-textShifter+scaleRefLat2,"1","FontSize",textFontSize);
% 
% plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
% text(textShifter+lonout5,-textShifter+scaleRefLat2,"5 km","FontSize",textFontSize);
% 
% %plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-');
% %text(textShifter+lonout10,-textShifter+scaleRefLat2,"10 km","FontSize",textFontSize);

[latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);
[latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);
[latout20,lonout20] = reckon(scaleRefLat0,scaleRefLon0,20*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + 0.004;
scaleRefLon2 = scaleRefLon0;
plot(ax(laxes),[scaleRefLon0 lonout20],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -0.004;
textFontSize = 12;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout5,-textShifter+scaleRefLat2,"5","FontSize",textFontSize);

plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout10,-textShifter+scaleRefLat2,"10","FontSize",textFontSize);

plot(ax(laxes),[lonout20 lonout20],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout20,-textShifter+scaleRefLat2,"20 km","FontSize",textFontSize);

%plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-');
%text(textShifter+lonout10,-textShifter+scaleRefLat2,"10 km","FontSize",textFontSize);

%%
labelFontSize = 22;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);

pointColormap = 'parula';
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
laxes = laxes + 1;
ax(laxes) = axes;

SS = bubblechart(ax(laxes),eqlon,eqlat,magFact*exp(eqmag),eqdepth); %,'filled');
SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.8;
SS.MarkerEdgeAlpha = 0.8;
SS.LineWidth = 1;
% blgd = bubblelegend('Magnitude','Location','northwest');
% blgd.LimitLabels = {'0.8','6'};
% blgd.Style = 'vertical';

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

cbar2 = colorbar(ax(laxes));
colormap(cbar2,pointColormap);
clim([10 20]);
cbar2.Label.String = 'Depth [km]';
%cbar2.Visible = 'off';

StationMarkerSize = 18;
[stlaH,stloH] = metaDataFromStationList(["HB10";"HB12";"HB15";"HB17";"HB18";"HB19";"HB56"]);
h1 = plot(ax(laxes),stloH,stlaH,'s','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);

allStationKstnms = ["AAT1","AATC","AES1","AES2","ALIT","ALOR","ANGU",...
    "ANTG","ANTI","ANTM","ANTS","AOTA","APED","ARA2","ASAM","ASDO","ATON",...
    "BBAC","BMAS","BMOR","BNAS","BONI","BPAT","BREF","BTAM","BTER","BVC2",...
    "CAB1","CAYR","CERN","CHL1","CHL2","CHMA","CHSH","CIVI","COTA","CUIC",...
    "CUSE","CUSW","ECEN","ESM1","FLF1","GGPC","GGPT","ILLI","IMBA",...
    "JAMA","JIPI","JSCH","JUA2","LITA","LNGL","MONB","NAS2","OTAV","PAC1",...
    "PAS1","PECV","PINO","PIS1","PITA","PORT","PTGL","PUEM","PULU","REVS",...
    "SLOR","SNLR","SRAM","SUCR","TAMB","TAMH","TING","TOMA","URCU",...
    "VC1","VCES","YAHU","YANA"];

[stla,stlo] = metaDataFromStationList(allStationKstnms);
h2 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);

[~,~,tRegional,eqlatRegional,eqlonRegional,eqdepthRegional,eqmagRegional] = readCat1(datetime(2013,01,01),datetime(2022,03,25),4);
plot(ax(laxes),eqlonRegional,eqlatRegional,'k.');

C1lon = -80.1509364;
C1lat = 1.0545203;
C2lon = -79.3009142; 
C2lat = 0.6332050;
%plot(ax(laxes),[C1lon C2lon],[C1lat C2lat],'k','linewidth',2);

GPS_STLA(1) = 0.78153;
GPS_STLO(1) = -80.03036;
GPS_STNM(1) = "PTGL";

GPS_STLA(2,1) = 0.86165;
GPS_STLO(2,1) = -79.85933;
GPS_STNM(2,1) = "TRAT";

GPS_STLA(3,1) = 0.97878;
GPS_STLO(3,1) = -79.62541;
GPS_STNM(3,1) = "TRTH";

h3 = plot(ax(laxes),GPS_STLO,GPS_STLA,'d','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);
legend([h1;h2;h3],'Seismic (Temp.)','Seismic (Perm.)','GPS','Location','NorthWest')

linkaxes(ax);
ax(laxes).Visible = 'off';
axis(ax(laxes),[minLon maxLon minLat maxLat]);

