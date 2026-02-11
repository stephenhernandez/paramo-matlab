clear; close all; %clc;

[minLon,maxLon,minLat,maxLat] = deal(-78.52,-78.34,-0.75,-0.62);
%[minLon,maxLon,minLat,maxLat] = deal(-78.26,-77.32,-0.55,0.46);
%[minLon,maxLon,minLat,maxLat] = deal(-78.7,-77.3,-0.6,0.6);
reve_lat = -0.0809;
reve_lon = -77.65835;

load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
%[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
lat = ncread('~/igdata/dem/Cotopaxi.grd','lat');
lon = ncread('~/igdata/dem/Cotopaxi.grd','lon');
demData = ncread('~/igdata/dem/Cotopaxi.grd','/z');
demData = double(demData);
kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];

%%
contrast = 1;
ilumAz = -90;
wMark = 2; %"brightness"... higher number means more washed out tones
elevAngle = 70;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);


I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot',...
    'Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1]); %flipud(summer(512)));
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [0.77,0.55,0.0083333333333333,0.35];

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

shortRuler = true;
if shortRuler
    refShifter = (maxLat-minLat)*(0.0095/0.3);
    scaleRefLat0 = minLat + refShifter;
    scaleRefLon0 = minLon + refShifter;
    refEllipse = referenceEllipsoid('wgs84');

    [latout1,lonout1] = reckon(scaleRefLat0,scaleRefLon0,1*1e3,90,refEllipse);
    [latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);
    [latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);

    scaleRefLat2 = scaleRefLat0 + refShifter/2;
    scaleRefLon2 = scaleRefLon0;
    plot(ax(laxes),[scaleRefLon0 lonout5],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

    textShifter = -refShifter/2;
    textFontSize = 24;
    plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

    plot(ax(laxes),[lonout1 lonout1],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+lonout1,-textShifter+scaleRefLat2,"1","FontSize",textFontSize);

    plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+lonout5,-textShifter+scaleRefLat2,"5 km","FontSize",textFontSize);

else
    % intermediate scaleruler
    refShifter = (maxLat-minLat)*(0.008/0.3);
    scaleRefLat0 = minLat + refShifter;
    scaleRefLon0 = minLon + refShifter;
    refEllipse = referenceEllipsoid('wgs84');

    [latout10,lonout10] = reckon(scaleRefLat0,scaleRefLon0,10*1e3,90,refEllipse);
    [latout20,lonout20] = reckon(scaleRefLat0,scaleRefLon0,20*1e3,90,refEllipse);
    [latout40,lonout40] = reckon(scaleRefLat0,scaleRefLon0,40*1e3,90,refEllipse);

    scaleRefLat2 = scaleRefLat0 + refShifter/2;
    scaleRefLon2 = scaleRefLon0;
    plot(ax(laxes),[scaleRefLon0 lonout40],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

    textShifter = -refShifter/2;
    textFontSize = 18;
    plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

    plot(ax(laxes),[lonout10 lonout10],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+lonout10,-textShifter+scaleRefLat2,"10","FontSize",textFontSize);

    plot(ax(laxes),[lonout20 lonout20],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+lonout20,-textShifter+scaleRefLat2,"20","FontSize",textFontSize);

    plot(ax(laxes),[lonout40 lonout40],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
    text(textShifter+lonout40,-textShifter+scaleRefLat2,"40 km","FontSize",textFontSize);
end

labelFontSize = 28;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);

%% permanent stations
StationMarkerSize = 22;
[stla,stlo] = metaDataFromStationList(kstnms);
h2 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);
text(stlo,stla,kstnms,'FontSize',20);


ax(laxes).XAxisLocation = 'top';
contour(ax(laxes),lon,lat,demData',(2000:250:5900)','k-','linewidth',0.2);
contour(ax(laxes),lon,lat,demData',(2000:500:5900)','k-','linewidth',1);
ax(laxes).XAxis.FontSize = labelFontSize;
ax(laxes).YAxis.FontSize = labelFontSize;
ax(laxes).YTickLabelRotation = -20;
ax(laxes).XTickLabelRotation = -20;
ax(laxes).LineWidth = 1;
