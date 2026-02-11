clear; close all; %clc;

shortRuler = true;
if shortRuler
    [minLon,maxLon,minLat,maxLat] = deal(-77.73,-77.58,-0.16,-0.01);
else
    [minLon,maxLon,minLat,maxLat] = deal(-78.95,-77.05,-0.6,1);
    %[minLon,maxLon,minLat,maxLat] = deal(-78.26,-77.32,-0.55,0.46);
end

reve_lat = -0.0809;
reve_lon = -77.65835;

load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
% lat = ncread('~/igdata/dem/Reventador2.grd','lat');
% lon = ncread('~/igdata/dem/Reventador2.grd','lon');
% demData = ncread('~/igdata/dem/Reventador2.grd','/z');
demData = double(demData);

%%
contrast = 0.75; %1/3;
ilumAz = -90; %-90;
wMark = 1.7; %2; %"brightness"... higher number means more washed out tones
elevAngle = 45; %70;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

dem_poster_palette = [4900 255 255 255 255;...
    4000    255 255 255 255;...
    2800    110 110 110 255;...
    1700    158   0   0 255;...
    1200    161  67   0 255;...
    500     232 215 125 255;...
    50       16 122  47 255;...
    0         0  97  71 255];

elevs = dem_poster_palette(:,1);
rgb = dem_poster_palette(:,2:end-1);
landcolormap = interp1(elevs,rgb,(1000:4000)')./255;

I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot',...
    'Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'LandColor',demcmap([500 4000])); %flipud(summer(512)));
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [0.77,0.55,0.0083333333333333,0.35];

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

%geoshow(ax(laxes),'~/igdata/ZonaUrbana/ZonaUrbana.shp');

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;


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

if shortRuler
    allBBStationKstnms = ["REVN";"REVS"];
    [stla,stlo] = metaDataFromStationList(allBBStationKstnms);
    h2 = plot(ax(laxes),stlo,stla,'v','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);

    [stla,stlo] = metaDataFromStationList("REVN");
    h6 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);
    h3 = plot(ax(laxes),reve_lon,reve_lat,'o','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);

    allSPStationKstnms = "LAV4";
    [stla,stlo] = metaDataFromStationList(allSPStationKstnms);
    h4 = plot(ax(laxes),stlo,stla,'<','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);
    allSPStationKstnms = "LAV4";
    [stla,stlo] = metaDataFromStationList(allSPStationKstnms);
    h5 = plot(ax(laxes),stlo,stla,'>','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);

    refEllipse = referenceEllipsoid('wgs84');
    [d_,az] = distance(reve_lat,reve_lon,stla,stlo,refEllipse);
    [~,baz] = distance(stla,stlo,reve_lat,reve_lon,refEllipse);
    d_ = d_*1e-3;

    [legh,legIcons] = legend([h2; h4; h5; h6],'Broadband Seismometer',....
        'Short Period Seismometer','Ashmeter','DOAS','Location','NorthEastOutside');
    legIcons = findobj(legIcons, '-not', 'Marker', 'none');
    legIcons = findobj(legIcons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
    legIcons(1).MarkerSize = StationMarkerSize-4;

    %texti = 1;
    %text(ax(laxes),refShifter+stlo(texti),refShifter+stla(texti),"Ashmeter",'FontSize',labelFontSize);



    texti = 1;
    text(ax(laxes),refShifter+reve_lon(texti),refShifter+reve_lat(texti),"Vent",'FontSize',labelFontSize);

    allStationKstnms = ["REVN";"REVS";"LAV4"];
    [stla,stlo] = metaDataFromStationList(allStationKstnms);
    text(ax(laxes),0.0035+stlo,0.001+stla,allStationKstnms,'FontSize',labelFontSize);
    ax(laxes).XAxisLocation = 'top';
    ax(laxes).YAxisLocation = 'left';
else

    allStationKstnms = ["CASC";"BONI";"ANTS";"ANTG"];
    [stla,stlo] = metaDataFromStationList(allStationKstnms);
    h2 = plot(ax(laxes),stlo,stla,'v','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize,'linewidth',2);
    %-77.72,-77.58,-0.16,0);
    plot(ax(laxes),[-77.73;-77.58;-77.58;-77.73;-77.73],...
        [-0.16;-0.16;-0.01;-0.01;-0.16],'k-','linewidth',2);
    quito_lat = -0.22;
    quito_lon = -78.523;
    plot(ax(laxes),quito_lon,quito_lat,'d','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);
    texti = 1;
    text(ax(laxes),0.01+quito_lon(texti),-0.02+quito_lat(texti),"Quito",'FontSize',labelFontSize);
    texti = 1;
    text(ax(laxes),-0.1+stlo(texti),0.04+stla(texti),allStationKstnms(texti),'FontSize',labelFontSize); %CASC

    texti = 2;
    text(ax(laxes),-0.1+stlo(texti),0.04+stla(texti),allStationKstnms(texti),'FontSize',labelFontSize); %BONI

    texti = 3;
    text(ax(laxes),0.03+stlo(texti),0.01+stla(texti),allStationKstnms(texti),'FontSize',labelFontSize);%ANTS

    texti = 4;
    text(ax(laxes),-0.1+stlo(texti),0.04+stla(texti),allStationKstnms(texti),'FontSize',labelFontSize);%ANTG
    ax(laxes).XAxisLocation = 'top';
end

% contour(ax(laxes),lon,lat,demData',(1000:500:4000)','k-','linewidth',0.2);
% contour(ax(laxes),lon,lat,demData',(1000:1000:4000)','k-','linewidth',1);

ax(laxes).XAxis.FontSize = labelFontSize;
ax(laxes).YAxis.FontSize = labelFontSize;
ax(laxes).YTickLabelRotation = 30;
ax(laxes).XTickLabelRotation = 30;
ax(laxes).LineWidth = 1;

% cd ~/FigCatchAll/
% print -depsc reve_large_v4

cd ~/FigCatchAll/reventador/
print -depsc reve_short_v4
