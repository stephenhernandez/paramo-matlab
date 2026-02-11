clear; close all; clc;

[minLon,maxLon,minLat,maxLat] = deal(-91.7,-91.3,-0.55,-0.25);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
demData = double(demData);

%%
colorbarRightShift = 0.93;
demContrast = 1;
ilumAz = -90;
wMark = 1.8; %"brightness"... higher number means more washed out tones
elevAngle = 70;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',demContrast,'Azimuth',ilumAz,'noplot',...
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

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

refEllipse = referenceEllipsoid('wgs84');
shortRuler = true;
if shortRuler
    refShifter = (maxLat-minLat)*(0.0095/0.3);
    scaleRefLat0 = minLat + refShifter;
    scaleRefLon0 = minLon + refShifter;

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

allBBStationKstnms = ["FER1";"FER2"];
[stla,stlo] = metaDataFromStationList(allBBStationKstnms);
allBBStationKstnms = [allBBStationKstnms; "FER3";"CBDG";"CBHM"];
stla = [stla; -0.3508; -0.30626; -0.46960];
stlo = [stlo; -91.38539; -91.65126; -91.60975];
h2 = plot(ax(laxes),stlo,stla,'v','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);
text(ax(laxes),0.0035+stlo,0.001+stla,allBBStationKstnms,'FontSize',labelFontSize);

contour(ax(laxes),lon,lat,demData',(50:50:2000)','w-','linewidth',0.2);
contour(ax(laxes),lon,lat,demData',(50:100:2000)','k-','linewidth',1);

ax(laxes).XAxis.FontSize = labelFontSize;
ax(laxes).YAxis.FontSize = labelFontSize;
ax(laxes).YTickLabelRotation = 30;
ax(laxes).XTickLabelRotation = 30;
ax(laxes).LineWidth = 1;

readList = ["~/igdata/Fernandina/Fernandina2018.shp";...
    "~/igdata/Fernandina/Fernandina2017.shp";...
    "~/igdata/Fernandina/Fernandina2020_3mar"];
zones = [15;15;-15];
hold(ax(laxes),'on')
for k = 1:length(readList)
    fName = readList(k);
    S = shaperead(fName,'UseGeoCoords',true);
    if length(S)==1
        [slat,slon] = utm2ll(S.Lon,S.Lat,zones(k));
        slat = slat'; slon = slon';
        plot(ax(laxes),slon,slat,'r','linewidth',2);
    else
        for i = 1:length(S)
            [slat,slon] = utm2ll(S(i).Lon,S(i).Lat,zones(k));
            slat = slat'; slon = slon';
            badI = ~isfinite(slon) | ~isfinite(slat);
            if ~sum(badI)
                plot(ax(laxes),slon,slat,'r','linewidth',2);
                zoom on; grid on;
            else
                badI = find(badI);
                si = [1; badI(1:end-1) + 1];
                ei = badI - 1;
                for j = 1:length(si)
                    plot(ax(laxes),slon(si(j):ei(j)),slat(si(j):ei(j)),'r','linewidth',2);
                    zoom on; grid on;
                end
            end
        end
    end
end

%%
tStart = datetime(2013,01,01);
tEnd = datetime(2024,03,05);
magFact = 2;
minMag = 2;
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

% cd ~/FigCatchAll/
% print -depsc Fer_v1
