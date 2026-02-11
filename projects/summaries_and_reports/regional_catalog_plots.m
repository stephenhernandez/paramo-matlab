function [t,eqlat,eqlon,eqdepth,eqmag,eqids,...
    ax2_1,ax2_2,ax2_3,ax2_4] = regional_catalog_plots(varargin)
% Inputs:
% regionName = 'cotopaxi';
% tbeg = datenum(2016,01,01);
% tend = datenum(2016,04,01);
% typeLabel = 'VT';
% markerSymbol = '+';
% language = 'spanish'
% Nfilt = 11;
% magFact1 = 1.5;
% magFact2 = 1.;
% colorPalette = 'jet';
% maxDistCrat = 10;
% maxDepth = 20;
% nDays = 7;
% contourFlag = false; %to plot (or not) contours


%% parse inputs, deal variables
nVarargin = length(varargin);
functionDefaults = {'cotopaxi',datetime(2016,01,01),datetime(2017,01,01),'ALL',10,-30,70,'o',...
    'spanish',21,'default',-10,50,7,false};
optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[regionName,tStart,tEnd,typeLabel,magFact,minMag,maxMag,~,... %`markerSymbol' in tilde...
    language,Nfilt,colorPalette,minDepth,maxDepth,nDays,plotContour] = deal(optsToUse{:});
rotAngle = 20;

% beginning and end time
fprintf("Begin time: %s; End time: %s\n",tStart,tEnd);

if strcmp(language,'english')
    if strcmp(typeLabel,'ALL')
        figure1Title = ['All catalog event types between: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3Title = ['Cumulative number of events: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3SubTitle = ['Event rate per day: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    else
        figure1Title = ['Catalog ',typeLabel, 's : ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3Title = ['Cumulative number of catalog ',typeLabel,'s: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3SubTitle = ['Event rate per day: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    end
    figure4Title = ['Event Magnitudes: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure5Title = ['Cumulative Energy: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure5SubTitle = ['Energy Rate: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure6Title = ['Depths: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    magLabel = 'Magnitude';
    lonLabel = 'Longitude';
    latLabel = 'Latitude';
    depthLabel = 'Depth [km.]';
    energyLabel = 'Energy [Joules]';
    cumulativeLabel = 'Cumulative Number';
    legendLabel = [num2str(Nfilt),'pt. Med.'];
else %default is spanish
    if strcmp(typeLabel,'ALL')
        figure1Title = ['Eventos en Catalogo Entre: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3Title = ['Numero acumulativo de eventos: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3SubTitle = ['Tasa diaria de sismicidad: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    else
        figure1Title = ['Eventos en catalogo de tipo ',typeLabel, ': ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3Title = ['Numero acumulativo de eventos de tipo ',typeLabel,'s: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
        figure3SubTitle = ['Tasa diaria de sismicidad: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    end
    figure4Title = ['Magnitud de Eventos: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure5Title = ['Energia Acumulada: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure5SubTitle = ['Tasa de Energia: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    figure6Title = ['Profundidad vs Tiempo: ',datestr(tStart,1),' - ',datestr(tEnd,1)];
    magLabel = 'Magnitud';
    lonLabel = 'Longitud';
    latLabel = 'Latitud';
    depthLabel = 'Profundidad [km.]';
    energyLabel = 'Energia [Joules]';
    cumulativeLabel = 'Numero Acumulativo';
    legendLabel = [num2str(Nfilt),'pt. Med.'];
end

%% read sc3 catalog and filter
parentDir = fullfile('~/FigCatchAll/',lower(regionName));
shape_dir = fullfile("~","igdata","shape_files");
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');

if ~exist(parentDir,'dir')
    mkdir(parentDir);
end
cd(char(parentDir));

boundaryBox = getRegionSpatialDimensions(regionName);
minLon = boundaryBox(1);
maxLon = boundaryBox(2);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);

%% dem stuff
[lon,lat,demData] = cutDEM(boundaryBox,true);
demData = medfilt2(demData',[3,3]);
numElements = numel(demData);
if numElements < 1e6
    [Xq,Yq] = meshgrid(linspace(min(lon),max(lon),3*length(lon)+1),linspace(min(lat),max(lat),3*length(lat)+1));
    demData = interp2(lon,lat,demData,Xq,Yq,"makima");
    lon = Xq(1,:)';
    lat = Yq(:,1);
end
maxValue = max(demData);
[~,maxLonI] = max(maxValue);

maxValue = max(demData,[],2);
[~,maxLatI] = max(maxValue);

ewindex = maxLatI;
nsindex = maxLonI;
center_lat = lat(maxLatI);
center_lon = lon(maxLonI);

[t,eqlat,eqlon,eqmag,eqdepth,rms,type,azgap,eqids,~,~,deptherr,~,used_phases] = ...
    readSC5Catalog(regionName,tStart,tEnd,typeLabel);
t = datenum(t);

%% read topo
ewcut = demData(ewindex,:);
nscut = demData(:,nsindex);

eI = isnan(ewcut);
ewcut(eI) = 0;
eI = isnan(nscut);
nscut(eI) = 0;

% some parameteres for hillshade option
ilumAz = -60;
wMark = 2;
elevAngle = 60;
contrast = 7/4;

%% read coastline
clear d;
refEllipse = referenceEllipsoid('wgs84');
d = distance(center_lat,center_lon,eqlat,eqlon,refEllipse);
d = d/1000;

%% various filters
rmsI = rms < 50;
nPhasesThresh = 3;
maxDepthError = 10;
phaseI = used_phases >= nPhasesThresh;
fprintf('minimum number of phases: %d\n',nPhasesThresh);
fprintf('max error in depth: %d\n',maxDepthError);
gapI = azgap <= 360;
magI = eqmag >= minMag & eqmag < maxMag;
dI = d <= 400*1e3;
depthI = eqdepth > minDepth & eqdepth <= maxDepth;
deI = deptherr <= maxDepthError & deptherr >= -999;
bI = eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon;

%% apply filters
sI = rmsI & phaseI & gapI & magI & dI & deI & depthI & bI;
t = t(sI);
eqlat = eqlat(sI);
eqlon = eqlon(sI);
eqmag = eqmag(sI);
eqids = eqids(sI);
eqdepth = eqdepth(sI);
type = type(sI);
deptherr = deptherr(sI);

%%
maxdepth = ceil(max(max(demData/1000)))+1;
mindepth = floor(min(-eqdepth));
minMag = floor(min(eqmag));
maxMag = ceil(max(eqmag));
fprintf('MinDepth: %f, MaxDepth: %f\n',mindepth,maxdepth);
fprintf('MinMag: %f, MaxMag: %f\n',minMag,maxMag);

%%
vt = contains(type,"vt",'ignorecase',true); fprintf('Numero de VTs: %d\n',sum(vt));
expl = contains(type,"exp",'ignorecase',true); fprintf('Numero de EXPs: %d\n',sum(expl));
vlp = contains(type,"vlp",'ignorecase',true); fprintf('Numero de VLPs: %d\n',sum(vlp));
lp = contains(type,"lp",'ignorecase',true); fprintf('Numero de LPs: %d\n',sum(lp));
hb = contains(type,"hb",'ignorecase',true); fprintf('Numero de HBs: %d\n',sum(hb));
trem = contains(type,"tre",'ignorecase',true); fprintf('Numero de TREMs: %d\n',sum(trem));
% reg = contains(type,"tect",'ignorecase',true) | ...
%    contains(type,"unk",'ignorecase',true);
reg = true(size(type)) & ~(vt | expl | vlp | lp | hb | trem); fprintf('Numero de REGs: %d\n',sum(reg));

fprintf('Numero Total: %d\n',sum(vt));

%% FIGURE 01
allPresent = [lp vt expl vlp hb reg trem];
allSum = sum(double(allPresent),1); %column-wise sum (greater than zero when at least one event found)
fprintf('Numero Total: %d\n',sum(allSum));
allI = find(allSum); %find which typelabels have triggered
allPresent = allPresent(:,allI);
lall = length(allI);
legH = gobjects(lall,1);
allMarkers = {'d','s','p','^','x','o','+'};
markerStr = {'LP','VT','EXP','VLP','HB','Regional','TREM'};
allMarkers = allMarkers(allI);
markerStr = markerStr(allI);

fig1outerposition = [0 0 0.8 1];
fig(1) = figure('units','normalized','outerposition',fig1outerposition);
ax1 = axes;
ax1.Parent = fig(1);
hold(ax1,"on");

I = dem(lon,lat,demData,'Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'NoDecim');
imagesc(ax1,I.x,I.y,I.rgb);
hold(ax1,'on');
axis(ax1,'xy');
plot(ax1,lonEC,latEC,'k','linewidth',2);
plot(ax1,lon_noec,lat_noec,'-','linewidth',2,'Color',[0.5 0.5 0.5]);
axis(ax1,'equal');

%%
minElev = min(min(demData));
maxElev = max(max(demData));
elevRange = maxElev - minElev;
if elevRange > 1200
    zInc = 500;
else
    zInc = 200;
end

if plotContour
    [~,cc] = contour(ax1,lon,lat,demData,(minElev:zInc:maxElev)');
    cc.LineWidth = 0.04;
    cc.Color = 'k';
end

%%
for i = 1:lall
    figure(fig(1));
    tmpLabel = allPresent(:,i);
    sI = find(tmpLabel);
    hh(i) = scatter(ax1,eqlon(sI),eqlat(sI),magFact*exp(eqmag(sI)),datenum(t(sI)),'filled',...
        'MarkerEdgeColor','k','MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75,'Marker',allMarkers{i});
    axis(ax1,'equal');
end

%
legend(hh,markerStr,'Location','NorthEastOutside','AutoUpdate','off');
cb3 = colorbar(ax1);
clim(ax1,datenum([tStart tEnd]));
cTickOrig = cb3.Ticks;

if length(cTickOrig) >= 7
    newTick = cTickOrig(1:2:length(cTickOrig));
else
    newTick = cTickOrig;
end

cb3.Ticks = newTick;
cb3.TickLabels = datestr(cb3.Ticks);
colormap(colorPalette);
sgtitle(figure1Title);

S = shaperead(fullfile(shape_dir,"volcanes_2.shp"));
plot(ax1,[S.X],[S.Y],'k--','linewidth',1);

% if strcmp(regionName,'chiles')
%     load('~/research/now/chiles/PuntosTermas');
%     plot(ax1_1,ptlon,ptlat,'d','MarkerEdgeColor','k','MarkerFaceColor','r');
%     text(ax1_1,ptlon,ptlat,ptstnm);
%     geoshow(ax1_1,'~/igdata/ZonaUrbana/ZonaUrbana.shp');
% end

% if strcmp(regionName,'fernandina') || strcmp(regionName,'galapagos')
%     S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2018.shp','UseGeoCoords',true);
%     [slat,slon] = utm2ll(S.Lon,S.Lat,15);
%     plot(ax1,slon,slat,'r','linewidth',2); zoom on; grid on;
%     S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2017.shp','UseGeoCoords',true);
%     [slat,slon] = utm2ll(S.Lon,S.Lat,15);
%     %figure(1); hold on;
%     plot(ax1,slon,slat,'r','linewidth',2); zoom on; grid on;
%     axis equal;
%     S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2020_3mar','UseGeoCoords',true);
%     for i = 1:length(S)
%         [slat,slon] = utm2ll(S(i).Lon,S(i).Lat,-15);
%         badI = find(~isfinite(slon))';
%         if ~sum(badI)
%             fill(ax1,slon,slat,'r'); zoom on; grid on;
%         else
%             si = [1; badI(1:end-1) + 1]; ei = badI - 1;
%             for j = 1:length(si)
%                 fill(ax1,slon(si(j):ei(j)),slat(si(j):ei(j)),'r'); zoom on; grid on;
%             end
%         end
%     end
% end

geoshow(ax1,fullfile(shape_dir,"ZonaUrbana.shp"));
geoshow(ax1,fullfile(shape_dir,"fallas2008completas.shp"),'Color','r','linewidth',0.5);
xlabel(ax1,lonLabel);
ylabel(ax1,latLabel);
linkaxes(ax1,'xy');
axis(ax1,boundaryBox);
ax1.Box = "on";
zoom on;
pause(1);

%% FIGURE 02
set(0,'units','pixels'); %reset to pixels
SCREENSIZE = get(0,'screensize');
toolbarFraction = 73/80; %this is fixed... may be a different ratio for different monitors/screen resolutions, do not apply blindly
S_HEIGHTPX = SCREENSIZE(4);
WIDTHKM = 111.19*(maxLon - minLon);     %deg2km
HEIGHTKM = 111.19*(maxLat - minLat);    %deg2km
DATAASPECTRATIO = WIDTHKM/HEIGHTKM;
D_WIDTHPX = round(toolbarFraction*S_HEIGHTPX*DATAASPECTRATIO);
fig2outerposition = [1 1 D_WIDTHPX S_HEIGHTPX];
fig(2) = figure('units','pixels','outerposition',fig2outerposition);
nrows = 48;
ncols = round(DATAASPECTRATIO*nrows);
TL = tiledlayout(nrows,ncols, 'Padding', 'compact', 'TileSpacing', 'none');
ax2_1(1) = nexttile(1,[round(2*nrows/3) round(2*ncols/3)]);

I = dem(lon,lat,demData,'Contrast',1,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'NoDecim');
imagesc(ax2_1,I.x,I.y,I.rgb);
hold(ax2_1,'on');
axis(ax2_1,'xy');
plot(ax2_1,lonEC,latEC,'k','linewidth',2);
plot(ax2_1,lon_noec,lat_noec,'-','linewidth',2,'Color',[0.5 0.5 0.5]);
if plotContour
    [~,cc] = contour(ax2_1,lon,lat,demData,(minElev:zInc:maxElev)');
    cc.LineWidth = 0.1;
    cc.Color = 'k'; %[0.5 0.5 0.5];
end
grid(ax2_1,'on');
ax2_1.XAxisLocation = 'top';
axis(ax2_1,boundaryBox);
axis(ax2_1,'equal');
box(ax2_1,'on');
ax2_1.YLabel.String = latLabel; % 'latitude';
ax2_1.XLabel.String = lonLabel; %'longitude';

ax2_2 = nexttile([nrows-round(2*nrows/3) round(2*ncols/3)]);%subplot(3,3,[3 6]);
plot(ax2_2,lon,ewcut/1000,'k','linewidth',2); %hold on;
hold(ax2_2,'on');
grid(ax2_2,'on');
box(ax2_2,'on');
ax2_2.XLabel.String = lonLabel; %'longitude';
ax2_2.YLabel.String = depthLabel; % 'depth [km.]';

ax2_3 = nexttile([round(2*nrows/3) ncols-round(2*ncols/3)]); %subplot(3,3,[7 8]);
plot(ax2_3,nscut/1000,lat,'k','linewidth',2);
hold(ax2_3,'on');
ax2_3.XAxis.Direction = 'reverse';
ax2_3.YAxis.Visible = 'off';
ax2_3.XAxisLocation = 'top';
box(ax2_3,'on');
grid(ax2_3,'on');
ax2_3.XLabel.String = depthLabel; %'depth [km.]';

figure(fig(2));
ax2_1(2) = axes(TL);
ax2_1(2).Layout.Tile = 1;
ax2_1(2).Layout.TileSpan = [round(2*nrows/3) round(2*ncols/3)];
hold(ax2_1(2),'on');
clear hh2;
for i = 1:lall
    tmpLabel = allPresent(:,i);
    sI = find(tmpLabel);
    hh = scatter(ax2_1(2),eqlon(sI),eqlat(sI),magFact*exp(eqmag(sI)),t(sI),'filled',...
        'MarkerEdgeColor','k');
    hh.Marker = allMarkers{i};
    hh.MarkerFaceAlpha = 0.75;
    hh.MarkerEdgeAlpha = 0.75;

    hh2(i) = scatter(ax2_2,eqlon(sI),-eqdepth(sI),magFact*exp(eqmag(sI)),t(sI),'filled',...
        'MarkerEdgeColor','k');
    hh2(i).Marker = allMarkers{i};
    hh2(i).MarkerFaceAlpha = 0.75;
    hh2(i).MarkerEdgeAlpha = 0.75;

    hh3 = scatter(ax2_3,-eqdepth(sI),eqlat(sI),magFact*exp(eqmag(sI)),t(sI),'filled',...
        'MarkerEdgeColor','k');
    hh3.Marker = allMarkers{i};
    hh3.MarkerFaceAlpha = 0.75;
    hh3.MarkerEdgeAlpha = 0.75;
end

ax2_1(2).Visible = 'off';
axis(ax2_1,boundaryBox);
axis(ax2_1,'equal')
linkaxes(ax2_1,'xy');

ax2_4 = nexttile([nrows-round(2*nrows/3) ncols-round(2*ncols/3)]);
h224 = histogram(eqmag,'BinLimits',[minMag,maxMag]);
if h224.NumBins < 2
    histogram(eqmag,2);
end
ax2_4.YAxisLocation = 'right';
box(ax2_4,'on');
cbar = colorbar;
cbar.Visible = 'off';
ax2_4.XLim = [floor(minMag),ceil(maxMag)];
ax2_4.YLabel.Interpreter = 'tex';
ax2_4.YLabel.String = '\propto P';
ax2_4.XLabel.String = magLabel;

linkaxes([ax2_1 ax2_2],'x');
axis(ax2_1(1),'equal')
axis(ax2_1(1),'tight')
axis(ax2_1,boundaryBox);
linkaxes([ax2_1 ax2_3],'y');

colormap(ax2_1(2),colorPalette);
colormap(ax2_2,colorPalette);
cbar = colorbar(ax2_3);
cbar.TickLabels = datestr(cbar.Ticks);
cbar.Label.String = "UTC";
colormap(ax2_3,colorPalette);
legend(hh2,markerStr,'Location','SouthWest');

% S = shaperead(fullfile(shape_dir,"volcanes_2.shp"));
% hold(ax2_1(end),'on');
% plot(ax2_1(end),[S.X],[S.Y],'k--','linewidth',1);
% 
% geoshow(ax2_1(end),fullfile(shape_dir,"fallas2008completas.shp"),'Color','r','linewidth',0.5);
% if strcmp(regionName,'chiles')
%     load('~/research/now/chiles/PuntosTermas');
%     plot(ax2_1(end),ptlon,ptlat,'d','MarkerEdgeColor','k','MarkerFaceColor','r');
%     text(ax2_1(end),ptlon,ptlat,ptstnm);
%     geoshow(ax2_1(end),fullfile(shape_dir,"ZonaUrbana.shp"));
% end
% set(gca,'Box','on');

linkaxes([ax2_1 ax2_2],'x');
axis(ax2_1(1),'equal')
axis(ax2_1(1),'tight')
axis(ax2_1,boundaryBox);
zoom on;
pause(1);

%% FIGURE 03 (CUMULATIVE NUMBER AND DAILY RATE)
fig(3) = figure('units','normalized','outerposition',[0.1 0.05 0.55 0.9]);
Ntot = length(t);
Ncum = (1:Ntot)';
rate = t2r(dn2dt(t),days(nDays));
if Ntot <= 50
    plot(dn2dt(t),Ncum,'k');
    ax = gca;
    hold on;
    for i = 1:lall
        tmpLabel = allPresent(:,i);
        sI = find(tmpLabel);
        legH(i) = plot(dn2dt(t(sI)),Ncum(sI),allMarkers{i});
    end
    ax.XTickLabelRotation = rotAngle;
    ylabel(cumulativeLabel);
    title(figure3Title);
    legend(legH,markerStr,"Location","NorthWest");
    grid on;
else
    tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');
    clear ax3;
    ax3(1) = nexttile();
    hold(ax3,"on");
    plot(ax3,dn2dt(t),Ncum,'k');
    for i = 1:lall
        tmpLabel = allPresent(:,i);
        sI = find(tmpLabel);
        legH(i) = plot(dn2dt(t(sI)),Ncum(sI),allMarkers{i});
    end
    ax3.XTickLabelRotation = rotAngle;
    ylabel(cumulativeLabel);
    legend(legH,markerStr,"Location","NorthWest");
    grid on; zoom on;
    %ax3(1).XAxis.Visible = "off";
    ax3(1).Box = "on";

    ax3(2,1) = nexttile();
    plot(ax3(2,1),dn2dt(t),rate/nDays,'.');
    ax3(2,1).XTickLabelRotation = rotAngle;
    ylabel('Event rate per day');
    grid on; zoom on;
    ax3(2,1).YAxisLocation = "right";
    grid on;
    linkaxes(ax3,'x');
end
pause(1);

%% FIGURE 04 (MAGNITUDES VS. TIME)
fig(4) = figure('units','normalized','outerposition',[0.1 0.05 0.55 0.9]);
tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');
clear ax4;
ax4(1) = nexttile();
hold(ax4,"on");

for i = 1:lall
    tmpLabel = allPresent(:,i);
    sI = find(tmpLabel);
    legH(i) = scatter(ax4,dn2dt(t(sI)),eqmag(sI),magFact*exp(eqmag(sI)),-eqdepth(sI),allMarkers{i},...
        'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7);
end
%ax4.XAxis.Visible = "off";
ylabel(magLabel);
grid on; zoom on;
cbar = colorbar;
cbar.Label.String = depthLabel;
legend(legH,markerStr,"Location","NorthWest");

ax4(2,1) = nexttile();
h = errorbar(ax4(2),dn2dt(t),-eqdepth,deptherr,'k.');
h.Marker = 'none';
hold(ax4(2),"on");
for i = 1:lall
    tmpLabel = allPresent(:,i);
    sI = find(tmpLabel);
    colormap(colorPalette);
    legH(i) = scatter(dn2dt(t(sI)),-eqdepth(sI),magFact*exp(eqmag(sI)),-eqdepth(sI),allMarkers{i},...
        'filled','MarkerEdgeColor','k','MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7);
end

ax4(2).XTickLabelRotation = rotAngle;
ylabel(depthLabel);
ax4(2).YAxisLocation = "right";
if Ntot >= Nfilt
    legH(lall+1) = plot(ax4(2),dn2dt(medfiltSH(t,Nfilt)),medfiltSH(-eqdepth,Nfilt),'k');
end
zoom on;
grid on;
%whos
%lall
%lall+1
%legH
%legendLabel
legend(legH(lall),legendLabel,"Location","NorthWest");
linkaxes(ax4,'x');
pause(1);

%% FIGURE 05
[~,meanmag,medmag,sumenergy] = t2r(dn2dt(t),days(nDays),eqmag);
fig(5) = figure('units','normalized','outerposition',[0.1 0.05 0.55 0.9]);
disp(['Length of eqmag: ',num2str(length(eqmag))]);

ax(1) = subplot(211);
plot(dn2dt(t),meanmag,'.'); hold on;
plot(dn2dt(t),medmag,'o'); grid on;
legend(ax(1),"Mean Magnitude","Median Magnitude","Location","NorthEast");
title(sprintf("Mean and Median Magnitude over %d-day window",nDays));
pause(1);

ax(2) = subplot(212);
semilogy(ax(2),dn2dt(t),sumenergy/nDays,'.');
linkaxes(ax,'x');
zoom on; grid on;
title(sprintf("Daily Energy Release [Joules] averaged over %d-day window",nDays));
pause(1);

%%
t = dn2dt(t);
pos1 = ax2_1.Position;
pos3 = ax2_3.Position;
stretch3 = ((pos3(4)./pos1(4))-1)/2;
ylim_old = ax2_3.YLim;
ylim_new = ylim_old;
ylim_new(2) = ylim_old(1) + (1+stretch3)*range(ylim_old);
ylim_new(1) = ylim_old(2) - (1+stretch3)*range(ylim_old);
ax2_3.YLim = ylim_new;

%% print files
printFlag = true;
if printFlag
    fprintf('printing figures to jpg files...\n');
    parentDir = fullfile('~/quick_plots/',lower(regionName));
    if ~exist(parentDir,'dir')
        mkdir(parentDir);
    end
    cd(char(parentDir));

    %imageParent = ['/opt/lampp/htdocs/images/',regionName,'/',typeLabel];
    %imageParent = ['/Users/stephen/recent_plots/',regionName,'/',typeLabel];
    %imageParent = ['/home/shernandez/Desktop/mario/',regionName,'/',typeLabel];
    %imageParent = strcat('~/public_html/prueba/images/',regionName,'/',typeLabel);
    imageParent = fullfile(parentDir,typeLabel);
    if ~exist(imageParent,'dir')
        disp('Directory doesnt exist. Creating...')
        mkdir(imageParent)
    end

    totFigs = length(fig);
    for i = totFigs:-1:1
        disp(i);
        %fname = strcat(imageParent,'/0',num2str(i),'_',regionName,'_',typeLabel);
        fname = strcat(imageParent,'/0',num2str(i),'_',regionName,'_',datestr(tStart,'yyyy_mm_'),typeLabel);
        print(fig(i),'-djpeg',fname);
    end
end
