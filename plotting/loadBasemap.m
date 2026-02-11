function [fig,boundaryBox,Zlimits] = loadBasemap(boundaryBox,cmap,urbanFlag,overwriteFlag)
if nargin < 2; cmap = []; end
if nargin < 3; urbanFlag = false; end
if nargin < 4; overwriteFlag = false; end

%%
boundaryBoxOrig = boundaryBox;
if isstring(boundaryBoxOrig) || ischar(boundaryBoxOrig)
    boundaryBox = getRegionSpatialDimensions(boundaryBoxOrig);
end
minLon = boundaryBox(1);
maxLon = boundaryBox(2);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);

%%
linewidth = 0.5;
if isempty(cmap)
    pinkcmap = pink(256);
    cmap = demcmap([-5476; 5957],256,(winter(64)),(pinkcmap(128:end,:)));
end

if isempty(urbanFlag)
    urbanFlag = false;
end

%%
if isempty(minLat)
    disp('something went wrong');
    fig = gobjects(1);
    return;
end

%%
disp(boundaryBox');
roundN = 2;
minLon = round(minLon,roundN);
maxLon = round(maxLon,roundN);
minLat = round(minLat,roundN);
maxLat = round(maxLat,roundN);
cbarLocation = 'eastoutside';

if urbanFlag
    thisName = strcat("~/igdata/basemaps/basemap_",num2str(abs(minLat)),"_",num2str(abs(maxLat)),"_",...
        num2str(abs(minLon)),"_",num2str(abs(maxLon)),"_urbanTrue.fig");
else
    thisName = strcat("~/igdata/basemaps/basemap_",num2str(abs(minLat)),"_",num2str(abs(maxLat)),"_",...
        num2str(abs(minLon)),"_",num2str(abs(maxLon)),"_urbanFalse.fig");
end

if isfile(thisName) && ~overwriteFlag
    disp('file already exists, opening it...')
    fig = open(char(thisName));
    Zlimits = [];
    return;
end

%% doesnt exist, so make it
load('~/igdata/ec_boundaries.mat','lonEC','latEC');
load('~/igdata/soam_noec.mat','lat_noec','lon_noec');

%%
screenSize = get(0,'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

lr = maxLon - minLon;
ud = maxLat - minLat;
initFigSize = [screenHeight/screenWidth 1];
newFigSize = initFigSize.*[lr ud];
newFigSize = newFigSize/max(newFigSize);

fig = figure('units','normalized','outerposition',[0 0 newFigSize]);
laxes = 1;
ax = newplot();
fig.Visible = 'off';

if isstring(boundaryBoxOrig) || ischar(boundaryBoxOrig)
    [lon,lat,demData] = loadNamedDEM(boundaryBoxOrig);
    demData= demData';
else
    [lon,lat,demData] = cutDEM(boundaryBoxOrig);
end

%demData(demData>0) = NaN; %experimental
minElev = floor(min(min(demData)));
maxElev = ceil(max(max(demData)));
Zlimits = [minElev; maxElev];
squareDegs = ud*lr;

if squareDegs > 10
    zinc = 300;
elseif squareDegs > 1
    zinc = 150;
else
    zinc = 75;
end

zlevs = (minElev:zinc:maxElev)';
zlevs(~isfinite(zlevs));
contour(ax(laxes),lon,lat,demData',zlevs,'linewidth',linewidth);
colormap(cmap);
axis(ax(laxes),'xy');

axis(ax(laxes),'equal');
cb1 = colorbar(ax(laxes),cbarLocation);
cb1.Visible = 'off';
axis(ax(laxes),[minLon maxLon minLat maxLat]);

hold(ax(laxes),'on');

plot(ax(laxes),lonEC,latEC,'k-','linewidth',2);
axis(ax(laxes),[minLon maxLon minLat maxLat]);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,'color',[0.5 0.5 0.5]);
%geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);
zoom on;

% plot urban are
if urbanFlag
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
end

% now sav
savefig(fig,thisName); %,'compact');

