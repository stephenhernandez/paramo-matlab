%function weeklySeismicityStatistics()
clear; close all;
load('~/igdata/ec_boundaries.mat','lonEC','latEC');
load('~/igdata/soam_noec.mat','lat_noec','lon_noec');

vFlag = true;
magFact = 6;
minMag = 1;
MarkerFaceAlpha = 0.3;
maxRMS = 2;
maxAzGap = 180;
cbarLocation = 'eastoutside';
cmap1 = pink(512);
MarkerEdgeColor = 'k';
ScatterLineWidth = 0.2;
maxDepth = 60;
minNMLV = 5;
printFlag = 0;

%% continental ecuador
boundaryBox = getRegionSpatialDimensions('ecuador-continental');
minLon = boundaryBox(1);
maxLon = boundaryBox(2);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);

%%
screenSize = get(0,'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);

lr = maxLon - minLon;
ud = maxLat - minLat;
initFigSize = [screenHeight/screenWidth 1];
newFigSize = initFigSize.*[lr ud];
newFigSize = newFigSize/max(newFigSize);

fig = figure('units','normalized','outerposition',[0 0 [newFigSize]]);
laxes = 1;
ax = newplot();
h = gcf;
if ~vFlag
    h.Visible = 'off';
end
disp(ax(1).TightInset);
[lon,lat,demData] = cutDEM(boundaryBox);

minElev = min(min(demData));
maxElev = max(max(demData));
squareDegs = ud*lr;

if squareDegs > 10
    zinc = 200;
elseif squareDegs > 1
    zinc = 80;
else
    zinc = 20;
end
zlevs = (minElev:zinc:maxElev)';
[~,hcontour] = contour(ax(laxes),lon,lat,demData,zlevs,'linewidth',1);
colormap(cmap1);
axis(ax(laxes),'xy');

axis(ax(laxes),'equal');
cb1 = colorbar(ax(laxes),cbarLocation);
cb1.Visible = 'off';
axis(ax(laxes),boundaryBox);

hold(ax(laxes),'on');

hcoast = plot(ax(laxes),lonEC,latEC,'k-','linewidth',2);
axis(ax(laxes),[minLon maxLon minLat maxLat]);
hsoam = plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,'color',[0.5 0.5 0.5]);
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);

%%
tRef = dn2dt(ceil(now)):-1:datetime(2019,10,18); %datetime(2019,10,18):-1:datetime(2011,05,20); %datetime(2019,10,15)];
for i = 1%:length(tRef)
    disp(tRef(i));
    tEnd = tRef(i);
    tStart = tRef(i)-7;
    disp(tStart);
    fName = strcat('~/products/summaries/EcuadorWeeklySeismicity_',datestr(tStart,'yyyymmdd'));
    [t,eqlat,eqlon,eqdepth,eqmag,ids,rmssec,azgap,~,nMLv] = readCat1(tStart,tEnd,minMag);
    
    cd ~;
    tI = eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon & rmssec <= maxRMS & (azgap <= maxAzGap | nMLv >= minNMLV);
    idsOrig = ids;
    t = t(tI);
    eqlat = eqlat(tI);
    eqlon = eqlon(tI);
    eqdepth = eqdepth(tI);
    eqmag = eqmag(tI);
    ids = ids(tI);
    rmssec = rmssec(tI);
    azgap = azgap(tI);
    
    %%
    laxes = 2;
    ax(laxes) = axes;
    S = scatter(ax(laxes),eqlon,eqlat,magFact*exp(eqmag),eqdepth,'o','filled');
    zoom on;
    S.MarkerFaceAlpha = MarkerFaceAlpha;
    S.MarkerEdgeColor = MarkerEdgeColor;
    S.LineWidth = ScatterLineWidth;
    c = colorbar(ax(laxes),cbarLocation);
    axis(ax(1:laxes),'equal');
    
    
    colormap(ax(laxes),parula);
    c.Label.String = 'depth [km.]';
    c.Label.Interpreter = 'latex';
    caxis([0 maxDepth]);
    h = title(ax(1),strcat(datestr(tStart,'yyyy/mm/dd'),' - ',datestr(tEnd,'yyyy/mm/dd')));
    h2 = xlabel(ax(1),'Longitude');
    h3 = ylabel(ax(1),'Latitude');
    
    %%
    ax(laxes).Visible = 'off';
    axis(ax(1:laxes),[minLon maxLon minLat maxLat]);
    linkaxes(ax(1:end),'xy');
    
    %%
    print('-dpng',fName,'-r70');
    if i > 1
        delete(ax(2));
        delete(h);
    end
end

