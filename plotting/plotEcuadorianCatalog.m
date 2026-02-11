function ax = plotEcuadorianCatalog(t,eqlat,eqlon,eqdepth,eqmag,varargin)
% Example: [] =
% generateSeismicityAnimationFrames(catalogData,boundaryBox,tStart,tEnd,regionName,demName,minMag,
% depthPanel,cumPanel,maxDepth,addSimbologia,allFrames,hillShadeFlag,magFact,fontSize,fstring,minLevInc);

% define defaults
functionOptionDefaults = {...
    [],...                      % boundaryBox
    [],...                      % tStart
    [],...                      % tEnd
    -2,...                       % minMag
    false,...                   % plot depth panel?
    true,...                    % plot cum panel?
    100,...                     % maxDepth
    false,...                   % hillshading?
    20,...                      % mag amplification factor
    18,...                      % text font size
    [],...                      % minLevInc
    false,...                   % timeColorCode
    false};                     % LT (local time?)

% deal variables
optsToUse = functionOptionDefaults;
nVarargin = length(varargin);
optsToUse(1:nVarargin) = varargin;
[boundaryBox,tStart,tEnd,minMag,depthPanel,cumPanel,maxDepth,hillShadeFlag,...
    magFact,fontSize,minLevInc,timeColorCode,LT] = deal(optsToUse{:});

%%
% t = catalogData(:,1);
% eqlat = catalogData(:,2);
% eqlon = catalogData(:,3);
% eqdepth = catalogData(:,4);
% eqmag = catalogData(:,5);
% id = catalogData(:,6);

if isempty(boundaryBox)
    minLon = floor(min(eqlon));
    maxLon = ceil(max(eqlon));
    minLat = floor(min(eqlat));
    maxLat = ceil(max(eqlat));
else
    minLon = boundaryBox(1);
    maxLon = boundaryBox(2);
    minLat = boundaryBox(3);
    maxLat = boundaryBox(4);
end

if isempty(tStart)
    tStart = dateshift(min(t),'start','day')
end

if isempty(tEnd)
    tEnd = dateshift(max(t),'end','day')
end

%%
centerLat = (maxLat - minLat)/2;
centerLon = (maxLon - minLon)/2;

refEllipse = referenceEllipsoid('wgs84');
D = distance(eqlat,eqlon,centerLat,centerLon,refEllipse)*1e-3;

dI = t >= tStart & t < tEnd ...
    & eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon ...
    & eqdepth <= maxDepth & eqmag >= minMag;

[stnms,~,stla,stlo] = stationListFromPoint(centerLat,centerLon);
dIs = stlo >= minLon & stlo <= maxLon & stla >= minLat & stla <= maxLat; % & length(stnms) < 5;
stnms = stnms(dIs);
stla = stla(dIs);
stlo = stlo(dIs);


%%
if LT
    tStart = tStart + hours(LT);
    tEnd = tEnd + hours(LT);
    t = t + hours(LT);
end

%%
if ismember('TERV',stnms)
    [~,locb] = ismember('TERV',stnms);
    stnms(locb) = [];
    stla(locb) = [];
    stlo(locb) = [];
end

%%
% removeI = [];
% for kk = 1:length(stla)
%     stnm_ = char(stnms(kk));
%     if strcmp(stnm_(1:2),'SN')
%         removeI = [removeI; kk];
%     end
% end
% stnms(removeI) = [];
% stla(removeI) = [];
% stlo(removeI) = [];

%%
removeI = [];
for kk = 1:length(stla)
    stnm_ = char(stnms(kk));
    if strcmp(stnm_(1:2),'GV')
        removeI = [removeI; kk];
    end
end
stnms(removeI) = [];
stla(removeI) = [];
stlo(removeI) = [];

maxMag = ceil(max(eqmag(dI)))+0.5;
if depthPanel
    maxYLim = 10*ceil(max(eqdepth(dI))/10);
    minYLim = floor(min(eqdepth(dI)));
    if minYLim < 0
        minYLim = -5;
    else
        minYLim = 0;
    end
    logFlag = false;
elseif cumPanel
    logFlag = false;
else
    maxYLim = 10*ceil(max(D(dI))/10);
    minYLim = floor(log10(min(D(dI))));
    minYLim = 10^minYLim;
    logFlag = false;
end

% load important data
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);

% some parameteres for hillshade option
ilumAz = -70;
wMark = 2;
elevAngle = 60;

%% commence figure making
figNumber = 100;
laxes = 1;

fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);
subplot(2,2,[1 3],ax(laxes));

if hillShadeFlag
    I = dem(lon,lat,demData,'Contrast',2,'Azimuth',ilumAz,'Interp','noplot','Watermark',wMark,'Elevation',elevAngle);
    imagesc(ax(laxes),I.x,I.y,I.rgb);
else
    minElev = min(min(demData));
    maxElev = max(max(demData));
    
    if isempty(minLevInc)
        elevRange = maxElev - minElev;
        if elevRange > 1000
            contour(ax(laxes),lon,lat,demData',(minElev:50:maxElev)');
        else
            contour(ax(laxes),lon,lat,demData',(minElev:50:maxElev)');
        end
    else
        contour(ax(laxes),lon,lat,demData',(minElev:minLevInc:maxElev)');
    end
end
hold(ax(laxes),'on');
axis(ax(laxes),'xy');

stationOutlineColor = 'k';
plot(ax(laxes),stlo,stla,'v','markeredgecolor',stationOutlineColor,'markerfacecolor','w','markersize',15,'linewidth',1);
plot(ax(laxes),lonEC,latEC,'k','linewidth',3);

ylabel(ax(laxes),'Latitud');
xlabel(ax(laxes),'Longitud');
axis(ax(laxes),'equal');
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',2);
S = shaperead('~/igdata/volcanes/volcanes_2.shp');
plot(ax(laxes),[S.X],[S.Y],'k--','linewidth',1);

cb1 = colorbar(ax(laxes));
cb1.Visible = 'off';
axis(ax(laxes),[minLon maxLon minLat maxLat]);

cmap = 'copper';
if ~hillShadeFlag
    colormap(ax(laxes),cmap);
end

dI = find(dI);
lDI = length(dI);
markerFaceAlpha = 0.5;

disp(strcat('Total number of frames to produce: ',num2str(max(lDI))));
for i = lDI
    disp(i);
    ax(laxes+1) = axes;
    ax(laxes+1).Parent = fig(figNumber);
    subplot(2,2,[1 3],ax(laxes+1));
    if timeColorCode
        S = scatter(ax(laxes+1),eqlon(dI(1:i)),eqlat(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.02,'markeredgecolor','k');
        caxis(ax(laxes+1),datenum([tStart tEnd]));
        c = colorbar;
        c.Label.Interpreter = 'latex';
        c.TickLabels = datestr(c.Ticks);
    else
        S = scatter(ax(laxes+1),eqlon(dI(1:i)),eqlat(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.02,'markeredgecolor','k');
        caxis(ax(laxes+1),[-5 2]);
        c = colorbar;
        c.Label.Interpreter = 'latex';
    end
    
    hold(ax(laxes+1),'on');
    axis(ax(laxes+1),'equal');
    
    text(ax(laxes+1),stlo+0.001,stla+0.001,stnms,'fontsize',fontSize,'FontName','Castellar');
    titleStr = [datestr(t(dI(i))),', N=',num2str(i)];
    %h =
    title(ax(laxes),titleStr);
    S.MarkerFaceAlpha = markerFaceAlpha;
    
    linkaxes(ax(laxes:laxes+1));
    ax(laxes+1).Visible = 'off';
    axis(ax(laxes+1),[minLon maxLon minLat maxLat]);
    
    % upper right panel
    if depthPanel
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));
        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),-eqdepth(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),datenum([tStart tEnd]));
        else
            S = scatter(ax(laxes+2),t(dI(1:i)),-eqdepth(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),[-5 2]);
        end
        
        c = colorbar;
        xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        disp(-[maxYLim minYLim]);
        %ylim(ax(laxes+2),-[maxYLim minYLim]);
        ylabel(ax(laxes+2),'Profundidad [km]');
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    elseif cumPanel
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));
        
        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),1:i,magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),datenum([tStart tEnd]));
        else
            S = scatter(ax(laxes+2),t(dI(1:i)),1:i,magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),[-5 2]);
        end
        
        c = colorbar;
        %xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        %ylim(ax(laxes+2),[0 max(lDI)+1]);
        ylabel(ax(laxes+2),'Numero Acumulativo');
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    else
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));
        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),D(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),datenum([tStart tEnd]));
        else
            S = scatter(ax(laxes+2),t(dI(1:i)),D(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            caxis(ax(laxes+2),[-5 2]);
        end
        
        c = colorbar;
        xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        ylim(ax(laxes+2),[minYLim maxYLim]);
        ylabel(ax(laxes+2),'Distancia a Diamante Rojo [km]');
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    end
    
    % lower right panel
    ax(laxes+3) = axes;
    ax(laxes+3).Parent = fig(figNumber);
    subplot(2,2,4,ax(laxes+3));
    if timeColorCode
        S = scatter(ax(laxes+3),t(dI(1:i)),eqmag(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
        caxis(ax(laxes+3),datenum([tStart tEnd]));
    else
        S = scatter(ax(laxes+3),t(dI(1:i)),eqmag(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
        caxis(ax(laxes+3),[-5 2]);
    end
    c = colorbar;
    ylim(ax(laxes+3),[minMag maxMag]);
    ylabel(ax(laxes+3),'Magnitud');
    c.Visible = 'off';
    S.MarkerFaceAlpha = markerFaceAlpha;
    linkaxes(ax(laxes+(2:3)),'x');
end
zoom on;

%%
% if nargout > 1
%     %% output filtered catalog
%     catalogData = [t(dI),eqlat(dI),eqlon(dI),eqdepth(dI),eqmag(dI),id(dI)];
% end
