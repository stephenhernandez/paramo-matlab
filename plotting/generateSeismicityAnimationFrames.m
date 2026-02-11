function [t,eqlat,eqlon,eqdepth,eqmag,id,D,ax] = generateSeismicityAnimationFrames(varargin)
% Example: [t,eqlat,eqlon,eqdepth,eqmag,id,D] =
% generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,SC5Flag,minMag
% depthPanel,cumPanel,maxDepth,addSimbologia,allFrames,hillShadeFlag,magFact,fontSize,fstring,minLevInc);

% define defaults
functionDefaults = {...
    datetime('now'),...       % start time
    datetime('now'),...        % end time
    'chiles',...                % regionName
    [],...                      % boundaryBox
    true,...                    % SC5Flag
    2,...                       % min mag
    true,...                    % plot depth panel?
    false,...                   % plot cum panel?
    0,...                       % min depth
    20,...                      % max depth
    true,...                    % add legend?
    false,...                   % all frames?
    false,...                   % hillshading?
    20,...                      % mag amplification factor
    18,...                      % text font size
    '-djpeg',...                % output format
    [],...                      % minLevInc
    false,...                   % diasFlag
    [16 21],...                 % depthRange
    0};                         % depth correction

% deal variables
optsToUse = functionDefaults;
nVarargin = length(varargin);
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,regionName,boundaryBox,SC5Flag,minMag,depthPanel,cumPanel,...
    minDepth,maxDepth,addSimbologia,allFrames,hillShadeFlag,magFact,...
    fontSize,fstring,minLevInc,diasFlag,depthRange,refVMdepth] = ...
    deal(optsToUse{:});

%% execute
if allFrames
    dayNoExist = ~exist('~/animations','dir');
    if dayNoExist
        mkdir('~/animations');
    end
end

%%
if ~isempty(boundaryBox)
    minLon = boundaryBox(1);
    maxLon = boundaryBox(2);
    minLat = boundaryBox(3);
    maxLat = boundaryBox(4);
    centerLat = mean([maxLat minLat]);
    centerLon = mean([maxLon minLon]);
else
    [centerLat,centerLon,minLat,maxLat,minLon,maxLon] = get_region_dimensions(regionName);
end

if SC5Flag
    % read all data, filter later
    if diasFlag
        load('~/research/now/sierra_negra/SierraNegraCatalog30Aug2019_2','t','eqlat','eqlon','eqdepth','eqmag','id');
        [t,sortI] = sort(t);
        eqlat = eqlat(sortI);
        eqlon = eqlon(sortI);
        eqdepth = eqdepth(sortI);
        eqmag = eqmag(sortI);
        id = id(sortI);
        clear sortI;
    else
        % [evDescription,evStatus,t,eqlat,eqlon,eqdepth,~,id,...
        %     ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,eqmag] = readCat1(tStart);
        minMLv = 0; %7;
        nPhases = -2; %12;
        [~,~,t,eqlat,eqlon,eqdepth,~,id,~,~,phaseTot,nMLv,~,~,~,~,~,~,~,~,...
            ~,~,~,~,eqmag] = readCat1(tStart);
        %disp(unique(evType))
    end
else
    % read all data, filter later
    [t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD(regionName);
    table(t,eqlat,eqlon,eqdepth,eqmag,id)
    %refVMdepth = 0; %1.1 for SN, other value for other volcanoes; % 0 for iasp91
    eqdepth = eqdepth-refVMdepth;
end

%%
versionNumber = version('-release');
versionNumber = char(versionNumber);
versionNumber = versionNumber(1:4);
versionNumber = str2double(versionNumber);

timeColorCode = ~true;
LT = false;
if LT
    tStart = tStart - hours(5);
    tEnd = tEnd - hours(5);
    t = t - hours(5);
end

refEllipse = referenceEllipsoid('wgs84');
D = distance(eqlat,eqlon,centerLat,centerLon,refEllipse)*1e-3;

if strcmp(regionName,'cotopaxi')
    dI = t >= tStart & t < tEnd ...
        & eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon ...
        & eqdepth <= maxDepth & eqmag >= minMag;
else
    if ~SC5Flag
        dI = t >= tStart & t < tEnd ...
            & eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon ...
            & eqdepth <= maxDepth & eqdepth >= minDepth & eqmag >= minMag;
    else
        %serror = sqrt(eqlaterr.^2 + eqlonerr.^2 + eqdeptherr.^2);
        dI = t >= tStart & t < tEnd ...
            & eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon ...
            & eqdepth >= minDepth & eqdepth <= maxDepth & eqmag >= minMag & phaseTot >= nPhases & nMLv >= minMLv;% & eqmagerr < 1 & ...
        %eqmagerr > 0 & serror <= 50 & eqmag >= -3 & azgap < 300 & stderr <= 1.5;
        %lI = eqmagerr < 1 & eqmagerr > 0 & eqlaterr <= 15 & eqlonerr <= 15 & eqmag >= 3 & azgap < 300;
    end
end

dI = find(dI);
lDI = length(dI);
if ~allFrames
    %put smaller magnitudes on top of larger magnitudes
    [~,sI] = sort(eqmag(dI),'descend');
    dI = dI(sI);
    sIorig = sI;
    %clear sI;
end

[stnms,~,stla,stlo] = stationListFromPoint(centerLat,centerLon);
dIs = stlo >= minLon & stlo <= maxLon & stla >= minLat & stla <= maxLat; % & length(stnms) < 5;
stnms = stnms(dIs);
stla = stla(dIs);
stlo = stlo(dIs);

% disp(unique(evType(dI)))
% disp(unique(evStatus(dI)))

%% slab stuff
if strcmp(regionName,'esmeraldas') || strcmp(regionName,'bahia') || ...
        strcmp(regionName,'plata')
    plotSlabContours = true;
else
    plotSlabContours = false;
end

if plotSlabContours
    load ~/igdata/ecuador_slab_model.mat
    latSlabI = latSlab >= minLat & latSlab <= maxLat;
    slabDepth = slabDepth(latSlabI,:);
    latSlab = latSlab(latSlabI);
    nans = nans(latSlabI,:);

    lonSlabI = lonSlab >= minLon & lonSlab <= maxLon;
    slabDepth = slabDepth(:,lonSlabI);
    lonSlab = lonSlab(lonSlabI);
    nans = nans(:,lonSlabI);
end

%%
if ismember('TERV',stnms)
    [~,locb] = ismember('TERV',stnms);
    stnms(locb) = [];
    stla(locb) = [];
    stlo(locb) = [];
end

%%
removeI = [];
for kk = 1:length(stla)
    stnm_ = char(stnms(kk));
    if strcmp(stnm_(1:2),'SN')
        removeI = [removeI; kk];
    end
end
stnms(removeI) = [];
stla(removeI) = [];
stlo(removeI) = [];

%%
removeI = [];
for kk = 1:length(stla)
    stnm_ = char(stnms(kk));
    if strcmp(stnm_(1:2),'GV') || strcmp(stnm_(1:2),'EC') || strcmp(stnm_(1:end),'BUCE')
        removeI = [removeI; kk];
    end
end
stnms(removeI) = [];
stla(removeI) = [];
stlo(removeI) = [];

%%
maxMag = ceil(max(eqmag(dI)))+0.5;
if depthPanel
    maxYLim = 10*ceil(max(eqdepth(dI))/10);
    minYLim = floor(min(eqdepth(dI)));
    if minYLim < 0
        minYLim = -5;
    else
        minYLim = 0;
    end
    legLoc = 'SouthWest';
    logFlag = false;
elseif cumPanel
    legLoc = 'NorthEast';
    logFlag = false;
else
    maxYLim = 10*ceil(max(D(dI))/10);
    minYLim = floor(log10(min(D(dI))));
    minYLim = 10^minYLim;
    legLoc = 'NorthEast';
    logFlag = false;
end

% load important data
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);

% some parameters for hillshade option
ilumAz = -45;
wMark = 1.65;
elevAngle = 70;
demContrast = 2; %7/8;

%% commence figure making
figNumber = 100;
laxes = 1;

fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);
subplot(2,2,[1 3],ax(laxes));

geoshowFlag = true;
if hillShadeFlag
    disp('applying hillshading')
    I = dem(lon,lat,demData','Contrast',demContrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1]);
    imagesc(ax(laxes),I.x,I.y,I.rgb);
    hold(ax(laxes),'on');
    axis(ax(laxes),'xy');

    if strcmp(regionName,'esmeraldas')

        contour(ax(laxes),lon,lat,demData',(-(50:100:4000))','-','linewidth',0.1);
        cbar12 = colorbar(ax(laxes));
        colorbarRightShift = 0.93; 
        %cbar12.Position = [colorbarRightShift,0.55,0.0083333333333333,0.35];
        %axis(ax(laxes),'equal');
        %hold(ax(laxes),'on');
        %axis(ax(laxes),[minLon maxLon minLat maxLat]);
        bonemap = flipud(bone(512));
        colormap(ax(laxes),bonemap);
        % ax(laxes).Visible = 'off';
        cbar12.Visible = 'off';
        % hold(ax(laxes),'on');
        % axis(ax(laxes),'equal');
    end

    if strcmp(regionName,'guayaquil') || ...
            strcmp(regionName,'pululahua') || strcmp(regionName,'fernandina') ...
            || strcmp(regionName,'sierra_negra') || strcmp(regionName,'tungurahua') ...
            || strcmp(regionName,'wolf') || strcmp(regionName,'cerro_azul') ...
            || strcmp(regionName,'alcedo') || strcmp(regionName,'darwin') ...
            || strcmp(regionName,'bahia')
        minElev = min(min(demData));
        maxElev = max(max(demData));

        if geoshowFlag
            elevRange = maxElev - minElev;
            fprintf("Elevation Range in DEM: %d\n",elevRange);
            if elevRange > 1200
                disp('generating 400m contour levels');
                [~,cc] = contour(ax(laxes),lon,lat,demData',(minElev:400:maxElev)');
                cc.LineWidth = 0.04;
                cc.Color = 'k'; %[0.5 0.5 0.5];
            else
                disp('generating 200m contour levels');
                [~,cc] = contour(ax(laxes),lon,lat,demData',(minElev:200:maxElev)');
                cc.LineWidth = 0.04;
                cc.Color = 'k';
            end
        end
    end

    if plotSlabContours
        zInc = 3;
        zlevs = 1000*(ceil(min(min(slabDepth)/1000)):zInc:floor(max(max(slabDepth)/1000)));
        for iii = 1:length(zlevs)
            Ctmp = contour(lonSlab,latSlab,slabDepth,[zlevs(iii) zlevs(iii)]);
            [~,~,xstart,ystart] = plotContourMatrix(ax(laxes),Ctmp);
            HH = text(xstart(1),ystart(1),[num2str(zlevs(iii)/1000),' km.'],'FontSize',15);
            for jj = 1:length(HH)
                HH(jj).Rotation = 63; %.Interpreter = 'none';
            end
        end
    end
    %plot(ax(laxes),stlo,stla,'v','MarkerEdgeColor','k','MarkerFaceColor','w');
    %text(ax(laxes),stlo,stla,stnms,'FontSize',15);
    % if strcmp(regionName,'esmeraldas')
    %     PWD = pwd;
    %     cd ~/research/now/esmeraldas/shapefiles/rupture-slowslip/
    %     %close all; 
    %     %fig = figure();
    %     files = dir('*.shp');
    %     lFiles = length(files);
    %     for i = 4:lFiles-1
    %         T = readgeotable(files(i).name);
    %         try
    %             hold on; 
    %             H = geoshow(T);
    %         catch
    %             try
    %                 hold on; 
    %                 H = mapshow(T); 
    %             catch
    %                 continue; 
    %             end
    %         end
    %     end
    %     axis equal
    %     cd(PWD);
    % end
else
    minElev = min(min(demData));
    maxElev = max(max(demData));

    if isempty(minLevInc)
        elevRange = maxElev - minElev;
        fprintf("Elevation Range in DEM: %d\n",elevRange);
        if elevRange > 1000
            [~,cc] = contour(ax(laxes),lon,lat,demData',(minElev:80:maxElev)');
            cc.LineWidth = 0.03;
            %cc.Color = [0.5 0.5 0.5];
        else
            [~,cc] = contour(ax(laxes),lon,lat,demData',(minElev:30:maxElev)');
            cc.LineWidth = 0.03;
            %cc.Color = [0.5 0.5 0.5];
        end
    else
        if geoshowFlag
            elevRange = maxElev - minElev;
            fprintf("Elevation Range in DEM: %d\n",elevRange);
            [~,cc] = contour(ax(laxes),lon,lat,demData',(minElev:minLevInc:maxElev)');
            cc.LineWidth = 0.01;
            %cc.Color = [0.5 0.5 0.5];
        end
    end
    hold(ax(laxes),'on');
    axis(ax(laxes),'xy');
end

%%
stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);

ylabel(ax(laxes),'Latitud','FontSize',fontSize);
xlabel(ax(laxes),'Longitud','FontSize',fontSize);
axis(ax(laxes),'equal');

if geoshowFlag
    geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);
end

if strcmp(regionName,'fernandina') || strcmp(regionName,'galapagos')
    S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2018.shp','UseGeoCoords',true);
    [slat,slon] = utm2ll(S.Lon,S.Lat,15);
    plot(ax(laxes),slon,slat,'r','linewidth',2); zoom on; grid on;
    S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2017.shp','UseGeoCoords',true);
    [slat,slon] = utm2ll(S.Lon,S.Lat,15);
    %figure(1); hold on;
    plot(ax(laxes),slon,slat,'r','linewidth',2); zoom on; grid on;
    axis equal;
    S = shaperead('~/research/now/fernandina/Fernandina/Fernandina2020_3mar','UseGeoCoords',true);
    for i = 1:length(S)
        [slat,slon] = utm2ll(S(i).Lon,S(i).Lat,-15);
        badI = find(~isfinite(slon))';
        if ~sum(badI)
            fill(slon,slat,'r'); zoom on; grid on;
        else
            si = [1; badI(1:end-1) + 1]; ei = badI - 1;
            for j = 1:length(si)
                fill(slon(si(j):ei(j)),slat(si(j):ei(j)),'r'); zoom on; grid on;
            end
        end
    end
    %plot(ax(laxes),stlo,stla,'kv','linewidth',2,'MarkerFaceColor','w'); zoom on; grid on;
end

if strcmp(regionName,'sierra_negra') || strcmp(regionName,'galapagos')
    S = shaperead('~/research/now/sierra_negra/sierraNegraLavaFlowSHP/lava_SN2018_31jul.shp','UseGeoCoords',true);
    for i = 1:length(S)
        [slat,slon] = utm2ll(S(i).Lon,S(i).Lat,-15);
        badI = find(~isfinite(slon))';
        if ~sum(badI)
            fill(slon,slat,'r'); zoom on; grid on;
        else
            si = [1; badI(1:end-1) + 1]; ei = badI - 1;
            for j = 1:length(si)
                fill(slon(si(j):ei(j)),slat(si(j):ei(j)),'r'); zoom on; grid on;
            end
        end
    end
    %plot(ax(laxes),stlo,stla,'kv','linewidth',2,'MarkerFaceColor','w'); zoom on; grid on;
end

if strcmp(regionName,'wolf') || strcmp(regionName,'galapagos')
    S = shaperead('~/research/now/wolf/Wolf_LF_25May22.shp','UseGeoCoords',true);
    for i = 1:length(S)
        [slat,slon] = utm2ll(S(i).Lon,S(i).Lat,-15);
        badI = find(~isfinite(slon))';
        if ~sum(badI)
            fill(slon,slat,'r'); zoom on; grid on;
        else
            si = [1; badI(1:end-1) + 1]; ei = badI - 1;
            for j = 1:length(si)
                fill(slon(si(j):ei(j)),slat(si(j):ei(j)),'r'); zoom on; grid on;
            end
        end
    end
    %plot(ax(laxes),stlo,stla,'kv','linewidth',2,'MarkerFaceColor','w'); zoom on; grid on;
end

if strcmp(regionName,'chiles')
    load('~/research/now/chiles/PuntosTermas');
    plot(ptlon,ptlat,'d','MarkerEdgeColor','k','MarkerFaceColor','r');
    text(ptlon,ptlat,ptstnm);
end

% if geoshowFlag
%     geoshow(ax,'~/igdata/ZonaUrbana/ZonaUrbana.shp');
% end

% if strcmp(regionName,'esmeraldas')
%     text(ax(laxes),-79.61,+0.96,'Esmeraldas','fontsize',fontSize+0,'FontName','Castellar');
%     text(ax(laxes),-79.8,+0.88,'Atacames/','fontsize',fontSize-1,'FontName','Castellar');
%     text(ax(laxes),-79.8,+0.85,'Tonsupa','fontsize',fontSize-1,'FontName','Castellar');
% end
%S = shaperead('~/igdata/volcanes/volcanes_2.shp');
%plot(ax(laxes),[S.X],[S.Y],'k--','linewidth',1);

cb1 = colorbar(ax(laxes));
cb1.Visible = 'off';
axis(ax(laxes),[minLon maxLon minLat maxLat]);

pointColormap = 'parula';
if strcmp(pointColormap,'turbo')
    if dn2dt(datenum(version('-date'))) < datetime(2018,01,01)
        pointColormap = 'parula';
    end
end

cmap = 'copper';
if ~hillShadeFlag
    colormap(ax(laxes),cmap);
end

if allFrames
    lDI = 1:lDI;
end
markerFaceAlpha = 0.75;

fprintf("Total number of frames to produce: %d\n",max(lDI));
for i = lDI
    fprintf("frame number: %d\n",i);
    ax(laxes+1) = axes;
    ax(laxes+1).Parent = fig(figNumber);
    subplot(2,2,[1 3],ax(laxes+1));
    if timeColorCode
        S = scatter(ax(laxes+1),eqlon(dI(1:i)),eqlat(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled');
        clim(ax(laxes+1),datenum([tStart tEnd]));
        c = colorbar;
        c.Label.Interpreter = 'latex';
        c.TickLabels = datestr(c.Ticks);
        colormap(ax(laxes+1),pointColormap);
    else
        S = scatter(ax(laxes+1),eqlon(dI(1:i)),eqlat(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled');
        %[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]
        %ceil(max(eqdepth(dI)))
        %min([floor(min(eqdepth(dI))) 0])

        if isempty(depthRange)
            %ceil(max(eqdepth(dI)))
            %floor(min(eqdepth(dI)))
            clim(ax(laxes+1),-[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]);
        else
            clim(ax(laxes+1),-depthRange);
        end
        c = colorbar;
        c.Label.Interpreter = 'latex';
        colormap(ax(laxes+1),pointColormap);
        disp(-[maxDepth ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0]) min(eqdepth(dI))]);
    end

    hold(ax(laxes+1),'on');
    laxes
    axis(ax(laxes+1),'equal');

    %text(ax(laxes+1),stlo+0.001,stla+0.001,stnms,'fontsize',fontSize,'FontName','Castellar');
    if LT
        titleStr = [datestr(t(dI(i))),' (Tiempo Local), N=',num2str(i)];
    else
        i
        dI(i)
        t(dI(i))
        titleStr = sprintf('%s (UTC), N=%d',t(dI(i)),i); %[datestr(t(dI(i))),' (UTC), N=',num2str(i)];
    end

    h = title(ax(laxes),titleStr,'FontSize',fontSize);
    S.MarkerFaceAlpha = markerFaceAlpha;
    S.MarkerEdgeColor = 'k';
    S.MarkerEdgeAlpha = markerFaceAlpha;

    if addSimbologia
        if strcmp(regionName,'chiles')
            hl = zeros(8,1);
            hl(2) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(1))); %,NaN,'filled');
            hl(4) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(2))); %,NaN,'filled');
            hl(6) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(3))); %,NaN,'filled');
            hl(8) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(4))); %,NaN,'filled');
            hl(1) = plot(ax(laxes),NaN,NaN,'v','markersize',15,'markeredgecolor',stationOutlineColor,'markerfacecolor','w','linewidth',1); %,NaN,'filled');
            hl(3) = plot(ax(laxes),[NaN,NaN],[NaN,NaN],'k-','linewidth',2);
            hl(5) = plot(ax(laxes),[NaN,NaN],[NaN,NaN],'k--','linewidth',1);
            hl(7) = plot(ax(laxes),[NaN,NaN],[NaN,NaN],'r-','linewidth',2);
            legendTexts = {'Sismometro','M1','Frontera','M2','Volcanes','M3','Fallas','M4'};
            lgd = legend(hl,legendTexts,'location',legLoc);
            title(lgd,'Simbologia');
            lgd.Orientation = 'horizontal';

            if versionNumber >= 2018
                lgd.NumColumns = 2;
            end
        end

        if strcmp(regionName,'sierra_negra')
            hl = zeros(6,1);
            hl(2) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(2))); %,NaN,'filled');
            hl(3) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(4))); %,NaN,'filled');
            hl(5) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(3))); %,NaN,'filled');
            hl(6) = plot(ax(laxes),NaN,NaN,'ko','markersize',sqrt(magFact*exp(5))); %,NaN,'filled');
            hl(1) = plot(ax(laxes),NaN,NaN,'v','markersize',15,'markeredgecolor',stationOutlineColor,'markerfacecolor','w','linewidth',1); %,NaN,'filled');
            hl(4) = plot(ax(laxes),[NaN,NaN],[NaN,NaN],'k-','linewidth',2);

            legendTexts = {'Sismometro','M2','M4','Costa','M3','M5'};
            lgd = legend(hl,legendTexts,'location',legLoc);
            title(lgd,'Simbologia');

            if versionNumber >= 2018
                lgd.NumColumns = 2;
            end
        end
    end

    linkaxes(ax(laxes:laxes+1));
    ax(laxes+1).Visible = 'off';
    axis(ax(laxes+1),[minLon maxLon minLat maxLat]);

    if i < 10
        fname = strcat('~/animations/frame_000',num2str(i));
    elseif i < 100
        fname = strcat('~/animations/frame_00',num2str(i));
    elseif i < 1000
        fname = strcat('~/animations/frame_0',num2str(i));
    else
        fname = strcat('~/animations/frame_',num2str(i));
    end

    % upper right panel
    if depthPanel
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));
        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),-eqdepth(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            clim(ax(laxes+2),datenum([tStart tEnd]));
            colormap(ax(laxes+2),pointColormap);
            grid on;
        else
            S = scatter(ax(laxes+2),t(dI(1:i)),-eqdepth(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            if isempty(depthRange)
                clim(ax(laxes+2),-[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]);
            else
                clim(ax(laxes+2),-depthRange);
            end
            colormap(ax(laxes+2),pointColormap);
            disp(-[maxDepth min([floor(min(eqdepth(dI))) 0]) min(eqdepth(dI))]);
        end

        c = colorbar;
        xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        disp(-[maxYLim minYLim]);
        %ylim(ax(laxes+2),-[maxYLim minYLim]);
        ylabel(ax(laxes+2),'Profundidad [km]','FontSize',fontSize);
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    elseif cumPanel
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));

        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),1:i,magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            clim(ax(laxes+2),datenum([tStart tEnd]));
            colormap(ax(laxes+2),pointColormap);
        else
            if ~allFrames
                %[~,sI] = sort(t(dI));
                cumDum = (1:length(t(dI)))';
                S = scatter(ax(laxes+2),t(dI),cumDum(sIorig),magFact*exp(eqmag(dI)),-eqdepth(dI),'filled','linewidth',0.2,'markeredgecolor','k');
            else
                S = scatter(ax(laxes+2),t(dI(1:i)),1:i,magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            end
            if isempty(depthRange)
                clim(ax(laxes+2),-[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]);
            else
                clim(ax(laxes+2),-depthRange);
            end
            colormap(ax(laxes+2),pointColormap);
            disp(-[maxDepth min([floor(min(eqdepth(dI))) 0]) min(eqdepth(dI))]);
        end

        c = colorbar;
        xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        %ylim(ax(laxes+2),[0 max(lDI)+1]);
        ylabel(ax(laxes+2),'Numero Acumulativo','FontSize',fontSize);
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    else
        ax(laxes+2) = axes;
        ax(laxes+2).Parent = fig(figNumber);
        subplot(2,2,2,ax(laxes+2));
        if timeColorCode
            S = scatter(ax(laxes+2),t(dI(1:i)),D(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
            clim(ax(laxes+2),datenum([tStart tEnd]));
            colormap(ax(laxes+2),pointColormap);
        else
            S = scatter(ax(laxes+2),t(dI(1:i)),D(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
            if isempty(depthRange)
                clim(ax(laxes+2),-[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]);
            else
                clim(ax(laxes+2),-depthRange);
            end
            colormap(ax(laxes+2),pointColormap);
            disp(-[maxDepth min([floor(min(eqdepth(dI))) 0]) min(eqdepth(dI))]);
        end

        c = colorbar;
        xlim(ax(laxes+2),[tStart tEnd]);
        if logFlag
            ax(laxes+2).YScale = 'log';
        end
        ylim(ax(laxes+2),[minYLim maxYLim]);
        ylabel(ax(laxes+2),'Distancia a Diamante Rojo [km]','FontSize',fontSize);
        c.Visible = 'off';
        S.MarkerFaceAlpha = markerFaceAlpha;
    end

    % lower right panel
    ax(laxes+3) = axes;
    ax(laxes+3).Parent = fig(figNumber);
    subplot(2,2,4,ax(laxes+3));
    if timeColorCode
        S = scatter(ax(laxes+3),t(dI(1:i)),eqmag(dI(1:i)),magFact*exp(eqmag(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
        clim(ax(laxes+3),datenum([tStart tEnd]));
        colormap(ax(laxes+3),pointColormap);
        grid on;
    else
        S = scatter(ax(laxes+3),t(dI(1:i)),eqmag(dI(1:i)),magFact*exp(eqmag(dI(1:i))),-eqdepth(dI(1:i)),'filled','linewidth',0.2,'markeredgecolor','k');
        if isempty(depthRange)
            clim(ax(laxes+3),-[ceil(max(eqdepth(dI))) min([floor(min(eqdepth(dI))) 0])]);
        else
            clim(ax(laxes+3),-depthRange);
        end
        colormap(ax(laxes+3),pointColormap);
        disp(-[maxDepth min([floor(min(eqdepth(dI))) 0]) min(eqdepth(dI))]);
    end
    c = colorbar;
    xlim(ax(laxes+3),[tStart tEnd]);
    ylim(ax(laxes+3),[minMag maxMag]);
    ylabel(ax(laxes+3),'Magnitud','FontSize',fontSize);
    c.Visible = 'off';
    S.MarkerFaceAlpha = markerFaceAlpha;
    linkaxes(ax(laxes+(2:3)),'x');
    zoom on; grid on;
    ax(laxes+3).Box = 'on';

    zoom(ax(laxes+2),'on');
    grid(ax(laxes+2),'on');
    ax(laxes+2).Box = 'on';
    if allFrames
        disp(fname);
        print(fstring,fname);
        if i < max(lDI)
            delete(h);
            delete(ax(laxes+1));
            delete(ax(laxes+2));
            delete(ax(laxes+3));
            if addSimbologia
                if strcmp(regionName,'sierra_negra')
                    delete(lgd);
                end
            end
        end
    end
end

%%
t = t(dI);
eqlat = eqlat(dI);
eqlon = eqlon(dI);
eqdepth = eqdepth(dI);
eqmag = eqmag(dI);
id = id(dI);
D = D(dI);

if ~allFrames
    [t,dI] = sort(t);
    eqlat = eqlat(dI);
    eqlon = eqlon(dI);
    eqdepth = eqdepth(dI);
    eqmag = eqmag(dI);
    id = id(dI);
    D = D(dI);
end

%%
%plot(ax(laxes+1),stlo,stla,'kv','linewidth',1,'MarkerFaceColor','w'); zoom on; grid on;
%HH = text(ax(laxes+1),stlo-0.06,stla,stnms,"FontSize",13); %,'linewidth',0.1); %,'MarkerFaceColor','w')
% for jj = 1:length(HH)
%     HH(jj).Rotation = 60; %.Interpreter = 'none';
% end

