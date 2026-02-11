function ax = plotEvents(E,visibleFlag,saveDir)
newDirFlag = false;
if nargin < 3
    saveDir = fullfile("~","masa","products","events","html");
    newDirFlag = true;
end

%%
MAXLON = -72;
MINLON = -92;
MAXLAT = 4;
MINLAT = -6;

%%
load('~/igdata/ec_boundaries.mat','lonEC','latEC');
load('~/igdata/soam_noec.mat','lat_noec','lon_noec');
lfiles = length(E);

figure('units','normalized','outerposition',[0 0 1 1]);
ax = newplot();
h = gcf;
if ~visibleFlag
    h.Visible = 'off';
end

for i = 1:lfiles
    disp(i);
    
    cd ~/igdata/dem/
    eqlon = E(i).lon;
    eqlat = E(i).lat;
    if eqlon >= MAXLON || eqlon <= MINLON || eqlat <= MINLAT || eqlat >= MAXLAT
        continue;
    end

    %%
    t = E(i).t;
    eqdepth = E(i).depth;
    eqmag = E(i).mag;
    id = E(i).id;
    thisSaveDir = saveDir;
    if newDirFlag
        thisSaveDir = strcat(saveDir,id);
    end
    disp(id);

    %%
    Pphases = E(i).Pphases;
    stnm_ = pull(Pphases,"stnm");
    uniqStnms = unique(stnm_);

    %%
    [stla,stlo] = metaDataFromStationList(uniqStnms);
    minLon = min([stlo; eqlon]) - 0.05;
    maxLon = max([stlo; eqlon]) + 0.05;
    minLat = min([stla; eqlat]) - 0.05;
    maxLat = max([stla; eqlat]) + 0.05;

    minLon = floor(10*max([minLon MINLON]))/10;
    maxLon = ceil(10*min([maxLon MAXLON]))/10;
    minLat = floor(10*max([minLat MINLAT]))/10;
    maxLat = ceil(10*min([maxLat MAXLAT]))/10;

    %%
    [lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
    demData = double(demData);

    contrast = 3/4;
    ilumAz = -45;
    wMark = 1.75;
    elevAngle = 45;

    I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'nodecim');
    I = I.rgb;
    imagesc(ax,lon,lat,I);

    %%
    axis xy;
    axis equal;
    %colorbar;
    hold on;

    LW = 1; 
    hepi = plot(ax,eqlon,eqlat,'p','linewidth',LW,'markerfacecolor','w','MarkerSize',24);
    hs = plot(ax,stlo,stla,'^','MarkerFaceColor','w','linewidth',LW,'MarkerSize',12);
    hcoast = plot(ax,lonEC,latEC,'k-','linewidth',2);
    xlabel(ax,'Longitude');
    ylabel(ax,'Latitude');

    topStr = string(strcat('$t_{0}$: ',datestr(t),', id: ',id));
    bottomStr = string(strcat('lat: ',num2str(eqlat),', lon: ',num2str(eqlon),', depth: ',num2str(eqdepth),', mag: ',num2str(eqmag)));
    title(ax,{topStr;bottomStr});

    axis(ax,[minLon maxLon minLat maxLat]);
    hsoam = plot(ax,lon_noec,lat_noec,'-','linewidth',2,'color',[0.5 0.5 0.5]);
    shape_dir = fullfile("~","igdata","shape_files"); 
    geoshow(fullfile(shape_dir,"fallas2008completas.shp"),'Color',[0.5, 0.5, 0.5],'linewidth',0.5);

    %%
    if ~exist(thisSaveDir,'dir')
        mkdir(thisSaveDir);
    end
    cd(thisSaveDir);
    print('-djpeg',id);

    %%
    if ~visibleFlag
        %delete(hcontour);
        delete(hcoast);
        delete(hs);
        delete(hepi);
        delete(hsoam);
    else
        zoom on;
    end
    ax.ColorOrderIndex = 1;
end
