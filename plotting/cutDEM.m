function [lon,lat,demData] = cutDEM(boundaryBox,verboseFlag)
if nargin < 2
    verboseFlag = false;
end

lon = [];
lat = [];
demData = [];

%%
boundaryBoxOrig = boundaryBox;
if isstring(boundaryBoxOrig) || ischar(boundaryBoxOrig)
    boundaryBox = getRegionSpatialDimensions(boundaryBoxOrig);
end

%%
minLon = boundaryBox(1);
maxLon = boundaryBox(2);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);

%%
if verboseFlag
    disp([minLon maxLon minLat maxLat]);
end

%%
PWD = pwd;
cd ~/igdata/dem/
if exist('~/mambaforge/envs/paramo/bin/gmt','file')
    gmtBin = '~/mambaforge/envs/paramo/bin/gmt';
elseif exist('~/mambaforge/bin/gmt','file')
    gmtBin = '~/mambaforge/bin/gmt';
elseif exist('~/miniforge3/bin/gmt','file')
    gmtBin = '~/miniforge3/bin/gmt';
else
    gmtBin = 'gmt';
end

%%
if minLon >= -92 && maxLon <= -88.62 && minLat >= -1.8 && maxLat <= 0.75
    % GALAPAGOS
    %load('~/igdata/dem/galapagos_dem.mat','lat','lon','demData'); %demData = demData'; 
    runCmd = strcat(gmtBin,' grdcut -R',num2str(minLon),'/',...
        num2str(maxLon),'/',num2str(minLat),'/',num2str(maxLat),...
        ' ~/igdata/dem/Galapagos1s.nc -Gec.cut.grd');
    unix(runCmd);
    try
        lon = ncread('~/igdata/dem/ec.cut.grd','lon');
    catch
        lon = ncread('~/igdata/dem/ec.cut.grd','x');
    end

    try
        lat = ncread('~/igdata/dem/ec.cut.grd','lat');
    catch
        lat = ncread('~/igdata/dem/ec.cut.grd','y');
    end
    demData = ncread('~/igdata/dem/ec.cut.grd','/z');
    demData = double(demData);
elseif minLon >= -78 && maxLon <= -77 && minLat >= 0 && maxLat <= 1
    % CHILES
    load('~/igdata/dem/chiles_dem_highRes.mat','lat','lon','demData');
    lonI = lon >= minLon & lon <= maxLon;
    latI = lat >= minLat & lat <= maxLat;
    lon = lon(lonI);
    lat = lat(latI);
    demData = double(demData(latI,lonI))';
elseif minLon <= -84 && minLat <= -5 && maxLat >= 1
    runCmd = strcat(gmtBin,' grdcut -R',num2str(minLon),'/',...
        num2str(maxLon),'/',num2str(minLat),'/',num2str(maxLat),...
        ' ~/igdata/dem/topo_downsample.grd -Gec.cut.grd');
    disp(runCmd);
    unix(runCmd);

    try
        lon = ncread('~/igdata/dem/ec.cut.grd','lon');
    catch
        lon = ncread('~/igdata/dem/ec.cut.grd','x');
    end

    try
        lat = ncread('~/igdata/dem/ec.cut.grd','lat');
    catch
        lat = ncread('~/igdata/dem/ec.cut.grd','y');
    end
    demData = ncread('~/igdata/dem/ec.cut.grd','/z');
    demData = double(demData);
else
    runCmd = strcat(gmtBin,' grdcut -R',num2str(minLon),'/',...
        num2str(maxLon),'/',num2str(minLat),'/',num2str(maxLat),...
        ' ~/igdata/dem/Ecuador1s.grd -Gec.cut.grd');
    disp(runCmd);
    unix(runCmd);

    try
        lon = ncread('~/igdata/dem/ec.cut.grd','lon');
    catch
        lon = ncread('~/igdata/dem/ec.cut.grd','x');
    end

    try
        lat = ncread('~/igdata/dem/ec.cut.grd','lat');
    catch
        lat = ncread('~/igdata/dem/ec.cut.grd','y');
    end
    demData = ncread('~/igdata/dem/ec.cut.grd','/z');
    demData = double(demData);
end

%%
numelData = numel(demData);
while numelData > 80e6
    decimationFactor = 2;
    fprintf('numel: %d; decimating the data in each direction by: %d...\n',...
        numelData,decimationFactor);
    lonq = lon(1:decimationFactor:end);
    latq = lat(1:decimationFactor:end);
    Vq = interp2(lon,lat,demData',lonq',latq);
    demData = Vq';
    lon = lonq;
    lat = latq;
    numelData = numel(demData);
end
fprintf('done decimating, number of elements is: %d\n',numelData); 
cd(PWD);
