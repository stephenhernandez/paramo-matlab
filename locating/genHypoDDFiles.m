function [ax,E,E1,E2] = ...
    genHypoDDFiles(tStart,tEnd,regionName,boundaryBox,depth_correction,...
    diasFlag,minMag,maxDepth,min_total_phases,min_s_phases)
%[E,t,eqlat,eqlon,eqdepth,eqmag,id,...
%    tOrig,eqlatOrig,eqlonOrig,eqdepthOrig,eqmagOrig,idOrig] = ...
%    genHypoDDFiles(tStart,tEnd,regionName,boundaryBox,depth_correction,...
%    diasFlag,minMag,maxDepth,depthRange,minPhases,minSPhases)

if nargin < 1; tStart = datetime(2017,01,01); end
if nargin < 2; tEnd = datetime(2017,07,01); end
if nargin < 3; regionName = "ecuador"; end
if nargin < 4
    boundaryBox = getRegionSpatialDimensions(regionName);
end
if nargin < 5; depth_correction = 0; end
if nargin< 6; diasFlag = false; end
if nargin < 7; minMag = 2; end
if nargin < 8; maxDepth = 1000; end
if nargin < 9; min_total_phases = 3; end
if nargin < 10; min_s_phases = 0; end

%%
if depth_correction ~= 0
    T = readtable("estaciones.dat");
    stelv = T.Var4;
    stelv = stelv - (depth_correction*1000);
    T.Var4 = stelv;
    writetable(T,"estaciones_mod.dat","WriteVariableNames",false,"Delimiter"," ");
else
    !\cp -f estaciones.dat estaciones_mod.dat
end
magFact = 5;
[minLon,maxLon,minLat,maxLat] = deal(boundaryBox(1),boundaryBox(2),boundaryBox(3),boundaryBox(4));
hillShadeFlag = true;
depthPanel = true;
cumPanel = false;
minLevInc = 100;
minDepth = -depth_correction; %-6;

if isdatetime(tStart)
    SC5Flag = true;
    [~,~,t,eqlat,eqlon,eqdepth,~,ids,stderr,azgap,phaseTot,nMLv,...
        timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,~,...
        ~,~,~,~,~,eqmag,~,~] = readCat1(tStart,tEnd);
    eqdepthOrig = eqdepth;
    eqdepth = eqdepthOrig + depth_correction;

    totErr = sqrt(eqdeptherr.^2 + eqlonerr.^2 + eqlaterr.^2);
    depthConstraint = totErr./abs(eqdepth);

    maxMag = 8;
    errorDepthRatio1 = 10; %0.5;
    errorDepthRatio2 = 10;
    %minSPhases = 2; %5; %this criteria gets enforced further down in "sc2hypodd"
    minNmlv = 1; %7;

    magerrThresh = 2;
    azGapThresh = 355;
    spatialErrThresh = 20;
    timeErrThresh = 3;
    rmsThresh = 3;

    minDeptherr = -999; %0; %0.1; %0;

    goodI = t >= tStart & t <= tEnd & eqlon >= minLon & eqlon <= maxLon & ...
        eqlat >= minLat & eqlat <= maxLat & nMLv >= minNmlv & phaseTot >= min_total_phases & ...
        eqmag <= maxMag & eqmag >= minMag & azgap <= azGapThresh & ...
        abs(eqmagerr) <= magerrThresh & strcmpi(magType,"mlv") & ...
        timerr <= timeErrThresh & stderr <= rmsThresh & ...
        eqdepthOrig >= minDepth & eqdepthOrig <= maxDepth & ...
        eqdeptherr >= minDeptherr & eqlonerr >= minDeptherr & ...
        eqlaterr >= minDeptherr & (totErr < spatialErrThresh | depthConstraint <= errorDepthRatio1) & ...
        depthConstraint <= errorDepthRatio2;% & ...
    %~(strcmp(evDescription,"lp") | contains(evDescription,"vlp") | strcmp(evDescription,"hb") | contains(evDescription,"exp"));

    eqlatOrig = eqlat(goodI);
    eqlonOrig = eqlon(goodI);
    eqdepthOrig = eqdepthOrig(goodI);
    idOrig = ids(goodI);
    phaseTotOrig = phaseTot(goodI);
    nMLvOrig = nMLv(goodI);
    disp([min(eqdepthOrig) max(eqdepthOrig) sum(goodI)])

    %%
    lid = length(idOrig);
    E = populateSeisCompStructure(lid);
    tic;

    for i = 1:lid
        thisID = idOrig(i);
        fprintf("read file number %d/%d: %s\n",i,lid,thisID);
        E_ = readSCBulletin(thisID,diasFlag);
        lat_ = E_.lat;
        lon_ = E_.lon;
        eqdepth_ = E_.depth;
        eqlatOrig(i) = lat_;
        eqlonOrig(i) = lon_;
        eqdepthOrig(i) = eqdepth_;
        E_.lat = lat_;
        E_.lon = lon_;
        E_.depth = eqdepth_;
        E(i) = E_;
    end
    toc;

    horErr = sqrt(pull(E,'laterr').^2 + pull(E,'lonerr').^2);
    hI = horErr <= 50;
    E = E(hI);
    % eqlatOrig = eqlatOrig(hI);
    % eqlonOrig = eqlonOrig(hI);
    % eqdepthOrig = eqdepthOrig(hI);
    % phaseTotOrig = phaseTotOrig(hI);
    % nMLvOrig = nMLvOrig(hI);
    tOrig = pull(E,'t');
else
    E = tStart;
    tOrig = pull(E,'t');
    tStart = dateshift(min(tOrig),"start","day");
end
eqmagOrig = pull(E,'mag');
idOrig = pull(E,'id');
eqlatOrig = pull(E,'lat');
eqlonOrig = pull(E,'lon');
eqdepthOrig = pull(E,'depth');
nMLvOrig = pull(E,'nMLv');
nP = pull(E,'nPphases');
nS = pull(E,'nSphases');
phaseTotOrig = nP + nS;

%%
if strcmp(regionName,"chilesExpanded")
    regionName = "chiles";
end
baseDir = fullfile("~","masa","relocation");
WDIR = fullfile(baseDir,regionName,"hypodd");
PHASE_FILE = sprintf("phase_%s.dat",regionName);
outfile = fullfile(WDIR,PHASE_FILE);
cd(WDIR);

%%
sc2hypodd(E,outfile,depth_correction,min_total_phases,min_s_phases);
phase_input_file = sprintf("ph2dt_%s.inp",regionName); % input parameters, no data
phase_data_file = PHASE_FILE; %actual PHASE _data_!

%%
cd(WDIR);
unix('rm reloc_cat.txt');
unix('rm orig_cat.txt');
unix('rm no_reloc_cat.txt');
unix('rm hypoDD.reloc*');

%%
fprintf("phase input file: %s\n",phase_input_file);
try
    unix(sprintf("~/soft/hypodd/src/ph2dt/ph2dt %s",phase_input_file));
catch
    ax = [];
    E1 = [];
    E2 = E1;
    return;
end
unix('grep -v "\*\*\*\*" event.sel > tmp.txt');
unix("mv tmp.txt event.sel");
hypoStatus = unix("~/soft/hypodd/src/hypoDD/hypoDD hypoDD.inp");
if hypoStatus
    fprintf("hypoDD failed, check your code.\n");
    return;
end

%%
unix(sprintf("bash ~/scripts/separateCatalogs.sh %s",phase_data_file));
unix("wc reloc_cat.txt");

[t,lat,lon,depth,mag] = readHypoDD("reloc_cat.txt",true);
T = table(t,lat,lon,depth,mag);
E1 = table2struct(T);
[t,lat,lon,depth,mag] = readHypoDD("orig_catalog_events_also_relocated.txt",true);
T = table(t,lat,lon,depth,mag);
E2 = table2struct(T);
ax = compare_catalogs(E1,E2,[],depth_correction);