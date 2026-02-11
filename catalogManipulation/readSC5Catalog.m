function [t,eqlat,eqlon,eqmag,eqdepth,rms,evDescription,azgap,id,eqlaterr,eqlonerr,...
    eqdeptherr,magType,used_phases,locMethod,earthModel] = readSC5Catalog(varargin)

% Inputs:
% regionName = 'ecuador';
% tStart = datetime(2010,01,01)
% tEnd = datetime(2030,01,01)
% typeLabel = 'ALL'

%% parse inputs, deal variables
nVarargin = length(varargin);
functionDefaults = {'ecuador',datetime(2010,01,01),datetime(2030,01,01),'ALL',-999};
optsToUse = functionDefaults;
if nVarargin
    optsToUse(1:nVarargin) = varargin;
end

[regionName,tStart,tEnd,typeLabel,minMag] = deal(optsToUse{:});
typeLabel = lower(typeLabel);

%% read sc3 catalog and filter
boundaryBox = getRegionSpatialDimensions(regionName);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);
minLon = boundaryBox(1);
maxLon = boundaryBox(2);

%%
[evDescription,~,t,eqlat,eqlon,eqdepth,~,id,rms,azgap,used_phases,~,...
    ~,eqlaterr,eqlonerr,eqdeptherr,~,magType,~,...
    locMethod,earthModel,~,~,~,eqmag] = readCat1(tStart);

%% sort by time
[t,sI] = sort(t);
locMethod = locMethod(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
evDescription = evDescription(sI);
id = id(sI);
earthModel = earthModel(sI);

%% sort by geography
sI = t >= tStart & t <= tEnd & eqmag >= minMag & eqlon >= minLon & eqlon <= maxLon & ...
    eqlat >= minLat & eqlat <= maxLat;

t = t(sI);
locMethod = locMethod(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
evDescription = evDescription(sI);
id = id(sI);
earthModel = earthModel(sI);

%% filter according to type label
vt = contains(evDescription,"vt",'ignorecase',true);
expl = contains(evDescription,"exp",'ignorecase',true);
vlp = contains(evDescription,"vlp",'ignorecase',true);
lp = contains(evDescription,"lp",'ignorecase',true);
hb = contains(evDescription,"hb",'ignorecase',true);
trem = contains(evDescription,"tre",'ignorecase',true) | contains(evDescription,"tra",'ignorecase',true) | contains(evDescription,"trm",'ignorecase',true);
reg = true(size(evDescription)) & ~(vt | expl | vlp | lp | hb | trem);

if strcmp(typeLabel,'lp')
    sI = lp;
elseif strcmp(typeLabel,'vt')
    sI = vt;
elseif strcmp(typeLabel,'exp')
    sI = expl;
elseif strcmp(typeLabel,'vlp')
    sI = vlp;
elseif strcmp(typeLabel,'hb')
    sI = hb;
elseif strcmp(typeLabel,'trem')
    sI = trem;
elseif strcmp(typeLabel,'regional') || strcmp(typeLabel,'ecuador') || strcmp(typeLabel,"reg")
    sI = reg;
elseif strcmp(typeLabel,'all')
    sI = true(size(t)); %everything!
    table(unique(evDescription),groupcounts(evDescription))
else
    fprintf(2,'innapropriate option: %s\n',typeLabel);
    fprintf(2,'must choose from: lp, vt, vlp, exp, hb, trem, regional, or all\n');
    return;
end

%%
t = t(sI);
locMethod = locMethod(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
evDescription = evDescription(sI);
id = string(id(sI));
earthModel = earthModel(sI);

%table(t,eqlon,eqlonerr,eqlat,eqlaterr,eqdepth,eqdeptherr,eqmag,magType,used_phases,rms,azgap,evDescription,id,locMethod,earthModel)