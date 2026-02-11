function [t,eqlat,eqlon,eqmag,eqdepth,rms,type,azgap,id,eqlaterr,eqlonerr,...
    deptherr,magType,used_phases,velocityModel] = readSC3Catalog(varargin)

% Inputs:
% regionName = 'ecuador';
% sy = startYear (2010)
% sm = startMonth (01)
% sd = startDay (01)
% ey = endYear (2020)
% em = endMonth (01)
% ed = endDay (01)
% typeLabel = 'ALL'

%% parse inputs, deal variables
nVarargin = length(varargin);
functionDefaults = {'ecuador',2010,01,01,2030,01,01,'ALL'};
optsToUse = functionDefaults;
if nVarargin
optsToUse(1:nVarargin) = varargin;
end
[regionName,sy,sm,sd,ey,em,ed,typeLabel] = deal(optsToUse{:});

tbeg = datetime(sy,sm,sd);
tend = datetime(ey,em,ed);

%% read sc3 catalog and filter
boundaryBox = getRegionSpatialDimensions(regionName);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);
minLon = boundaryBox(1);
maxLon = boundaryBox(2);

%%
[type,~,t,eqlat,eqlon,eqdepth,eqmag,id,rms,azgap,used_phases,~,...
    ~,eqlaterr,eqlonerr,deptherr,~,magType,~,...
    ~,velocityModel] = readCat1();

[t,sI] = sort(t);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
deptherr = deptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
type = type(sI);
id = id(sI);
velocityModel = velocityModel(sI);

%%
sI = t >= tbeg & t <= tend & eqmag >= -30 & eqlon >= minLon & eqlon <= maxLon & ...
    eqlat >= minLat & eqlat <= maxLat;

t = t(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI); 
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
deptherr = deptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
type = type(sI); 
id = id(sI);
velocityModel = velocityModel(sI);

%% filter according to type label
if strcmp(typeLabel,'LP')
    sI = contains(type,"lp",'ignorecase',true);
elseif strcmp(typeLabel,'VT')
    sI = contains(type,"vt",'ignorecase',true);
elseif strcmp(typeLabel,'EXP')
    sI = contains(type,"exp",'ignorecase',true);
elseif strcmp(typeLabel,'VLP')
    sI = contains(type,"vlp",'ignorecase',true);
elseif strcmp(typeLabel,'HB')
    sI = contains(type,"hb",'ignorecase',true);
elseif strcmp(typeLabel,'TREM')
    sI = contains(type,"tre",'ignorecase',true);
elseif strcmp(typeLabel,'Regional') || strcmp(typeLabel,'Ecuador')
    sI = contains(type,"tect",'ignorecase',true) | ...
        contains(type,"unk",'ignorecase',true) | contains(type,"ecu",'ignorecase',true);
elseif strcmp(typeLabel,'ALL')
    sI = true(size(t)); %everything!
    disp(unique(type))
end

%%
t = t(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
deptherr = deptherr(sI);
eqmag = eqmag(sI);
magType = magType(sI);
used_phases = used_phases(sI);
rms = rms(sI);
azgap = azgap(sI);
type = type(sI); 
id = string(id(sI));
velocityModel = velocityModel(sI);
