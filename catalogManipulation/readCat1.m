function [evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1(varargin)

% Inputs:
% tStart = start time (01-01-2000)
% tEnd = start end (01-01-2050)
% minMag = min. mag (-999)
% maxMag = max. mag (10)

%% parse inputs, deal variables
nVarargin = length(varargin);
functionDefaults = {datetime(2000,01,01),datetime(2050,01,01),-999,10,true};
optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,minMag,maxMag,keepSCBulletin] = deal(optsToUse{:});

%% read sc3 catalog and filter
catalogFile = fullfile("~","phaseInformationSC5","cat1.txt");
fileID = fopen(catalogFile);
frmt = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s %f %f %f %s %s %s %f %f %f %f %f %f %s %s';
C1 = textscan(fileID,frmt);
fclose(fileID);

%%
[y,m,d,h,mm,sec,timerr,...
    eqlat,eqlaterr,eqlon,eqlonerr,eqdepth,eqdeptherr,...
    stderr,azgap,phaseTot,...
    ids,locMethod,earthModel,evMode,evStatus,agencyID,...
    scbullmag,eqmagerr,nMLv,...
    magType,meanMethod,evDescription,...
    y2,m2,d2,hh2,mm2,sec2,authorID,evType] = deal(C1{:});

%%
magType = string(magType);
evMode = string(evMode);
evStatus = string(evStatus);
agencyID = string(agencyID);
locMethod = string(locMethod);
ids = string(ids);
earthModel = string(earthModel);
meanMethod = string(meanMethod);
evDescription = string(evDescription);
authorID = string(authorID);
evType = string(evType);

%%
t = datetime([y,m,d,h,mm,sec]);
creationTime = datetime([y2,m2,d2,hh2,mm2,sec2]);

%%
[t,sI] = sort(t);
creationTime = creationTime(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
scbullmag = scbullmag(sI);
eqmagerr = eqmagerr(sI);
timerr = timerr(sI);
magType = magType(sI);
phaseTot = phaseTot(sI);
nMLv = nMLv(sI);
stderr = stderr(sI);
azgap = azgap(sI);
evMode = evMode(sI);
evStatus = evStatus(sI);
agencyID = agencyID(sI);
evDescription = evDescription(sI);
ids = ids(sI);
locMethod = locMethod(sI);
earthModel = earthModel(sI);
meanMethod = meanMethod(sI);
authorID = authorID(sI);
evType = evType(sI);

%%
sI = t >= tStart & t < tEnd & ...
    ~(evStatus == "rejected") & ~(evType == "notlocatable");
t = t(sI);
creationTime = creationTime(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
scbullmag = scbullmag(sI);
eqmag = scbullmag;
eqmagerr = eqmagerr(sI);
timerr = timerr(sI);
magType = magType(sI);
phaseTot = phaseTot(sI);
nMLv = nMLv(sI);
stderr = stderr(sI);
azgap = azgap(sI);
evMode = evMode(sI);
evStatus = evStatus(sI);
agencyID = agencyID(sI);
evDescription = evDescription(sI);
ids = ids(sI);
locMethod = locMethod(sI);
earthModel = earthModel(sI);
meanMethod = meanMethod(sI);
authorID = authorID(sI);
evType = evType(sI);

%%
sI = eqmag >= minMag & eqmag <= maxMag;
eqmag = eqmag(sI);
t = t(sI);
creationTime = creationTime(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
eqmagerr = eqmagerr(sI);
timerr = timerr(sI);
magType = magType(sI);
phaseTot = phaseTot(sI);
nMLv = nMLv(sI);
stderr = stderr(sI);
azgap = azgap(sI);
evMode = evMode(sI);
evStatus = evStatus(sI);
agencyID = agencyID(sI);
evDescription = evDescription(sI);
ids = ids(sI);
locMethod = locMethod(sI);
earthModel = earthModel(sI);
meanMethod = meanMethod(sI);
scbullmag = scbullmag(sI);
authorID = authorID(sI);
evType = evType(sI);

%%
fprintf("running code to remove dups...\n");
difft = diff(t);
diffh = sqrt(diff(eqlat).^2 + diff(eqlon).^2);
dupi = seconds(difft) <= 1 & diffh < 1e-3;
sI = true(size(t)); %keep
sI(dupi) = false;

eqmag = eqmag(sI);
t = t(sI);
creationTime = creationTime(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
eqmagerr = eqmagerr(sI);
timerr = timerr(sI);
magType = magType(sI);
phaseTot = phaseTot(sI);
nMLv = nMLv(sI);
stderr = stderr(sI);
azgap = azgap(sI);
evMode = evMode(sI);
evStatus = evStatus(sI);
agencyID = agencyID(sI);
evDescription = evDescription(sI);
ids = ids(sI);
locMethod = locMethod(sI);
earthModel = earthModel(sI);
meanMethod = meanMethod(sI);
scbullmag = scbullmag(sI);
authorID = authorID(sI);
evType = evType(sI);

%%
if keepSCBulletin
    return;
end

%% keepSC5 instead...
[magMain,~,~,idMain] = keepSC5OrigMags();
eqmag = magMain; %<-- eqmag changed to SC5 magnitudes...
sI = ismember(ids,idMain);
t = t(sI);
creationTime = creationTime(sI);
eqlon = eqlon(sI);
eqlonerr = eqlonerr(sI);
eqlat = eqlat(sI);
eqlaterr = eqlaterr(sI);
eqdepth = eqdepth(sI);
eqdeptherr = eqdeptherr(sI);
scbullmag = scbullmag(sI);
eqmagerr = eqmagerr(sI);
timerr = timerr(sI);
magType = magType(sI);
phaseTot = phaseTot(sI);
nMLv = nMLv(sI);
stderr = stderr(sI);
azgap = azgap(sI);
evMode = evMode(sI);
evStatus = string(evStatus(sI));
agencyID = agencyID(sI);
evDescription = string(evDescription(sI));
ids = ids(sI);
locMethod = locMethod(sI);
earthModel = earthModel(sI);
meanMethod = meanMethod(sI);
authorID = authorID(sI);
evType = evType(sI);