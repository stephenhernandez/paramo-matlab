function [E,t,lat,lon,depth,mag,...
    timerr,lonerr,laterr,deptherr,magerr,...
    nPphases,nSphases,rms,azgap,...
    nMLv,nML,nMjma,nMsBB,nMwp,nmB,nmb,ids,usedPhases] = extractSC3CatalogFromList(listName,textFlag)
if nargin < 2; textFlag = true; end

cd ~/phaseInformationSC3/
if textFlag
    list_ids = string(importdata(listName));
else
    list_ids = listName;
end

lids = length(list_ids);
E = populateSeisCompStructure(lids);
n = 0;
for i = 1:lids
    disp(i);
    id_ = char(list_ids(i));
    %ystr = id_(6:9);
    fname = id_; %fullfile(cd,ystr,strcat(id_,'.txt'));
    if exist(fname,'file')
        n = n+1;
        disp(fname)
        E(n) = readSCBulletin(fname);
        residualVector = pull(E(n).Pphases,'res');
        if E(n).nSphases
            residualVector = [residualVector; pull(E(n).Sphases,'res')]; %#ok<AGROW>
        end
        origrms_ = mad(residualVector,1);
        nmlv_ = E(n).nMLv;
        if nmlv_
            magStruct = E(n).MLv;
            mags_ = pull(magStruct,'value');
            E(n).magerr = std(mags_);
        end
        E(n).rms = origrms_;
    else
        disp('File does not exist')
    end
end
E = E(1:n);
clearvars -except E n

t           = pull(E,'t');
lat     = pull(E,'lat');
lon     = pull(E,'lon');
depth   = pull(E,'depth');
timerr  = pull(E,'timerr');
laterr  = pull(E,'laterr');
lonerr  = pull(E,'lonerr');
deptherr  = pull(E,'deptherr');
magerr  = pull(E,'magerr');
%mag  = pull(E,'mag');
usedPhases = pull(E,'usedPhases');
nPphases    = pull(E,'nPphases');
nSphases    = pull(E,'nSphases');
rms     = pull(E,'rms');
azgap     = pull(E,'azgap');
nMLv        = pull(E,'nMLv');
nML         = pull(E,'nML');
nMjma       = pull(E,'nMjma');
nMsBB       = pull(E,'nMsBB');
nMwp        = pull(E,'nMwp');
nmB         = pull(E,'nmB');
nmb         = pull(E,'nmb');
ids         = string(pull(E,'id'));

[t,tI] = sort(t);
E = E(tI);
lat = lat(tI);
lon = lon(tI);
depth = depth(tI);
timerr = timerr(tI);
laterr = laterr(tI);
lonerr = lonerr(tI);
deptherr = deptherr(tI);
magerr = magerr(tI);
usedPhases = usedPhases(tI);
nPphases = nPphases(tI);
nSphases = nSphases(tI);
rms = rms(tI);
azgap = azgap(tI);
nMLv = nMLv(tI);
nML = nML(tI);
nMjma = nMjma(tI);
nMsBB = nMsBB(tI);
nMwp = nMwp(tI);
nmB = nmB(tI);
nmb = nmb(tI);
ids = ids(tI);

[magMaster,idMaster] = keepSC3OrigMags();
tI = ismember(ids,idMaster);
if sum(~tI)
    disp('toss');
    disp(find(~tI));
end
E = E(tI);
t = t(tI);
lat = lat(tI);
lon = lon(tI);
depth = depth(tI);
timerr = timerr(tI);
laterr = laterr(tI);
lonerr = lonerr(tI);
deptherr = deptherr(tI);
magerr = magerr(tI);
usedPhases = usedPhases(tI);
nPphases = nPphases(tI);
nSphases = nSphases(tI);
rms = rms(tI);
azgap = azgap(tI);
nMLv = nMLv(tI);
nML = nML(tI);
nMjma = nMjma(tI);
nMsBB = nMsBB(tI);
nMwp = nMwp(tI);
nmB = nmB(tI);
nmb = nmb(tI);
ids = ids(tI);

%%
[~,tI] = ismember(ids,idMaster);
mag = magMaster(tI);

%%
for i = 1:length(ids)
    E(i).t =        t(i);
    E(i).lat =      lat(i);
    E(i).lon =      lon(i);
    E(i).depth =    depth(i);
    E(i).mag =      mag(i);
    E(i).timerr =   timerr(i);
    E(i).laterr =   laterr(i);
    E(i).lonerr =   lonerr(i);
    E(i).deptherr = deptherr(i);
    E(i).magerr =   magerr(i);
    E(i).usedPhases =   usedPhases(i);
    E(i).nPphases =     nPphases(i);
    E(i).nSphases =     nSphases(i);
    E(i).rms =      rms(i);
    E(i).azgap =      azgap(i);
    E(i).nMLv =         nMLv(i);
    E(i).nML =          nML(i);
    E(i).nMjma =        nMjma(i);
    E(i).nMsBB =        nMsBB(i);
    E(i).nMwp =         nMwp(i);
    E(i).nmB =          nmB(i);
    E(i).nmb =          nmb(i);
    E(i).id =          ids(i);
end