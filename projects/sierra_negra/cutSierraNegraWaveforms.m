% this code is mostly to test that loadEventWaveforms works
clear; close all; clc;

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Friday, Aug 2, 2019

%%
threeFlag = true;
refFlag = true;
startTime = -5;
endTime = 60;

%%
cd ~/research/now/sierra_negra/IGUANA_PICKS_20180624_20180627/BULLETINS
files = dir('*.bulletin');
files = files(1:10);
lfiles = length(files);

diasFlag = true;
E = populateSeisCompStructure(lfiles);
stnms = [];
ntwks = [];
chans = [];
locids = [];
ptimes = [];
ids = [];
torig = [];

for i = 1:lfiles
    E(i) = readSCBulletin(files(i).name,diasFlag);
    t = E(i).t;
    id_ = string(E(i).id);
    
    Pphases = E(i).Pphases;
    nPphases = E(i).nPphases;
    stnm_ = struct2data(Pphases,'stnm');
    ntwk_ = struct2data(Pphases,'ntwk');
    chan_ = struct2data(Pphases,'chan');
    locid_ = struct2data(Pphases,'locid');
    pt_ = struct2data(Pphases,'t');
    
    stnms = [stnms; stnm_];
    ntwks = [ntwks; ntwk_];
    chans = [chans; chan_];
    locids = [locids; locid_];
    ptimes = [ptimes; pt_];
    ids = [ids; repmat(id_,nPphases,1)];
    torig = [torig; repmat(t,nPphases,1)];
end
idMaster = struct2data(E,'id');

%%
uniqStnms = unique(stnms);
lstnm = length(uniqStnms);
Sall = populateWaveforms([lfiles,lstnm]);
if threeFlag
    Sall = repmat(Sall,1,3);
end

%%
for j = 1:lstnm
    stnmTmp = uniqStnms(j);
    lia = ismember(stnms,stnmTmp);
    locb = find(lia);
    chanTmp = char(chans(locb(1)));
    locidTmp = locids(locb(1));
    idTmp = ids(lia);
    
    if refFlag
        refTime = ptimes(lia);
    else
        refTime = torig(lia);
    end
    
    if strcmp(stnmTmp,"SN07")
        rawDataDir = '~/data/iguana/NO_GPS/';
    else
        rawDataDir = '~/data/iguana/BROADBAND/';
    end
    
    if strcmp(stnmTmp,"VCH1")
        ntwkTmp = "9D";
    else
        ntwkTmp = ntwks(locb(1));
    end
    
    if threeFlag
        tic;
        [S_,successFlag] = extractWaveforms(refTime+seconds(startTime),endTime,stnmTmp,...
            [string([chanTmp(1:2),'Z']); string([chanTmp(1:2),'N']); string([chanTmp(1:2),'E'])],...
            ntwkTmp,locidTmp,false,false,rawDataDir);
        toc;
        if successFlag
            id_lia = ismember(idMaster,idTmp);
            Sall(id_lia,1+3*(j-1)) = S_(:,1);
            Sall(id_lia,2+3*(j-1)) = S_(:,2);
            Sall(id_lia,3+3*(j-1)) = S_(:,3);
        end
    else
        tic;
        [S_,successFlag] = extractWaveforms(refTime+seconds(startTime),endTime,stnmTmp,...
            string([chanTmp(1:2),'Z']),ntwkTmp,locidTmp,false,false,rawDataDir);
        toc;
        if successFlag
            id_lia = ismember(idMaster,idTmp);
            Sall(id_lia,j) = S_;
        end
    end
end
