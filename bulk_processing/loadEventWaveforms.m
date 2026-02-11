function [Sall,E,goodI] = loadEventWaveforms(fnames,varargin)
%
% Sall = loadEventWaveforms(fnames,startTime,endTime,refFlag,threeFlag,saveFlag,prsFlag,diasFlag,rawDataDir)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Saturday, Jul 27, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    -10,...                     % startTime
    60,...                      % endTime
    1,...                       % refFlag
    1,...                       % threeFlag
    0,...                       % save flag
    1,...                       % print record section?
    0};                         % diasFlag

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[startTime,endTime,refFlag,threeFlag,saveFlag,prsFlag,diasFlag] = deal(optsToUse{:});

%%
lfiles = length(fnames);
E = populateSeisCompStructure(lfiles);
stnms = [];
ntwks = [];
locids = [];
ptimes = [];
ids = [];
torig = [];
dists = [];
sncl = [];

%% read all sc3 files
for i = 1:lfiles
    E(i) = readSCBulletin(fnames(i),diasFlag);
    t = E(i).t;
    id_ = string(E(i).id);

    Pphases = E(i).Pphases;
    nPphases = E(i).nPphases;

    stnm_ = pull(Pphases,'stnm');
    ntwk_ = pull(Pphases,'ntwk');
    chan_ = char(pull(Pphases,'chan'));
    chan_ = chan_(:,1:2);
    locid_ = pull(Pphases,'locid');
    pt_ = pull(Pphases,'t');
    dist_ = pull(Pphases,'dist');

    stnms = [stnms; stnm_];                                         %#ok<AGROW>
    ntwks = [ntwks; ntwk_];                                         %#ok<AGROW>
    locids = [locids; locid_];                                      %#ok<AGROW>
    ptimes = [ptimes; pt_];                                         %#ok<AGROW>
    ids = [ids; repmat(id_,nPphases,1)];                            %#ok<AGROW>
    torig = [torig; repmat(t,nPphases,1)];                          %#ok<AGROW>
    dists = [dists; dist_];                                         %#ok<AGROW>
    sncl = [sncl; ...
        strcat(stnm_,".",ntwk_,".",string(chan_),".",locid_)];      %#ok<AGROW>
end

%%
refEllipse = referenceEllipsoid('wgs84');
idMaster = pull(E,'id');
eqlat = pull(E,'lat');
eqlon = pull(E,'lon');
eqdepth = pull(E,'depth');
eqmag = pull(E,'mag');

%%
uniqSNCL = unique(sncl);
splitStr = regexp(uniqSNCL,'\.','split');
uniqSNCLSplit = [];
for ii = 1:length(uniqSNCL)
    str_ = splitStr{ii};
    uniqSNCLSplit = [uniqSNCLSplit; str_];                          %#ok<AGROW>
end
[stla,stlo,stel] = metaDataFromStationList(uniqSNCLSplit(:,1),uniqSNCLSplit(:,2),...
    strcat(uniqSNCLSplit(:,3),"Z"),uniqSNCLSplit(:,4));

isfin = isfinite(stla);
stla = stla(isfin);
stlo = stlo(isfin);
stel = stel(isfin);
lsncls = sum(isfin);

%%
Sall = populateWaveforms([lfiles,lsncls]);
badSensors = uniqSNCL(~isfin);
uniqSNCL = uniqSNCL(isfin);
uniqSNCLSplit = uniqSNCLSplit(isfin,:);

for i = 1:length(badSensors)
    bs = sncl == badSensors(i);
    stnms(bs) = [];
    ntwks(bs) = [];
    locids(bs) = [];
    ptimes(bs) = [];
    ids(bs) = [];
    torig(bs) = [];
    dists(bs) = [];
    sncl(bs) = [];
end

if refFlag
    allTimes = NaT(lsncls,lfiles);
else
    allTimes = NaN(lsncls,lfiles);
end

if threeFlag
    Sall = repmat(Sall,1,3);
    allTimes = repmat(allTimes,3,1);
end

%%
for j = 1:lsncls
    sncl_ = uniqSNCL(j);
    lia = ismember(sncl,sncl_);
    idTmp = ids(lia);
    id_lia = ismember(idMaster,idTmp);
    sumEvents = sum(id_lia);

    %%
    stnmTmp = uniqSNCLSplit(j,1);
    chanTmp = char(uniqSNCLSplit(j,3));
    locidTmp = uniqSNCLSplit(j,4);

    %%
    [delKM,azs] = distance(eqlat(id_lia),eqlon(id_lia),stla(j),stlo(j),refEllipse);
    [~,bazs] = distance(stla(j),stlo(j),eqlat(id_lia),eqlon(id_lia),refEllipse);

    %%
    if refFlag
        refTime = ptimes(lia);
    else
        refTime = torig(lia);
    end

    if ismac
        if strcmp(stnmTmp,"VCH1") || strcmp(stnmTmp,"PVIL") || ...
                strcmp(stnmTmp,"FER1") || strcmp(stnmTmp,"FER2") || ...
                strcmp(stnmTmp,"CEAZ") || strcmp(stnmTmp,"ALCE")
            ntwkTmp = "9D";
        else
            ntwkTmp = uniqSNCLSplit(j,2);
        end
    else
        ntwkTmp = uniqSNCLSplit(j,2);
    end

    %%
    tic;
    if threeFlag
        [S_,successFlag] = extractWaveforms(refTime+seconds(startTime),endTime,stnmTmp,...
            [string([chanTmp(1:2),'Z']); string([chanTmp(1:2),'N']); string([chanTmp(1:2),'E'])],...
            ntwkTmp,locidTmp,false,false);
        if successFlag
            for k = 1:3
                S_(:,k) = push(S_(:,k),'evid',idTmp);
                S_(:,k) = push(S_(:,k),'evla',eqlat(id_lia));
                S_(:,k) = push(S_(:,k),'evlo',eqlon(id_lia));
                S_(:,k) = push(S_(:,k),'evdp',eqdepth(id_lia));
                S_(:,k) = push(S_(:,k),'eqmag',eqmag(id_lia));

                S_(:,k) = push(S_(:,k),'stla',repmat(stla(j),sumEvents,1));
                S_(:,k) = push(S_(:,k),'stlo',repmat(stlo(j),sumEvents,1));
                S_(:,k) = push(S_(:,k),'stel',repmat(stel(j),sumEvents,1));
                S_(:,k) = push(S_(:,k),'dist',delKM*1e-3);
                S_(:,k) = push(S_(:,k),'az',azs);
                S_(:,k) = push(S_(:,k),'baz',bazs);
                Sall(id_lia,k+3*(j-1)) = S_(:,k);
                if refFlag
                    allTimes(k+3*(j-1),id_lia) = refTime;
                else
                    allTimes(k+3*(j-1),id_lia) = dists(lia);
                end
            end
        end
    else
        [S_,successFlag] = extractWaveforms(refTime+seconds(startTime),endTime,stnmTmp,...
            string([chanTmp(1:2),'Z']),ntwkTmp,locidTmp,false,false);
        if successFlag
            S_ = push(S_,'evid',idTmp);
            S_ = push(S_,'evla',eqlat(id_lia));
            S_ = push(S_,'evlo',eqlon(id_lia));
            S_ = push(S_,'evdp',eqdepth(id_lia));
            S_ = push(S_,'eqmag',eqmag(id_lia));
            S_ = push(S_,'stla',repmat(stla(j),sumEvents,1));
            S_ = push(S_,'stlo',repmat(stlo(j),sumEvents,1));
            S_ = push(S_,'stel',repmat(stel(j),sumEvents,1));
            S_ = push(S_,'dist',delKM*1e-3);
            S_ = push(S_,'az',azs);
            S_ = push(S_,'baz',bazs);
            Sall(id_lia,j) = S_;
            if refFlag
                allTimes(j,id_lia) = refTime;
            else
                allTimes(j,id_lia) = dists(lia);
            end
        end
    end
    toc;
end
Sall = Sall';

%%
for i = 1:lfiles
    at = allTimes(:,i);
    [~,sI] = sort(at);
    Sall(:,i) = Sall(sI,i);
end

%%
sizeS = size(Sall);
refs = pull(Sall,'ref');
goodI = ~isnat(reshape(refs,sizeS));

%%
delI = sum(goodI,2) == 0;
goodI(delI,:) = [];
Sall(delI,:) = [];

%%
if saveFlag
    saveDir = "~/products/events/html/";
    if ~exist(saveDir,"dir")
        mkdir(saveDir)
    end
    cd(saveDir);
    for i = 1:lfiles
        idNow = char(idMaster(i));
        disp(idNow)
        if ~exist(idNow,'dir')
            mkdir(idNow);
        end
        cd(idNow);
        W = Sall(:,i);
        Wref = pull(W,'ref');
        isnatW = isnat(Wref);
        W(isnatW) = [];
        W = W(:);
        if length(W) == 1
            W = struct2table(W,'AsArray',true);
        else
            W = struct2table(W);
        end
        save(strcat(idNow,'_waveforms'),'W','-v7.3');
        cd ..;
    end
end

if prsFlag
    saveDir = "~/products/events/html/";
    if ~exist(saveDir,"dir")
        mkdir(saveDir)
    end
    cd(saveDir);
    for i = 1:lfiles
        idNow = char(idMaster(i));
        disp(idNow)
        if ~exist(idNow,'dir')
            mkdir(idNow);
        end
        cd(idNow);
        W = Sall(:,i);
        W = differentiateWaveforms(W);
        W = detrendWaveforms(W);
        W = taperWaveforms(W,0.02);
        W = filterWaveforms(W,0.25,2);
        W = intWaveforms(W);
        plotRecordSection(W,-inf,-inf,false,true,idNow,'Z',false);
        cd ..;
    end
end
