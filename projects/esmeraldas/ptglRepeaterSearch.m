%function esmeraldasRepeaterSearch()
clear; close all; %clc;

% version 1, 30 march 2022
threshold = 0.14;
recordLength = 32;
maxN = 2e3;
mpd = 15;
maxTemplates = 10;
linearccnorm = true;
plotFlag = true;
verboseFlag = true;
diffFlag = false;
saveFlag = false;

templateFileName = '~/research/now/esmeraldas/ptgl_event_2022/ptgl_2022_basisFunctions';
temporaryCatalogFile = '~/research/now/esmeraldas/esmeTmpUpdate_v1';
longTermCatalogFile = '~/research/now/esmeraldas/esmeraldasSubspaceDetectorResults_v1';
load(templateFileName,'kstnm','chan');

%%
dInc = 5;
lk = size(kstnm,1); % version 8

tStart = datetime(2022,03,26);
tEnd = datetime(2022,03,30);

dayRange = (tStart:dInc:tEnd)';
lDays = length(dayRange);

basisFunctions = load(templateFileName);
newFs = basisFunctions.newFs;

%%
for i = 1%:lDays
    dayStart = dayRange(i);

    S_ = loadWaveforms(dayStart,dInc,unique(kstnm),unique(chan)); % use this one for version 5 of basis functions
    refs = pull(S_,'ref');
    badref = isnat(refs);
    Sorig = S_(~badref);

    clear S_ %SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    lS = length(Sorig);

    if ~lS
        fprintf(2,'no data for day: %s\n',datestr(dayStart));
        fprintf('\n');
        continue;
    end

    %% gen new
    Sorig = resampleWaveforms(Sorig,newFs);
    [~,tabs,NCC,z2p,Neff,p2rms,kurt,ccnorm,t] = subspaceDetector(Sorig,...
        threshold,templateFileName,maxTemplates,recordLength,maxN,mpd,...
        diffFlag,linearccnorm,plotFlag,verboseFlag);

    %%
    if ~sum(~isnat(tabs))
        disp('NO EVENTS FOUND, NOT SAVING ANY DATA');
        continue;
    end

    if ~saveFlag
        disp('----------');
        disp('DATA NOT BEING SAVED TO FILE.');
        disp('----------');
        continue;
    end

    %%
    T = load(longTermCatalogFile,'tabs');
    lastT = min([min(pull(Sorig,'ref')) + seconds(60) max(T.tabs)]);
    clear T;

    %%
    save(temporaryCatalogFile,'tabs','p2rms','z2p','NCC','Neff','kurt');
    clear tabs p2rms templateIndex z2p NCC saveStds Neff kurt

    %% sync datasets
    load(longTermCatalogFile);

    %%
    tNew = load(temporaryCatalogFile,'tabs');
    tNew = tNew.tabs;
    prmsNew = load(temporaryCatalogFile,'p2rms');
    prmsNew = prmsNew.p2rms;
    z2pNew = load(temporaryCatalogFile,'z2p');
    z2pNew = z2pNew.z2p;
    NCCNew = load(temporaryCatalogFile,'NCC');
    NCCNew = NCCNew.NCC;
    kurtNew = load(temporaryCatalogFile,'kurt');
    kurtNew = kurtNew.kurt;
    neffNew = load(temporaryCatalogFile,'Neff');
    neffNew = neffNew.Neff;

    %%
    newGood = tNew > lastT;
    tKeep = tabs <= lastT;

    %%
    tabs = [tabs(tKeep); tNew(newGood)];
    p2rms = [p2rms(tKeep,:); prmsNew(newGood,:)];
    z2p = [z2p(tKeep,:); z2pNew(newGood,:)];
    NCC = [NCC(tKeep); NCCNew(newGood)];
    kurt = [kurt(tKeep,:); kurtNew(newGood,:)];
    Neff = [Neff(tKeep); neffNew(newGood)];

    %%
    clear tNew prmsNew z2pNew NCCNew kurtNew neffNew lastT;

    %%
    save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff');

    %%
    clear tabs p2rms z2p NCC kurt Neff
end
