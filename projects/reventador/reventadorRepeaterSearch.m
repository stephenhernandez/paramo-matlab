function reventadorRepeaterSearch(dayStart,dayEnd,dayInc)
%clear; close all; clc;

% version 8, 29 december 2021
threshold = 0.3;
recordLength = 90;
maxN = 2e3;
mpd = 20;
maxTemplates = 10;
linearccnorm = true;
plotFlag = false;
verboseFlag = false;
diffFlag = false;
writeFlag = true;

% templateFileName = '~/research/now/reventador/CASC_basisFunctions_13DEC2021_prepared';
% temporaryCatalogFile = '~/research/now/reventador/revTmpUpdate_v3';
% longTermCatalogFile = '~/research/now/reventador/ReventadorSubspaceDetectorResults_v7';
% load(templateFileName,'kstnm','chan');

%dataHome = '~/igdata';
dataHome = '~/research/now/reventador';

templateFileName = fullfile(dataHome,'CASC_BONI_basisFunctions_29DEC2021');
temporaryCatalogFile = fullfile(dataHome,'revTmpUpdate_v4');
longTermCatalogFile = fullfile(dataHome,'ReventadorSubspaceDetectorResults_v10');
load(templateFileName,'kstnm','chan');

%%
% dayStart = datetime(2024,05,13);
% dayEnd = datetime(2024,05,13);
% dayInc = 1;
lk = size(kstnm,1); % version 8

dayRange = (dayStart:dayInc:dayEnd)';
lDays = length(dayRange);

basisFunctions = load(templateFileName);
newFs = basisFunctions.newFs;

%%
for i = 1:lDays
    day_ = dayRange(i);

    S_ = loadWaveforms(day_,dayInc,unique(kstnm),unique(chan)); % use this one for version 5 of basis functions
    refs = pull(S_,'ref');
    badref = isnat(refs);
    Sorig = S_(~badref);

    clear S_ %SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    %Sorig = Sorig(1:nn,1);
    lS = length(Sorig);

    if ~lS
        fprintf(2,'no data for day: %s\n',day_);
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

    if ~writeFlag
        disp('----------');
        disp('DATA NOT BEING SAVED TO FILE.');
        disp('----------');
        continue;
    end

    %%
    T = load(longTermCatalogFile,'tabs');
    lastT = min([min(day_) + seconds(60) max(T.tabs)]);
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
    newI = tNew > lastT;
    keepI = tabs <= lastT;

    %%
    tabs = [tabs(keepI); tNew(newI)];
    p2rms = [p2rms(keepI,:); prmsNew(newI,:)];
    z2p = [z2p(keepI,:); z2pNew(newI,:)];
    NCC = [NCC(keepI); NCCNew(newI)];
    kurt = [kurt(keepI,:); kurtNew(newI,:)];
    Neff = [Neff(keepI); neffNew(newI)];

    %%
    clear tNew prmsNew z2pNew NCCNew kurtNew neffNew lastT;

    %%
    save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff','-v7.3');

    %%
    clear tabs p2rms z2p NCC kurt Neff
end