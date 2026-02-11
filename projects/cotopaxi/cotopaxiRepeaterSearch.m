function cotopaxiRepeaterSearch(dayStart,dayEnd,dayInc)
%clear; close all; %clc;

%% version 8
% basisFunctionFileName = '~/research/now/cotopaxi/cotopaxi_svd_basis_functions_v5.mat';
% 10 april 2021
threshold = 0.17;
recordLength = 35;
maxN = 1e4;
mpd = 10;
maxTemplates = 15;
linearccnorm = true;
plotFlag = false;
verboseFlag = false;
diffFlag = false;
saveFlag = true;

%% version 7
% threshold = 0.85;
% maxTemplates = 15;
% mpd = 15;
% recordLength = 20;
% maxN = 1e3;
% linearccnorm = true;
% plotFlag = false;
% verboseFlag = false;
% diffFlag = false;
% saveFlag = true;

%saveHome = "~/research/now/cotopaxi/";
saveHome = "~/igdata";
templateFileName = fullfile(saveHome,"cotopaxi_svd_basis_functions_v5.mat");
temporaryCatalogFile = fullfile(saveHome,"cotoTmpUpdate_v4")';
longTermCatalogFile = fullfile(saveHome,"cotopaxiSubspaceDetectorBREF_v8");
%longTermCatalogFile = '~/research/now/cotopaxi/OCT2022_cotopaxiSubspaceDetectorBREF_v8';
load(templateFileName,'kstnm','chan');

%%
% dayStart = datetime(2024,05,01);
% dayEnd = datetime(2024,05,13);
%dayInc = 1; % shouldnt be bigger than 4, else run into memory issues
%lk = size(cotoSNCLs,1); % version 7
lk = 18; %size(kstnm); % version 8

dayRange = (dayStart:dayInc:dayEnd)';
lDays = length(dayRange);

%%
for i = 1:lDays
    %% read data
    day_ = dayRange(i);

    %Sorig = populateWaveforms(lk);
    %nn = 0;
    %for j = 1:lk
    %         SNCL_ = cotoSNCLs(j,:);
    %         kstnm_ = SNCL_(1);
    %         kcmpnm_ = SNCL_(2);
    %         knetwk_ = SNCL_(3);
    %         khole_ = SNCL_(4);
    %         S_ = loadWaveforms(days(i),dInc,kstnm_,kcmpnm_,knetwk_,khole_,true,verboseFlag);

    S_ = loadWaveforms(day_,dayInc,unique(kstnm),unique(chan)); % use this one for version 5 of basis functions
    refs = pull(S_,'ref');
    badref = isnat(refs);
    Sorig = S_(~badref);

    %if ~isnat(S_.ref)
    %    nn = nn + 1;
    %    Sorig(nn,1) = S_;
    %end
    %end

    clear S_ %SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    %Sorig = Sorig(1:nn,1);
    lS = length(Sorig);

    if ~lS
        fprintf(2,'no data for day: %s\n',datestr(day_));
        fprintf('\n');
        continue;
    end

    %% gen new
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
    save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff','-v7.3');

    %%
    clear tabs p2rms z2p NCC kurt Neff
end
