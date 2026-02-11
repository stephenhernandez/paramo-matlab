%function sangayRepeaterSearch()
clear; close all; %clc;

sangaySNCLs = ["PUYO" "HHZ" "EC" "";...
    "TAIS" "HHZ" "EC" "";...
    "PKYU" "HHZ" "EC" "";...
    "TAMH" "HHZ" "EC" "";...
    "BPAT" "BHZ" "EC" "";...
    "BMAS" "BHZ" "EC" "";...
    "BULB" "BHZ" "EC" "";...
    "BRUN" "BHZ" "EC" "";...
    "PORT" "HHZ" "EC" ""];

dInc = 2; % shouldnt be bigger than 4, else run into memory issues
lk = size(sangaySNCLs,1);
%tStart = datetime(2019,05,04);
%tEnd = datetime(2019,05,07); %datetime(2020,09,30);
%tStart = datetime(2020,01,50);
%tEnd = datetime(2020,01,53); %datetime(2020,09,30);
%tStart = datetime(2020,09,27);
%tEnd = datetime(2020,09,30); %datetime(2020,09,30);

%tStart = datetime(2006,07,17);
%tEnd = datetime(2020,10,05); %datetime(2020,09,30);

%tStart = datetime(2020,09,19);
%tEnd = datetime(2020,10,05); %datetime(2020,09,30);

tStart = datetime(2016,04,16);
tEnd = datetime(2020,10,05); %datetime(2020,09,30);
days = (tStart:dInc:tEnd)';
lDays = length(days);

%%
threshold = 0.1;
maxTemplates = 5;
mpd = 30;
recordLength = 150;
maxN = 1e3;
linearccnorm = true;
plotFlag = true;
verboseFlag = true;
diffFlag = false;
saveFlag = 0;

%%
templateFileName = '~/research/now/sangay/sangay_svd_basis_functions_withSNR_v8';
temporaryCatalogFile = '~/research/now/sangay/sangayTenSensorsUpdate_v3';
longTermCatalogFile = '~/research/now/sangay/SangayRegionalAnalysis_v6.mat';

%%
for i = 1%:lDays
    %% read data
    Sorig = populateWaveforms(lk);
    nn = 0;
    for j = 1:lk
        SNCL_ = sangaySNCLs(j,:);
        kstnm_ = SNCL_(1);
        kcmpnm_ = SNCL_(2);
        knetwk_ = SNCL_(3);
        khole_ = SNCL_(4);
        
        S_ = loadWaveforms(days(i),dInc,kstnm_,kcmpnm_,knetwk_,khole_,true,verboseFlag);
        if ~isnat(S_.ref)
            nn = nn + 1;
            Sorig(nn,1) = S_;
        end
    end
    disp('done loading data');
    clear S_ SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    Sorig = Sorig(1:nn,1);
    lS = length(Sorig);
    
    %% gen new
    [~,tabs,NCC,z2p,Neff,p2rms,kurt,ccnorm,t] = subspaceDetector(Sorig,...
        threshold,templateFileName,maxTemplates,recordLength,maxN,mpd,...
        diffFlag,linearccnorm,plotFlag,verboseFlag);
    
    %%
    if isempty(tabs)
        disp('NO EVENTS FOUND, NOT SAVING ANY DATA');
        continue;
    end
    if ~saveFlag
        disp('----------');
        disp('DATA NOT BEING SAVED TO FILE.');
        disp('----------');
        continue;
    end
    
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
