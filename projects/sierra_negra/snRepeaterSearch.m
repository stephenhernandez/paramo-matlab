%function snRepeaterSearch()
clear; close all;

%% version 1
threshold = 0.2;
recordLength = 10;
maxN = 1e4;
mpd = 5;
maxTemplates = 20;
linearccnorm = true;
plotFlag = false;
verboseFlag = true;
smoothFlag = false;
diffFlag = false;
saveFlag = true;

templateFileName = '~/research/now/sierra_negra/sierra_negra_vch1_svd_basis_functions_noSNR_v1';
temporaryCatalogFile = '~/research/now/sierra_negra/snTmpUpdate_v1';
longTermCatalogFile = '~/research/now/sierra_negra/snSubspaceDetectorVCH1_v1';
load(templateFileName,'kstnm','chan');

%%
%dInc = 3;
%tStart = datetime(2018,06,25);

%dInc = 2;
%tStart = datetime(2016,04,16);

%dInc = 2;
%tStart = datetime(2020,01,12);

%dInc = 3;
%tStart = datetime(2017,03,20);

%dInc = 1;
%tStart = datetime(2014,09,01);

%
dInc = 1;
tStart = datetime(2012,06,11);
tEnd = datetime(2021,04,24); %datetime(2020,09,30);
dayRange = (tStart:dInc:tEnd)';
lDays = length(dayRange);

%%
for i = 1:lDays
    %% read data
    dayStart = dayRange(i);
    
    S_ = loadWaveforms(dayStart,dInc,unique(kstnm),unique(chan)); % use this one for version 5 of basis functions
    refs = pull(S_,'ref');
    badref = isnat(refs);
    Sorig = S_(~badref);
    if sum(badref) == 1
        S_ = loadWaveforms(dayStart,dInc,unique(kstnm),["BHZ";"BHN";"BHE"]); % use this one for version 5 of basis functions
        refs = pull(S_,'ref');
        badref = isnat(refs);
        if sum(badref) == 1
            fprintf(2,'no data for day: %s\n',datestr(dayStart));
        fprintf('\n');
            continue;
        end
        Sorig = S_(~badref);
        for kk = 1:sum(~badref)
            kcmpnm__ = char(Sorig(kk).kcmpnm);
            kcmpnm__(1) = 'H';
            Sorig(kk).kcmpnm = string(kcmpnm__);
        end
    end
    
    %%   
    clear S_ %SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    lS = length(Sorig);
    
    if ~lS
        fprintf(2,'no data for day: %s\n',datestr(dayStart));
        fprintf('\n');
        continue;
    end
    
    %%
    Sorig = interpolateWaveforms(Sorig);
    Sorig = resampleWaveforms(detrendWaveforms(Sorig),100);
    
    %% gen new
    [~,tabs,NCC,z2p,Neff,p2rms,kurt,ccnorm,t] = subspaceDetector(Sorig,...
        threshold,templateFileName,maxTemplates,recordLength,maxN,mpd,...
        diffFlag,linearccnorm,plotFlag,verboseFlag,smoothFlag);
    
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
    newGood = tNew > lastT & nanmedian(prmsNew,2) >= 4;
    tKeep = tabs <= lastT & nanmedian(p2rms,2) >= 4;
    fprintf('sum of newgood: %d\n',sum(newGood));
    
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
