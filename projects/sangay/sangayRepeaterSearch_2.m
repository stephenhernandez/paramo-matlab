%function reventadorRepeaterSearch()
clear; close all; %clc;

%% version 2, 14 january 2022
threshold = 0.25;
recordLength = 45;
maxN = 1e4;
mpd = 10;
maxTemplates = 10;
linearccnorm = true;
plotFlag = false;
verboseFlag = true;
diffFlag = false;
saveFlag = true;

templateFileName = '~/research/now/sangay/bdf_basis_functions_2';
temporaryCatalogFile = '~/research/now/sangay/sagaTmpUpdate_v21';
longTermCatalogFile = '~/research/now/sangay/SAGASubspaceDetectorResults_v21';
load(templateFileName,'kstnm','chan','ntwk','locID');

%%
dInc = 1;
lk = size(kstnm,1); % version 8


tStart = datetime(2022,01,31);
tEnd = dn2dt(ceil(datenum(dn2dt(now)+hours(5)))); %run on Feb 1st %datetime(2022,01,15); %datetime(2020,09,30);

dayRange = (tStart:dInc:tEnd)';
lDays = length(dayRange);

basisFunctions = load(templateFileName);
newFs = basisFunctions.newFs;

%%
for i = 1:lDays
    dayStart = dayRange(i);

    S_ = loadWaveforms(dayStart,dInc,unique(kstnm),unique(chan),"EC",locID); % use this one for version 5 of basis functions
    refs = pull(S_,'ref');
    badref = isnat(refs);
    Sorig = S_(~badref);

    clear S_ %SNCL_ kstnm_ kcmpnm_ knetwk_ khole_
    %Sorig = Sorig(1:nn,1);
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

    %%
    %Scut = extractWaveforms(tabs-seconds(0),seconds(recordLength),"SAGA","BDF","EC","01");
    %Scut = cutWaveforms(Sorig,tabs,0,seconds(recordLength));
    
    Scut = resampleWaveforms(intWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(Sorig)),basisFunctions.lfc,basisFunctions.hfc)),newFs);
    Sf2 = cutWaveforms(Scut,tabs,0,seconds(recordLength));

    d = double(pull(Sf2));
    winlen = 5*newFs; 
    Nmax = size(d,2); 
    snrCurve = fftfilt(ones(winlen,1)/winlen,abs(hilbert(d(:,1:Nmax))).^2);

    snrCurve = zpkFilter(abs(sqrt(snrCurve)),-inf,2/winlen,1,1,1);
    snrCurve2 = taper(snrCurve(winlen:end,:)./snrCurve(1:end-winlen+1,:),winlen/size(snrCurve,1));
    [snr,maxI] = max(snrCurve2); 
    snr = snr'; 
    maxI = maxI';

    z2p = 2.38e-6*z2p*1000/20;

    %%
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
    save(temporaryCatalogFile,'tabs','p2rms','z2p','NCC','Neff','kurt','snr');
    clear tabs p2rms templateIndex z2p NCC saveStds Neff kurt snr

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
    snrNew = load(temporaryCatalogFile,'snr');
    snrNew = snrNew.snr;

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
    snr = [snr(tKeep); snrNew(newGood)];

    %%
    clear tNew prmsNew z2pNew NCCNew kurtNew neffNew lastT snrNew;

    %%
    save(longTermCatalogFile,'tabs','p2rms','z2p','NCC','kurt','Neff','snr');

    %%
    clear tabs p2rms z2p NCC kurt Neff snr
end
