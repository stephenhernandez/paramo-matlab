clear; close all;
%dayStart = datetime(2010,01,01);
%dayEnd = datetime(2024,06,17);
dayStart = datetime(2018,06,16);
dayEnd = datetime(2018,06,16);


dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

thresh = 0.07;
basisFunctionFileName = '~/igdata/payg_multiplexed_basis_functions.mat';
writeFlag = ~true;
maxEvents = 2e3;
maxBasisFunctions = 20;
diffFlag = true;
linearccnorm = true;
customPrefix = 'PaygMultiplexed_v1';
basisFunctions = load(basisFunctionFileName);

%%
for i = 1:lDays
    tic;
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        "PAYG",...
        ["BHZ";"BH1";"BH2"],...
        "IU","00");
    toc;

    badI = isnat(pull(S,'ref'));
    S(badI) = [];
    lS = length(S);
    if ~lS
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end

    try
        [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
            obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
            multiplexedSubspaceDetector(S,thresh,basisFunctions,...
            maxBasisFunctions,maxEvents,writeFlag,diffFlag,linearccnorm,...
            customPrefix);
        toc;
    catch
        fprintf('something went wrong on day: %s',day_);
        continue;
    end
end

%%
plotFlag = true;
if ~plotFlag
    return;
end

tic;
clear; close all; clc;
customPrefix = 'PaygMultiplexed_v1';
cd(fullfile('~','subspaceDetector',customPrefix));
files = dir('*.mat');
lFiles = length(files);
maxN = 5e6;
tabs = NaT(maxN,1);
energyRatio = NaN(maxN,1);
Neff = energyRatio;
z2p = NaN(maxN,1);
z2pSynth = z2p;
p2rms = z2p;
kurt = z2p;
skew = z2p;
ReconstructCoefficients = NaN(20,maxN);
% obsAmpRatio = NaN(1,maxN);
% synthAmpRatio = NaN(1,maxN);

omitCO1V = ~true;
nTot = 1;
for i = 1:lFiles
    disp(i);
    S = load(files(i).name);

    tabs_ = S.tabs;
    energyRatio_ = S.energyRatio;
    Neff_ = S.Neff;
    z2p_ = S.z2p;
    z2pSynth_ = S.z2pSynth;
    p2rms_ = S.p2rms;
    kurt_ = S.kurt;
    skew_ = S.skew;
    ReconstructCoefficients_ = S.ReconstructCoefficients;
    obsAmpRatio_ = S.obsAmpRatio;
    synthAmpRatio_ = S.synthAmpRatio;
    
    ll = length(tabs_);
    if omitCO1V
        z2p(nTot:nTot+ll-1,2:4) = z2p_(:,2:4); % = cat(1,z2p,z2p_);
        z2pSynth(nTot:nTot+ll-1,2:4) = z2pSynth_(:,2:4); % = cat(1,z2pSynth,z2pSynth_);
        p2rms(nTot:nTot+ll-1,2:4) = p2rms_(:,2:4); % = cat(1,p2rms,p2rms_);
        kurt(nTot:nTot+ll-1,2:4) = kurt_(:,2:4); % = cat(1,kurt,kurt_);
        skew(nTot:nTot+ll-1,2:4) = skew_(:,2:4); % = cat(1,skew,skew_);
    else
        z2p(nTot:nTot+ll-1,:) = z2p_;
        z2pSynth(nTot:nTot+ll-1,:) = z2pSynth_;
        p2rms(nTot:nTot+ll-1,:) = p2rms_;
        kurt(nTot:nTot+ll-1,:) = kurt_;
        skew(nTot:nTot+ll-1,:) = skew_;
    end

    tabs(nTot:nTot+ll-1) = tabs_;
    energyRatio(nTot:nTot+ll-1) = energyRatio_;
    Neff(nTot:nTot+ll-1) = Neff_;
    ReconstructCoefficients(:,nTot:nTot+ll-1) = ReconstructCoefficients_;
    %obsAmpRatio(:,nTot:nTot+ll-1) = obsAmpRatio_;
    %synthAmpRatio(:,nTot:nTot+ll-1) = synthAmpRatio_;

    nTot = nTot + ll;
end

%
tabs(nTot:end) = [];
Neff(nTot:end) = [];
energyRatio(nTot:end) = [];
z2p(nTot:end,:) = [];
z2pSynth(nTot:end,:) = [];
p2rms(nTot:end,:) = [];
kurt(nTot:end,:) = [];
skew(nTot:end,:) = [];
ReconstructCoefficients(:,nTot:end) = [];
clear *_

z2pI = tabs < datetime(2018,09,30);
z2p(z2pI) = z2p(z2pI)/1.68136541126086;

Neff = sum(isfinite(z2p),2);
medZ2P = median(z2p,2,'omitnan');
maxZ2P = max(z2p,[],2,'omitnan');
minZ2P = min(z2p,[],2,'omitnan');
maxKurt = max(kurt,[],2,'omitnan');
minKurt = min(kurt,[],2,'omitnan');
maxP2rms = max(p2rms,[],2,'omitnan');
minP2rms = min(p2rms,[],2,'omitnan');
maxSkew = max(skew,[],2,'omitnan');
minSkew = min(skew,[],2,'omitnan');
toc; 

%%
tic;
minStations = 1;
minmaxAmp = 1e3; %5e2;
minAmp = 2e2;
tStart = datetime(2010,01,01);
goodI = tabs >= tStart & ...
    Neff >= minStations & energyRatio >= 0.15 ...
    & minZ2P >= minAmp & maxKurt <= 4; % & maxZ2P < 4e5 & maxZ2P >= minmaxAmp & ...
    %maxKurt <= 4 & maxP2rms >= 5;

close all;
axL3 = linkedPlot(tabs(goodI),maxZ2P(goodI),maxKurt(goodI),maxP2rms(goodI),energyRatio(goodI),maxSkew(goodI));
for i = 1:4
    axL3(i).YScale = 'log';
end
toc;

%figure(); plot(tabs(goodI),Neff(goodI),'.'); zoom on; grid on;
nDays = 1;
rate = t2r(tabs(goodI),days(nDays),[],true)/nDays;
axL4 = linkedPlot(tabs(goodI),medZ2P(goodI),rate);
axL4(1).YScale = 'log';