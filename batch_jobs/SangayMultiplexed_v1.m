clear; close all;
dayStart = datetime(2024,05,07);
dayEnd = datetime(2024,05,09);

dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

thresh = 0.07;
%basisFunctionFileName = '~/igdata/SangayMultiplexedBasisFunctions_v1.mat';
basisFunctionFileName = '~/igdata/SangayMultiplexedBasisFunctions_v2.mat';

sangay_kstnms = ["SAGA";"BPAT";"BMAS";"BRTU";"BULB";"BBIL";"BRUN";"PUYO";"TAMH";"PORT";"PKYU";"TAIS"];
%sangay_kstnms = ["BPAT";"BMAS";"BRTU";"BULB";"BBIL";"BRUN";"PUYO";"TAMH";"PORT";"PKYU";"TAIS"];

writeFlag = true;
maxEvents = 1e4;
maxBasisFunctions = 20;
linearccnorm = false;
diffFlag = true;
%customPrefix = 'SangayMultiplexed_v3_justSAGA';
customPrefix = 'SangayMultiplexed_v1';
basisFunctions = load(basisFunctionFileName);

%%
for i = 1:lDays
    tic;
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        sangay_kstnms,...
        ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"],...
        "EC",""); %,true,true);
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
tic;
plotFlag = true;
if ~plotFlag
    return;
end

clear; close all; clc;
customPrefix = "SangayMultiplexed_v1";
%customPrefix = 'SangayMultiplexed_v3_justSAGA';
cd(fullfile("~","masa","subspace_detector",customPrefix));
files = dir('*.mat');
lFiles = length(files);
maxN = 7e6;
tabs = NaT(maxN,1);
energyRatio = NaN(maxN,1);
Neff = energyRatio;
z2p = NaN(maxN,12);
z2pSynth = z2p;
p2rms = z2p;
kurt = z2p;
skew = z2p;
ReconstructCoefficients = NaN(240,maxN);
obsAmpRatio = NaN(66,maxN);
synthAmpRatio = NaN(66,maxN);
nTot = 1;

omitSAGA = ~true;
for i = 1:lFiles
    fName = files(i).name;
    disp(fName);

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
    if omitSAGA
        z2p(nTot:nTot+ll-1,2:end) = z2p_(:,2:end); % = cat(1,z2p,z2p_);
        z2pSynth(nTot:nTot+ll-1,2:end) = z2pSynth_(:,2:end); % = cat(1,z2pSynth,z2pSynth_);
        p2rms(nTot:nTot+ll-1,2:end) = p2rms_(:,2:end); % = cat(1,p2rms,p2rms_);
        kurt(nTot:nTot+ll-1,2:end) = kurt_(:,2:end); % = cat(1,kurt,kurt_);
        skew(nTot:nTot+ll-1,2:end) = skew_(:,2:end); % = cat(1,skew,skew_);
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
    obsAmpRatio(:,nTot:nTot+ll-1) = obsAmpRatio_;
    synthAmpRatio(:,nTot:nTot+ll-1) = synthAmpRatio_;
    nTot = nTot + ll;
end
toc; 

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
obsAmpRatio(:,nTot:end) = [];
synthAmpRatio(:,nTot:end) = [];

clear *_
z2p(z2p <= 10) = NaN;
maxZ2P = max(z2p,[],2,'omitnan');
% [tabs,maxZ2P,energyRatio,z2p,z2pSynth,p2rms,kurt,skew,...
%     ReconstructCoefficients,obsAmpRatio,synthAmpRatio] = ...
%     filterCatalog(tabs,maxZ2P,10,energyRatio,z2p,z2pSynth,p2rms,kurt,skew,...
%     ReconstructCoefficients',obsAmpRatio',synthAmpRatio');
% ReconstructCoefficients = ReconstructCoefficients';
% obsAmpRatio = obsAmpRatio';
% synthAmpRatio = synthAmpRatio';

Neff = sum(isfinite(z2p),2);
medZ2P = median(z2p,2,'omitnan');
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
minStations = 3;
minAmp = 1e1;
minmaxAmp = 5e1;
goodI = tabs >= datetime(2010,01,01) & ...
    Neff >= minStations & energyRatio >= 0.2 ...
    & minZ2P >= minAmp & maxZ2P < 1e7 & maxZ2P >= minmaxAmp; % & maxP2rms < 10;

close all;
axL3 = linkedPlot(tabs(goodI),[medZ2P(goodI),maxKurt(goodI),maxP2rms(goodI),energyRatio(goodI),maxSkew(goodI)],".","compact");
% figure(1);
% hold(axL3(1),'on');
% plot(axL3(1),tabs(goodI),minZ2P(goodI),'.'); zoom on; grid on;
for i = 1:4
    axL3(i).YScale = 'log';
end

%figure(); plot(tabs(goodI),Neff(goodI),'.'); zoom on; grid on;
nDays = 1/2;
rate = t2r(tabs(goodI),days(nDays),[],true)/nDays;
axL4 = linkedPlot(tabs(goodI),[medZ2P(goodI),rate],".","compact");
axL4(1).YScale = 'log';
% 
% axL5 = linkedPlot(tabs(goodI),minZ2P(goodI),rate);
% axL5(1).YScale = 'log';

toc;
