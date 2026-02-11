clear; close all;
dayStart = datetime(2025,12,18);
dayEnd = datetime(2026,01,26);

dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

thresh = 0.07;
basis_function_file_name = "cotopaxi";
basis_function_directory = fullfile("~","masa","subspace_detector","basis_functions");
basisFunctions = load(fullfile(basis_function_directory,basis_function_file_name));

writeFlag = true;
maxEvents = 1e4;
maxBasisFunctions = 20;
diffFlag = false;
linearccnorm = true;
customPrefix = "cotopaxi";

coto_kstnms = ["CO1V";"BREF";"BVC2";"BTAM"];
coto_channels = ["HHZ";"HHN";"HHE";...
    "BHZ";"BHN";"BHE"];

for i = 1:lDays
    tic;
    day_ = dayVec(i);
    if day_ >= datetime(2025,01,23)
        coto_kstnms = ["CO1V";"BVC2";"BTAM"];
    else
        coto_kstnms = ["CO1V";"BREF";"BVC2";"BTAM"];
    end

    S = loadWaveforms(day_,dayInc,coto_kstnms,coto_channels,...
        "EC","");
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
customPrefix = "cotopaxi"; %'CotopaxiMultiplexed_v1';
cd(fullfile("~","masa","subspace_detector",customPrefix));
files = dir('*.mat');
lFiles = length(files);
maxN = 5e6;
tabs = NaT(maxN,1);
energyRatio = NaN(maxN,1);
Neff = energyRatio;
z2p = NaN(maxN,4);
z2pSynth = z2p;
p2rms = z2p;
kurt = z2p;
skew = z2p;
ReconstructCoefficients = NaN(80,maxN);
obsAmpRatio = NaN(6,maxN);
synthAmpRatio = NaN(6,maxN);

omitCO1V = ~true;
nTot = 1;
for i = 1:lFiles
    name_ = files(i).name;
    fprintf("%d %s\n",i,name_);
    S = load(name_);

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
    obsAmpRatio(:,nTot:nTot+ll-1) = obsAmpRatio_;
    synthAmpRatio(:,nTot:nTot+ll-1) = synthAmpRatio_;

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
obsAmpRatio(:,nTot:end) = [];
synthAmpRatio(:,nTot:end) = [];
clear *_

%%
Neff = sum(isfinite(z2p),2);
medZ2P = median(z2p,2,'omitnan');
maxZ2P = max(z2p,[],2,'omitnan');
minZ2P = min(z2p,[],2,'omitnan');
meanZ2P = mean(z2p,2,'omitnan');
maxZ2P_synth = max(z2pSynth,[],2,'omitnan');
minZ2P_synth = min(z2pSynth,[],2,'omitnan');
meanZ2P_synth = mean(z2pSynth,2,'omitnan');
medZ2P_synth = median(z2pSynth,2,'omitnan');

maxKurt = max(kurt,[],2,'omitnan');
minKurt = min(kurt,[],2,'omitnan');
maxP2rms = max(p2rms,[],2,'omitnan');
minP2rms = min(p2rms,[],2,'omitnan');
maxSkew = max(skew,[],2,'omitnan');
minSkew = min(skew,[],2,'omitnan');
toc;

%%
tic;
MINTHRESH = 0.2;
MINSTATIONS = 2;
minmaxAmp = 4e2;
MINMED = 5e2;
MINMAX = 5e2;
P2RMS_THRESH = 8;
MINKURT = 10;
MINMEAN = 5e2;
tStart = datetime(2006,01,01);

goodI = maxZ2P_synth >= MINMAX & ...
    tabs >= tStart & ...
    energyRatio >= MINTHRESH & ...
    Neff >= MINSTATIONS & ...
    medZ2P_synth >= MINMED & ...
    maxP2rms >= P2RMS_THRESH & ...
    maxKurt >= MINKURT & ...
    meanZ2P_synth >= MINMEAN;

close all;
nDays = 7;
rate = t2r(tabs(goodI),days(nDays),[],true);
rate = rate/nDays;
axL3 = linkedPlot(tabs(goodI),...
    [maxZ2P_synth(goodI),...
    maxKurt(goodI),...
    maxP2rms(goodI),...
    energyRatio(goodI),...
    meanZ2P_synth(goodI),...
    rate]);

for i = 1:5
    axL3(i).YScale = "log";
end

%%
nDays = 30;
rate = t2r(tabs(goodI),days(nDays),[],true)/nDays;
axL4 = linkedPlot(tabs(goodI),[maxZ2P(goodI),rate]);
axL4(1).YScale = "log";

tGood = tabs(goodI);
ampGood = maxZ2P(goodI);

[tsheet,type,Tsp,amp,magnitude,codaDuration,period,refStation] = loadVolcanicEvents("cotopaxi",7);
[tsheet2,amp2,type2,Tsp2,magnitude2,codaDuration2,period2,refStation2] = ...
    filterCatalog(tsheet,amp,10,type,Tsp,magnitude,codaDuration,period,refStation);

%[tsheet,amp,type,Tsp,magnitude,codaDuration,period,refStation] = filterCatalog(tsheet,amp,5,type,Tsp,magnitude,codaDuration,period,refStation);
excelI = tsheet >= tStart & amp >= minmaxAmp & amp <= 1e7; % & contains(type,"LP");
tsheet(~excelI) = [];
type(~excelI) = [];
Tsp(~excelI) = [];
amp(~excelI) = [];
magnitude(~excelI) = [];
codaDuration(~excelI) = [];
period(~excelI) = [];
refStation(~excelI) = [];

[t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(tGood,tsheet,30,true);
length(just_t1I)
length(tGood)

figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(211);
semilogy(tGood(just_t1I),ampGood(just_t1I),'o'); zoom on; grid on; hold on; semilogy(tsheet,amp,'.'); title('all excel, unique subspace')
ax(2) = subplot(212);
semilogy(tGood,ampGood,'o'); zoom on; grid on; hold on; semilogy(tsheet(just_t2I),amp(just_t2I),'.'); title('all subspace, unique excel');
linkaxes(ax,'xy');

%%
tUniqExcel = tsheet(just_t2I);
tUniqExcel2 = tUniqExcel - dateshift(tUniqExcel,'start','day');

figure('units','normalized','outerposition',[0 0 1 1]);
plot(sort(hours(tUniqExcel2)),(0:length(tUniqExcel2)-1)'/length(tUniqExcel2),'.'); zoom on; grid on;
nDays = 1; rateExcel = t2r(tUniqExcel,days(nDays),[],true)/nDays;
nDays = 30; rateExcel30 = t2r(tUniqExcel,days(nDays),[],true)/nDays;
linkedPlot(tUniqExcel,[rateExcel,rateExcel30]);

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tGood(just_t1I),ampGood(just_t1I),'o'); zoom on; grid on; hold on;
semilogy(tsheet(just_t2I),amp(just_t2I),'.'); legend('unique subspace','unique excel');

tUniqSubspace = tGood(just_t1I);
nDays = 1; rateUniqSubspace = t2r(tUniqSubspace,days(nDays),[],true)/nDays;
nDays = 30; rateUniqSubspace30 = t2r(tUniqSubspace,days(nDays),[],true)/nDays;
linkedPlot(tUniqSubspace,[rateUniqSubspace,rateUniqSubspace30]);

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tGood(t1_commonI),ampGood(t1_commonI),'o'); zoom on; grid on; hold on;
semilogy(tsheet(t2_commonI),amp(t2_commonI),'.'); legend('common subspace','common excel');
toc;