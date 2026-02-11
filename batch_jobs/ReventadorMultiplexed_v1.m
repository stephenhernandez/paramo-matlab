clear; close all;
%dayStart = datetime(2018,04,24);
%dayEnd = datetime(2018,04,24);
dayStart = datetime(2022,03,20);
dayEnd = datetime(2022,03,20);
%dayStart = datetime(2023,09,13);
%dayEnd = datetime(2023,09,13);

dayInc = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

thresh = 0.1;
basisFunctionFileName = '~/igdata/ReventadorMultiplexedBasisFunctions_v1.mat';
reve_kstnms = ...%["REVN";"REVS";
    ["CAYR";"CASC";"BONI";"YAHU";"IMBA";"ANTS";...
    "ANTG";"CUSW";"URCU";"COTA";"CUIC";"TULM";"PULU";"GGPC";"GGPT"];

writeFlag = ~true;
maxEvents = 2e3;
maxBasisFunctions = 20;
linearccnorm = true;
diffFlag = true;
customPrefix = 'ReventadorMultiplexed_v1';
basisFunctions = load(basisFunctionFileName);
reve_channels = ["HHZ";"HHN";"HHE";...
    "BHZ";"BHN";"BHE";...
    "SHZ";"SHN";"SHE"];

%%
for i = 1:lDays
    tic;
    day_ = dayVec(i);
    [S,nSuccess] = loadWaveforms(day_,dayInc,...
        reve_kstnms,...
        reve_channels,...
        "EC","",~true,~true);
    if ~nSuccess
        S = loadWaveforms(day_,dayInc,...
            reve_kstnms,...
            reve_channels,...
            "EC","",~true,~true,'~/rawdata_old');
    end
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
basisFunctionFileName = '~/masa/subspace_detector/reventador/ReventadorMultiplexedBasisFunctions_v1.mat';
U = load(basisFunctionFileName,'U');
U = U.U;
[maxBasisFunctions,maxStations] = size(U);
clear U;

%customPrefix = 'ReventadorMultiplexed_v1';
customPrefix = "reventador";
cd(fullfile("~","masa","subspace_detector",customPrefix));
files = dir('daySubspace*.mat');
lFiles = length(files);
maxN = 1e6;
tabs = NaT(maxN,1);
energyRatio = NaN(maxN,1);
Neff = energyRatio;
z2p = NaN(maxN,maxStations);
z2pSynth = z2p;
p2rms = z2p;
kurt = z2p;
skew = z2p;
ReconstructCoefficients = NaN(maxStations*maxBasisFunctions,maxN);
maxAmpRatios = 0.5*maxStations*(maxStations-1);
obsAmpRatio = NaN(maxAmpRatios,maxN);
synthAmpRatio = NaN(maxAmpRatios,maxN);
nTot = 1;

sensitivities = [1256780000;... % REVN
    314195000;...   % REVS
    1;...           % CAYR
    503735000;...   % CASC
    503735000;...   % BONI
    314195000;...   % YAHU
    314195000;...   % IMBA
    503735000;...   % ANTS
    503735000;...   % ANTG
    1109410000;...  % CUSW
    315873000;...   % URCU
    554705000;...   % COTA 1168440000
    314195000;...   % CUIC
    314195000;...   % TULM
    314195000;...   % PULU 1256780000
    314195000;...   % GGPC
    314195000];     % GGPT

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
        z2p(nTot:nTot+ll-1,:) = z2p_;%./sensitivities';
        z2pSynth(nTot:nTot+ll-1,:) = z2pSynth_;%./sensitivities';
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
energyRatio(nTot:end) = [];
z2p(nTot:end,:) = [];
z2pSynth(nTot:end,:) = [];

% %boni correction
% tI = tabs < datetime(2020,12,03);
% z2p(tI,5) = z2p(tI,5)/4;
% z2pSynth(tI,5) = z2pSynth(tI,5)/4;
%
% %pulu correction
% tI = tabs < datetime(2019,12,16);
% z2p(tI,15) = z2p(tI,15)/4;
% z2pSynth(tI,15) = z2pSynth(tI,15)/4;
% clear tI;

p2rms(nTot:end,:) = [];
kurt(nTot:end,:) = [];
skew(nTot:end,:) = [];
ReconstructCoefficients(:,nTot:end) = [];
obsAmpRatio(:,nTot:end) = [];
synthAmpRatio(:,nTot:end) = [];

%
[tabs,energyRatio,z2p,z2pSynth,p2rms,kurt,skew,ReconstructCoefficients,obsAmpRatio,synthAmpRatio] = ...
    filterCatalog(tabs,energyRatio,60,z2p,z2pSynth,p2rms,kurt,skew,ReconstructCoefficients',obsAmpRatio',synthAmpRatio');
ReconstructCoefficients = ReconstructCoefficients';
obsAmpRatio = obsAmpRatio';
synthAmpRatio = synthAmpRatio';

clear *_
z2p(z2p <= 10) = NaN;
%z2p(z2p <= 1/314195000) = NaN;
maxZ2P = max(z2p,[],2,'omitnan');
Neff = sum(isfinite(z2p),2);
Neff_skipREV = sum(isfinite(z2p(:,4:end)),2);
Neff_afterCASC = sum(isfinite(z2p(:,5:end)),2);
Neff_afterCASC_synth = sum(isfinite(z2pSynth(:,5:end)),2);
CASCamp = z2p(:,4);
CASCamp_synth = z2pSynth(:,4);
medZ2P = median(z2p,2,'omitnan');
minZ2P = min(z2p,[],2,'omitnan');
maxKurt = max(kurt,[],2,'omitnan');
minKurt = min(kurt,[],2,'omitnan');
maxP2rms = max(p2rms,[],2,'omitnan');
minP2rms = min(p2rms,[],2,'omitnan');
maxSkew = max(skew,[],2,'omitnan');
minSkew = min(skew,[],2,'omitnan');
toc;

load ~/igdata/ecuadorSensorDataTable10.mat
reve_kstnms = ["REVN";"REVS";"CAYR";"CASC";"BONI";"YAHU";"IMBA";"ANTS";...
    "ANTG";"CUSW";"URCU";"COTA";"CUIC";"TULM";"PULU";"GGPC";"GGPT"];
[stla,stlo] = metaDataFromStationList(reve_kstnms);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(-0.080506,-77.657999,stla,stlo,refEllipse)*1e-3;
dReg = d_(4:end);

%%
% tic;
% close all;
% minStations = 2;
% minAmp = 1e1;
% minmaxAmp = 3e2;
% goodI = tabs >= datetime(2015,01,01) & Neff == 15 & ...
%     Neff >= minStations & energyRatio >= 0.2 ...
%     & minZ2P >= minAmp & medZ2P >= minmaxAmp & maxKurt < 60 & minP2rms < 7 & minSkew > -1; % & maxSkew > 0 & maxSkew < 1 & maxP2rms < 25;

tic;
close all;
minStations = 2;
MINTHRESH = 0.15;
%minAmp = 1e0/314195000; %1e1
minAmp = 1e1; %
%minmaxAmp = 1e-6; %7e1/314195000; %1e2;
minmaxAmp = 7e2;

% goodI = tabs >= datetime(2021,11,10) & Neff_afterCASC >= minStations & ...
%     Neff_afterCASC_synth >= minStations & CASCamp_synth >= 300 & ...
%     CASCamp >= 300 & energyRatio >= 0.2;% & isfinite(CASCamp);

goodI = tabs >= datetime(2016,01,01) & Neff >= minStations & energyRatio >= MINTHRESH ...
    & minZ2P >= minAmp & medZ2P >= minmaxAmp & (CASCamp_synth >= 200 | ~isfinite(CASCamp_synth)); % & maxKurt < 60 & minP2rms < 7 & minSkew > -1; % & maxSkew > 0 & maxSkew < 1 & maxP2rms < 25;

% goodI = tabs >= datetime(2010,01,01) & Neff_skipREV >= 14 & energyRatio >= MINTHRESH ...
%     & minZ2P >= minAmp & medZ2P >= minmaxAmp; % & maxKurt < 60 & minP2rms < 7 & minSkew > -1; % & maxSkew > 0 & maxSkew < 1 & maxP2rms < 25;
% fprintf("%d\n",sum(goodI));

% goodI = tabs >= datetime(2012,01,01) & Neff >= minStations & energyRatio >= MINTHRESH ...
%     & minZ2P >= minAmp & medZ2P >= minmaxAmp & isfinite(z2p(:,4));

% goodI = tabs >= datetime(2018,05,09) & tabs <= datetime(2018,05,11) & ...
%     Neff >= minStations & energyRatio >= 0.15 ...
%     & minZ2P >= minAmp & medZ2P >= minmaxAmp; % & maxP2rms < 10;

close all;
axL3 = linkedPlot(tabs(goodI),...
    [medZ2P(goodI),maxKurt(goodI),minP2rms(goodI),energyRatio(goodI),minSkew(goodI)],...
    ".","compact");
axL3(1).YScale = "log";
%axL3 = linkedPlot(tabs(goodI),...
%    [z2p(goodI,4),maxKurt(goodI),minP2rms(goodI),energyRatio(goodI),minSkew(goodI)]);

%%
figure();
plot(tabs(goodI),Neff(goodI),'.'); zoom on; grid on;
nDays = 1;
rateCaus = t2r(tabs(goodI),days(nDays),[],true)/nDays;
rateAcaus = t2r(tabs(goodI),days(nDays),[],false)/nDays;

axL4 = linkedPlot(tabs(goodI),[medZ2P(goodI),rateAcaus./rateCaus]);
axL4(1).YScale = 'log';
axL4(2).YScale = 'log';

axL5 = linkedPlot(tabs(goodI),[medZ2P(goodI),rateCaus]);
axL5(1).YScale = 'log';
axL5(1).YLabel.String = "amplitude proxy";
axL5(2).YLabel.String = "daily rate";

axL6 = linkedPlot(tabs(goodI),[CASCamp_synth(goodI),CASCamp(goodI)]);
axL6(1).YScale = 'log';
axL6(2).YScale = 'log';

%%
X = [energyRatio(goodI) log10([z2pSynth(goodI,5:end) z2p(goodI,5:end)])];
y = CASCamp_synth(goodI);
badI = ~isfinite(y);
X(badI,:) = [];
y(badI) = [];
[b,stats] = robustfit(X,log10(y));
ynew = 10.^([ones(size(X,1),1) X]*b);
t = tabs(goodI & isfinite(CASCamp_synth));
figure(); semilogy(t,[y ynew],'.'); zoom on; grid on;
h = stats.h;
length(t)

%%
% Nmed = 51;
% close all;
% for i = 5:17
%     figure();
%     semilogy(tabs(goodI),z2pSynth(goodI,4)./z2pSynth(goodI,i),'.');
%     zoom on; grid on;
%     hold on;
%     semilogy(tabs(goodI),medfiltSH(z2pSynth(goodI,4)./z2pSynth(goodI,i),Nmed,true),'.');
%     title(sprintf("%s-%s","CASC",reve_kstnms(i)));
% end
%
% close all;
% for i = 3:17
%     figure();
%     semilogy(tabs(goodI),z2pSynth(goodI,4)./z2pSynth(goodI,i),'.');
%     zoom on; grid on;
%     hold on;
%     semilogy(tabs(goodI),medfiltSH(z2pSynth(goodI,4)./z2pSynth(goodI,i),Nmed,true),'.');
%     title(sprintf("%s-%s","REVS",reve_kstnms(i)));
% end