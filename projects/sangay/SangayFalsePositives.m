clear; close all;

%
tic;
[tSAGA,aSAGA,nccSAGA] = filterUniqueEvents('~/research/now/sangay/sangaySubspaceDetectorSAGA_v2.mat');
[tSAGA,sI] = sort(tSAGA);
aSAGA = aSAGA(sI);
nccSAGA = nccSAGA(sI);
aSAGA(tSAGA < datetime(2018,11,11)) = aSAGA(tSAGA < datetime(2018,11,11))*4;

%
[tREG,aREG,nccREG,NeffReg] = filterUniqueEvents('~/research/now/sangay/SangayRegionalAnalysis_v6');
aREG = nanmedian(aREG,2);
Msangay = sangayMagnitudeCalculation(tREG,aREG);

%
[tREG,sI] = sort(tREG);
nccREG = nccREG(sI);
aREG = aREG(sI);
Msangay = Msangay(sI);

%
sI = tREG > tSAGA(1) & tREG < datetime(2021,02,01); %tSAGA(end);
%sI = tREG > tSAGA(1) & tREG < tSAGA(end); % & Msangay <= 3.15; %tSAGA(end);
tREG = tREG(sI);
aREG = aREG(sI);
nccREG = nccREG(sI);
NeffReg = NeffReg(sI);
Msangay = Msangay(sI);

%
tAll = [tSAGA; tREG];
iAll = [zeros(size(tSAGA)); ones(size(tREG))];
aAll = [aSAGA; aREG];
nccAll = [nccSAGA; nccREG];
neffAll = [zeros(length(tSAGA),1); NeffReg];
magAll = [zeros(length(tSAGA),1); Msangay];

%
[tAll,sI] = sort(tAll);
iAll = iAll(sI);
aAll = aAll(sI);
nccAll = nccAll(sI);
neffAll = neffAll(sI);
magAll = magAll(sI);
clear sI;

%
difft = seconds(diff(tAll));
t1 = tAll(1:end-1);
t2 = tAll(2:end);

regI = iAll(1:end-1);
a1 = aAll(1:end-1);     % SAGA amplitude
a2 = aAll(2:end);       % regional amplitude

ncc1 = nccAll(1:end-1);
ncc2 = nccAll(2:end);

neff1 = neffAll(1:end-1);
neff2 = neffAll(2:end);

mag1 = magAll(1:end-1);
mag2 = magAll(2:end);

%
% where i define SAGA true positives
tpSAGA = false(size(tAll));
tpREG = tpSAGA;
possibleTP = diff(iAll) == 1 & difft < 38 & difft >= 30 & a1 >= 1e4 & a1 <= 1e5 & ncc2 >= 0.1 & mag2 < 3 & neff2 >= 1 & a2./a1 <= 1e-2;
idum = find(possibleTP);
%idum = find(diff(iAll) == 1 & difft < 45 & difft >= 28 & a1 >= 1e4 & a1 <= 1e5 & ncc2 >= 0.1 & mag2 < 3 & neff2 >= 1);
tpSAGA(idum) = true;
tpREG(idum+1) = true;
clear idum;

%%
tRegTP = tAll(tpREG);
nccRegTP = nccAll(tpREG); %important, will serve as one of my most valuable features/predictors!
neffRegTP = neffAll(tpREG);

tSAGATP = tAll(tpSAGA);
nccSAGATP = nccAll(tpSAGA);
neffSAGATP = neffAll(tpSAGA);

%
tOrphans = tAll;
iOrphans = logical(iAll);
nccOrphans = nccAll;
neffOrphans = neffAll;

tOrphans(tpSAGA | tpREG) = [];
iOrphans(tpSAGA | tpREG) = [];
nccOrphans(tpSAGA | tpREG) = [];
neffOrphans(tpSAGA | tpREG) = [];

%
tRegOrphans = tOrphans(iOrphans);
nccRegOrphans = nccOrphans(iOrphans);
neffRegOrphans = neffOrphans(iOrphans);

tSAGAOrphans = tOrphans(~iOrphans);
nccSAGAOrphans = nccOrphans(~iOrphans);
neffSAGAOrphans = neffOrphans(~iOrphans);

%
lRegOrphans = length(tRegOrphans);
lSAGAOrphans = length(tSAGAOrphans)-1;
diffSAGAOrphans = diff(tSAGAOrphans);

%
tRegPossibleFP = NaT(lRegOrphans,1);
nccRegPossibleFP = NaN(lRegOrphans,1);
neffRegPossibleFP = NaN(lRegOrphans,1);

%
toc;
nn = 1;
thresh = seconds(40);
for i = 1:lSAGAOrphans
    tI = tRegOrphans < tSAGAOrphans(i+1);
    tRegPossibleFP_  = tRegOrphans(tI);
    nccRegPossibleFP_ = nccRegOrphans(tI);
    neffRegPossibleFP_ = neffRegOrphans(tI);
    
    %
    tRegPossibleFP_ = tRegPossibleFP_ - tSAGAOrphans(i);
    fpI = tRegPossibleFP_ > thresh & tRegPossibleFP_ < (diffSAGAOrphans(i) - thresh);
    nPossibleFP_thisBin = sum(fpI);
    
    %
    if nPossibleFP_thisBin %possible FPs found...
        tRegPossibleFP(nn:nn+nPossibleFP_thisBin-1) = tSAGAOrphans(i) + tRegPossibleFP_(fpI);
        nccRegPossibleFP(nn:nn+nPossibleFP_thisBin-1) = nccRegPossibleFP_(fpI);
        neffRegPossibleFP(nn:nn+nPossibleFP_thisBin-1) = neffRegPossibleFP_(fpI);
        nn = nn + nPossibleFP_thisBin;
        toc;
    end
    
    %
    tRegOrphans(tI) = [];       %clobber events already processed
    nccRegOrphans(tI) = [];     %clobber events already processed
    neffRegOrphans(tI) = [];
end

%
tRegPossibleFP = tRegPossibleFP(1:nn-1);
nccRegPossibleFP = nccRegPossibleFP(1:nn-1);
neffRegPossibleFP = neffRegPossibleFP(1:nn-1);

%
load('~/products/rsam/EC.SAGA..HHZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat','S');
d = S.d;
tRsam = getTimeVec(S);
gsi = find(diff(d<10 | ~isfinite(d)) == 1) + 1;
gei = find(diff(d<10 | ~isfinite(d)) == -1)+1;
gapDur = gei-gsi;
gapInfo = [gsi gapDur];

tStart = S.ref;
tEnd = S.ref + S.e;

keepI = true(size(tRegPossibleFP));
keepI(tRegPossibleFP < tStart | tRegPossibleFP > tEnd) = false;
%keepI(nccRegPossibleFP > 0.5) = false; %commented out on 07 march 2021

%
gI = gapDur >= 1;

%
gapStart = gapInfo(gI,1);
gapEnd = sum(gapInfo(gI,:),2);
lGaps = length(gapStart);
lSAGA_RSAM = S.npts;

%
disp('')
for i = 1:lGaps
    gapStart_ = gapStart(i);
    gapEnd_ = gapEnd(i);
    
    %
    if gapStart_ < lSAGA_RSAM && gapEnd_ <= lSAGA_RSAM
        tsI = tRegPossibleFP >= tRsam(gapStart_) & tRegPossibleFP < tRsam(gapEnd_);
        if sum(tsI)
            keepI(tsI) = false; %throw these
        end
    end
    toc;
end
tRegPossibleFP = tRegPossibleFP(keepI);
nccRegPossibleFP = nccRegPossibleFP(keepI);
neffRegPossibleFP = neffRegPossibleFP(keepI);

%
figure('units','normalized','outerposition',[0 0 1 1]);
ax__(1) = subplot(211);
semilogy(tRegPossibleFP,nccRegPossibleFP,'o'); zoom on; grid on;
hold on; semilogy(tRegTP,nccRegTP,'.');

ax__(2) = subplot(212);
plot(tRegPossibleFP,neffRegPossibleFP,'o'); zoom on; grid on;
hold on; plot(tRegTP,neffRegTP,'.');
linkaxes(ax__,'x');

%
timeMaster = [tRegTP; tRegPossibleFP];
labelMaster = [true(size(tRegTP)); false(size(tRegPossibleFP))];
nccMaster = [nccRegTP; nccRegPossibleFP];
neffMaster = [neffRegTP; neffRegPossibleFP];

istrain = false(size(timeMaster));
istrain([randsample(length(tRegTP),floor(length(tRegTP)/2)); ...
    length(tRegTP)+randsample(length(tRegPossibleFP),floor(length(tRegPossibleFP)/2))]) = true;

%
[timeMaster,sI] = sort(timeMaster);
labelMaster = labelMaster(sI);
nccMaster = nccMaster(sI);
neffMaster = neffMaster(sI);
istrain = istrain(sI);
toc;

% get possible features for both training/testing and tp/fp timepoints...
verboseFlag = true;
kstnms = ["PUYO";...
    "BULB";...
    "TAIS";...
    "TAMH";...
    "BMAS";...
    "BPAT";...
    "PKYU";...
    "PORT";...
    "BRUN"];

sensitivities = [3.141950e+08 4.872110e+08 3.141950e+08 3.141950e+08...
    4.872110e+08 4.872110e+08 3.141950e+08 3.141950e+08 4.872110e+08];
dists = [66700,63289,102391,71016,57654,56096,91971,76817,65576];

chans = ["HHZ";...
    "BHZ";...
    "HHZ";...
    "HHZ";...
    "BHZ";...
    "BHZ";...
    "HHZ";...
    "HHZ";...
    "BHZ"];

%
lKstnms = length(kstnms);
filterFlag = true;
filterObject = [0.6 1.2 false false];

%
lMaster = length(timeMaster);
secDur = 150;
newFs = 10;
extraFeatures = 2 + 5*lKstnms;
features = NaN(lMaster,lKstnms*(lKstnms-1)*0.5*2 + extraFeatures);
features(:,1) = nccMaster;
features(:,2) = neffMaster;

%clearvars -except labelMaster istrain timeMaster nccMaster neffMaster

%%
% cd ~/research/now/sangay/
% save('sangayRandomForest','labelMaster','istrain','timeMaster','nccMaster','neffMaster');

%%
% labels = labelMaster;
% si = 1;
% chunkSize = 200;
% tw = 0.1;
%
% clc;
% tic;
% while si <= lMaster - chunkSize + 1
%     %for i = 1
%     ei = si + chunkSize -1;
%     disp(timeMaster(ei));
%     tStart_ = timeMaster(si:ei);
%     disp([si ei]);
%     
%     %
%     data = [];
%     
%     %
%     for j = 1:lKstnms
%         kstnm_ = kstnms(j);
%         chan_ = chans(j);
%         disp(['Attempting to read data from: ',[char(kstnm_) char(chan_)]]);
%         S_ = extractWaveforms(tStart_,seconds(secDur),kstnm_,chan_,"EC","",...
%             true,false,'~/rawdata/',1,filterFlag,filterObject);
%         keepI = ~isnat(pull(S_,'ref'));
%         
%         %
%         data_ = NaN(secDur*newFs+1,chunkSize);
%         
%         %
%         if sum(keepI)
%             Skeep = S_(keepI);
%             dataKeep = detrend(double(pull(resampleWaveforms(Skeep,newFs))),1,'omitnan');
%             lkeep = size(dataKeep,1);
%             data_(1:lkeep,keepI) = taper(dataKeep,tw);
%         end
%         
%         %
%         data = [data; data_];
%         toc;
%     end
%     
%     %
%     for j = 1:chunkSize
%         data_ = data(:,j);
%         data_ = reshape(data_,[secDur*newFs+1,lKstnms]);
%         
%         %
%         [maxccp_,plags_,maxccn_,nlags_] = doccFreqCircShift(data_,false);
%         maxI = maxccp_ < maxccn_;
%         maxccp_(maxI) = maxccn_(maxI);
%         plags_(maxI) = nlags_(maxI);
%         
%         features_ = [max(abs(data_)) peak2rms(data_) kurtosis(data_,0)...
%             skewness(data_,0) rms(data_) maxccp_' plags_'];
%         features(si+j-1,3:end) = features_;
%         toc;
%     end
%     
%     %
%     labels_ = labelMaster(si:ei);
%     labels(si:ei) = labels_;
%     si = ei+1;
%     toc;
% end
