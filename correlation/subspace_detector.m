function [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
    obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
    subspace_detector(S,varargin)

tic;
warning off signal:findpeaks:largeMinPeakHeight

nVarargin = length(varargin);
functionDefaults = {...
    0.07,...
    '~/igdata/CotopaxiMultiplexedBasisFunctions.mat',...
    20,...
    5e3,...
    true,...
    false,...
    true,...
    "cvccn_drumbeats"};

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[threshold,basisFunctionFileName,maxBasisFunctions,maxEvents,writeFlag,...
    diffFlag,linearccnorm,customPrefix] = deal(optsToUse{:});

%%
customPrefixFlag = true;
if isempty(customPrefix)
    customPrefixFlag = false;
end

if isempty(diffFlag)
    diffFlag = false;
end

%%
if isstruct(basisFunctionFileName)
    basisFunctions = basisFunctionFileName;
else
    basisFunctions = load(basisFunctionFileName);
end
clear basisFunctionFileName;

lfc = basisFunctions.lfc;
hfc = basisFunctions.hfc;
B = basisFunctions.U;
B = B(1:maxBasisFunctions,:);

%%
Bknetwk = pull(B(1,:)','knetwk');
Bkstnm = pull(B(1,:)','kstnm');
Bkhole = pull(B(1,:)','khole');
BallSNCLs = strcat(Bknetwk,Bkstnm,Bkhole);
BFs = median(1./pull(B,'delta'),'all');

%%
if diffFlag
    S = differentiateWaveforms(S);
end
S = filterWaveforms(detrendWaveforms(S),lfc,hfc);
S = resampleWaveforms(S,BFs);
S = syncWaveforms(S,false,true,true);
S = nanGapWaveforms(S,NaN);

%%
npts_ = pull(S,'npts');
t = getTimeVec(S(find(isfinite(npts_),1))); %<-- hack
t = repmat(t,[1 3]);
t = t';
t = t(:);

%%
Sknetwk = pull(S,'knetwk');
Skstnm = pull(S,'kstnm');
Skhole = pull(S,'khole');
SallSNCLs = strcat(Sknetwk,Skstnm,Skhole);

%%
lB = length(BallSNCLs);
energyRatio = NaN(maxEvents,1);
Neff = energyRatio;
tabs = NaT(maxEvents,1);
ReconstructCoefficients = zeros(maxBasisFunctions*lB,maxEvents);

if lB == 1
    z2p = energyRatio;
    p2rms = energyRatio;
    kurt = energyRatio;
    skew = energyRatio;
else
    z2p = NaN(maxEvents,lB);
    p2rms = z2p;
    kurt = z2p;
    skew = z2p;
end

%%
eratio = NaN(median(npts_)*3,lB);
dLong = eratio;
goodI = false(lB,1);

%%
for i = 1:lB %loop over each monitoring site
    BSNCLs_ = BallSNCLs(i);
    lia = ismember(SallSNCLs,BSNCLs_);
    sumlia = sum(lia);
    if ~sumlia
        fprintf(1,'no data for %s\n',BSNCLs_);
        continue;
    end

    S_ = S(lia);
    dLong_ = pull(S_);
    thisNChan = size(dLong_);
    if thisNChan ~= 3
        fprintf(2,'not enough channels for %s, read more, error\n',BSNCLs_);
        % i can potentially fix this by adding zeros, but will be
        % complicated...
        continue;
    end

    theseU = pull(B(:,i));
    theseU = flipud(theseU);
    winlen = size(theseU,1);
    box = ones(winlen,1);

    dLong_ = dLong_';
    dLong_ = dLong_(:);

    badI = ~isfinite(dLong_);
    dLong_(badI) = 0;
    
    dLong2 = dLong_.^2;

    %%
    norms = fftfilt(box,dLong2); %sum of squares for each time window
    snorms = norms > 1;

    cc_ = fftfilt(theseU,dLong_); %do all basis functions at once
    cc_ = sum(cc_.^2,2); %sum energy controbutions across all nBF basis functions (this channel only)
    eratio_ = snorms.*cc_./norms; %akin to normalized cross correlation....

    badI = badI | ~isfinite(eratio_) | ~snorms;
    eratio_(badI) = 0;
    eratio_(badI) = NaN;
    eratio(:,i) = eratio_;

    dLong_(badI) = NaN;
    dLong(:,i) = dLong_;
    goodI(i) = true;
    
end
eratio = eratio(:,goodI);
eratioOrig = eratio;
mpd = 0.25*winlen;
ldf2 = size(dLong_,1);
maxIndex = ldf2 - winlen + 1;

%% compress ccnorm
notbadI = isfinite(eratio);
nEffective = sum(notbadI,2);

%% collapse eratio (weighted sum or median?)
if lB > 1
    % eratio is a vector after this step
    if linearccnorm
        %fprintf("sum weights: %f\n",sum(weights(isfinite(weights)),"all"));
        %eratio = sum(weights.*eratio,2,"omitnan");      %weighted sum, weights should (hopefully) sum to 1
        eratio = mean(eratio,2,"omitnan");
    else
        eratio = median(eratio,2,"omitnan");
    end
end

%%
eratio(~isfinite(eratio)) = 0; % fix NaNs

%%
if isempty(threshold)
    %threshold = median(eratio,"omitnan") + 8*mad(eratio,1);
    threshold = 10*mad(eratio,1);
end
[CC,locs_] = findpeaks(eratio,'MINPEAKHEIGHT',threshold,'MINPEAKDISTANCE',mpd); %,'Threshold',1e-4)
locs_ = locs_-winlen+1;

%%
lI = locs_ <= maxIndex & isfinite(CC) & locs_ > 0;
locs_ = locs_(lI);
CC = CC(lI);
tabs_ = t(locs_);
ll = length(locs_);
neff_ = nEffective(locs_+winlen-1);

%%
if ~ll
    obsAmpRatio = [];
    synthAmpRatio = [];
    fprintf(1,'0 events detected\n');
    return;
end

%%
fprintf('----\n');
fprintf('found: %d event(s)\n',ll);
fprintf('----\n');

%%
Neff(1:ll) = neff_;
tabs(1:ll) = tabs_;
energyRatio(1:ll) = CC;
z2pSynth = z2p;
for i = 1:lB
    theseU = pull(B(:,i));
    d_ = dLong(:,i);
    D = NaN(winlen,ll);
    for j = 1:ll
        tmp = d_(locs_(j):locs_(j)+winlen-1); %filtered data
        D(:,j) = tmp;
        tmp = reshape(tmp,[3,winlen/3])';
        medAmp_ = median(max(abs(tmp),[],1,'omitnan'));
        z2p(j,i) = medAmp_;
    end

    p2rms(1:ll,i) = peak2rms(D)';
    kurt(1:ll,i) = -3+kurtosis(D,0)';
    skew(1:ll,i) = skewness(D,0)';

    D = D';
    ReconstructCoefficients_ = D*theseU;
    ReconstructCoefficients_ = ReconstructCoefficients_';
    si = 1 + (i-1)*maxBasisFunctions;
    ei = i*maxBasisFunctions;
    ReconstructCoefficients(si:ei,1:ll) = ReconstructCoefficients_;
    synthD_ = theseU*ReconstructCoefficients_;
    z2pSynth(1:ll,i) = 0.5*peak2peak(synthD_)';
end

%%
tabs = tabs(1:ll);
energyRatio = energyRatio(1:ll);
Neff = Neff(1:ll);
z2p = z2p(1:ll,:);
z2pSynth = z2pSynth(1:ll,:);
p2rms = p2rms(1:ll,:);
kurt = kurt(1:ll,:);
skew = skew(1:ll,:);
ReconstructCoefficients = ReconstructCoefficients(:,1:ll);
obsAmpRatio = getDD(z2p',false);        %true is subtraction, false is division
synthAmpRatio = getDD(z2pSynth',false); %true is subtraction, false is division

%%
if ~writeFlag
    return;
end
currentDay = dateshift(S(1).ref,'start','day');
currentDay.Format = "yyyy.MM.dd";
rootDir = fullfile("~","masa","subspace_detector");

if customPrefixFlag
    prefix = fullfile(rootDir,customPrefix);
else
    kstnm = S(1).kstnm;
    prefix = fullfile(rootDir,kstnm);
end

if ~exist(prefix,"dir")
    mkdir(prefix);
end

fName = fullfile(prefix,sprintf("daySubspaceSearch_%s.%s",currentDay,"mat"));
try
    fprintf('attempting to write %d events on day: %s...\n',ll,currentDay);
    save(fName,'tabs','energyRatio','Neff','z2p','z2pSynth','p2rms',...
        'kurt','skew','ReconstructCoefficients','obsAmpRatio','synthAmpRatio','-v7.3');
catch
    fprintf("couldnt write file for day: %s, may need to reprocess...!\n",currentDay);
    return;
end
