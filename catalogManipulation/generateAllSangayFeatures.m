function [features,tREG] = generateAllSangayFeatures(tREG,nccREG)
warning off MATLAB:detrend:PolyNotUnique

if ~nargin
    tStart = datetime(2000,01,01);
    try
        [tREG,~,nccREG,~] = filterUniqueEvents('~/research/now/sangay/SangayRegionalAnalysis_v6');
    catch
        pause(5);
        try
            [tREG,~,nccREG,~] = filterUniqueEvents('~/research/now/sangay/SangayRegionalAnalysis_v6');
        catch
            fprintf(2,'something went wrong, couldnt read the catalog\n');
            features = [];
            tREG = [];
            return;
        end
    end

    %%
    tI = tREG > tStart;

    %%
    tREG = tREG(tI);
    nccREG = nccREG(tI);
    %NeffREG = NeffReg(tI);
elseif nargin < 2
    features = [];
    tREG = [];
    fprintf('not enough input arguments\n');
    return;
end

%% get possible features for both training/testing and tp/fp timepoints...
kstnms = ["PUYO";...
    "BULB";...
    "TAIS";...
    "TAMH";...
    "BMAS";...
    "BPAT";...
    "PKYU";...
    "PORT";...
    "BRUN"];

% chans = ["HHZ";...
%     "BHZ";...
%     "HHZ";...
%     "HHZ";...
%     "BHZ";...
%     "BHZ";...
%     "HHZ";...
%     "HHZ";...
%     "BHZ"];

%%
lKstnms = length(kstnms);

%%
secDur = 150;
newFs = 10;
extraFeatures = 2 + 5*lKstnms;
tw = 100;
lfc = 0.6;
hfc = 1.2;
filterFlag = true;
filterObject = [lfc hfc false false false newFs];
winlen = secDur*newFs+1;

%%
loopDays = (unique(dateshift(tREG,'start','day'))); %dn2dt(unique(floor(datenum(tREG))));
lDays = length(loopDays);
featuresCell = cell(lDays,1);
levents = length(tREG);

%%
fprintf('\n');
fprintf('extracting features for %d events\n',sum(levents));
fprintf('total number of days to read: %d\n',lDays);
fprintf('\n');

%%
for i = 1:lDays
    dayStart = loopDays(i);
    

    %%
    tcutI = tREG >= dayStart & tREG < dayStart + 1;
    nEvents = sum(tcutI);
    if ~nEvents
        continue;
    end

    %%
    fprintf('processing day %s with %d events\n',dayStart,nEvents);
    tStart_ = tREG(tcutI);
    featuresTmp = NaN(nEvents,lKstnms*(lKstnms-1)*0.5*2 + extraFeatures);
    featuresTmp(:,1) = nccREG(tcutI);
    %featuresTmp(:,2) = NeffREG(tcutI);

    %%
    data = NaN(winlen,nEvents,lKstnms);
    chans = true(nEvents,lKstnms);
    dayStart.Format = 'yyyyMMdd';
    tic;
    for j = 1:lKstnms
        kstnm_ = kstnms(j);

        %%
        S_ = extractWaveforms(tStart_,seconds(secDur),kstnm_,["BHZ";"HHZ"],"EC","",...
            true,false,1,filterFlag,filterObject);

        S_ = S_(:,1);
        chan_ = pull(S_,'kcmpnm');

        uniqChan = unique(chan_);
        uniqChan = uniqChan(end);
        if strcmp(uniqChan,"")
            fprintf('no data for %s, filling with nans\n',kstnm_)
            chans(:,j) = false;
        else
            chans(:,j) = strcmp(chan_,uniqChan);
        end

        keepI = ~isnat(pull(S_,'ref'));
        data_ = NaN(winlen,nEvents);

        %%
        if sum(keepI) %if at least 1 event for this station...
            Skeep = S_(keepI);
            dataKeep = detrend(double(pull(Skeep)),1,'omitnan');
            lkeep = size(dataKeep,1);
            data_(1:lkeep,keepI) = taper(dataKeep,tw);
        end

        %%
        data(:,:,j) = data_; 
    end

    %%
    Neff2 = sum(chans,2,'omitnan');
    for j = 1:nEvents
        data_ = squeeze(data(:,j,:));
        %data_ = reshape(data_,[secDur*newFs+1,lKstnms]);

        %
        [maxccp_,plags_,maxccn_,nlags_] = doccFreqCircShift(data_,false);
        maxI = maxccp_ < maxccn_;
        maxccp_(maxI) = maxccn_(maxI);
        plags_(maxI) = nlags_(maxI);

        features_ = [max(abs(data_)) peak2rms(data_) kurtosis(data_,0)...
            skewness(data_,0) rms(data_) maxccp_' plags_'];
        featuresTmp(j,3:end) = features_;
        featuresTmp(j,2) = Neff2(j); %NeffREG(tcutI);
    end

    %
    %fname = sprintf('featuresTmp_%s.mat',dayStart);
    %fname = fullfile('~','data','tmp',fname);
    %save(fname,'featuresTmp','-v7.3');
    featuresCell{i} = featuresTmp;
    toc;
end

%%
features = cat(1,featuresCell{:});
