function [features,tFeatures] = updateAllSangayFeatures(saveFlag)
if nargin < 1
    saveFlag = false;
end

%%
FEATUREDIR = fullfile("~","igdata");
%FEATUREDIR = fullfile("~","research","now","sangay");
CATALOGDIR = FEATUREDIR;
longtermfileName = fullfile(CATALOGDIR,"SangayRegionalAnalysis_v8");
featureFile = fullfile(FEATUREDIR,"AllFeatures");
AF = load(featureFile);
features = AF.features;
tFeatures = AF.tFeatures;

%%
try
    %load('tmp_cat_throw');
    [tREG,~,nccREG,NeffREG] = filterUniqueEvents(longtermfileName);
catch
    pause(5);
    try
        [tREG,~,nccREG,NeffREG] = filterUniqueEvents(longtermfileName);
    catch
        fprintf(2,'something went wrong, couldnt read the catalog\n');
        features = [];
        tFeatures = [];
        return;
    end
end

%%
refTime = dn2dt(min([min(datenum(tREG)) min(datenum(tFeatures))]));
lia = ismembertol(seconds(tFeatures - refTime),seconds(tREG - refTime),1,'DataScale',1);
if sum(~lia)
    tFeatures(~lia) = [];
    features(~lia,:) = [];
end
[lia2,~] = ismembertol(seconds(tREG - refTime),seconds(tFeatures - refTime),1,'DataScale',1);

%%
if ~sum(~lia2)
    fprintf(1,'The dataset is already up to date (as far as I can tell), doing nothing\n');
    features = [];
    tFeatures = [];
    return
end

%% generate "query" events that are NOT in the tFeatures set
tREGq = tREG(~lia2);
nccREGq = nccREG(~lia2);
NeffREGq = NeffREG(~lia2);

%%
[features_,tFeatures_] = generateAllSangayFeatures(tREGq,nccREGq); %,NeffREGq);

%%
tFeatures = [tFeatures; tFeatures_];
features = [features; features_];

[tFeatures,sI] = sort(tFeatures);
features = features(sI,:);

%%
if saveFlag
    save('~/igdata/AllFeatures','features','tFeatures','-v7.3');
end