function updateCotopaxiTremorCatalogBinaryClassifier()
%clear; close all;
%tremorHome = '~/igdata';
tremorHome = fullfile("~","masa","old","research","now","cotopaxi");

cd(tremorHome);
%saveFile = 'tremor_binary_classification/CotopaxiTremorCatalog_BinaryClassifier';
%modelFile = 'tremor_binary_classification/CotopaxiGentleBoostEnsembleCompact_v4';
%saveFileName = fullfile(tremorHome,"CotopaxiTremorCatalog_BinaryClassifier_M5L50");
saveFileName = fullfile(tremorHome,"CotopaxiTremorCatalog_BinaryClassifier_M5L50_2025BB");

%%
classifierModelFileName = "CotopaxiM5L50Classifier_GentleBoostEnsembleCompact.mat"; %'test04'; %<-- test04 == M5L50

%%
tic;
load(saveFileName);
toc; 

%%
%saveFileName = fullfile(tremorHome,'CotopaxiTremorCatalog_BinaryClassifier_M5L50');
saveFileName = fullfile(tremorHome,"CotopaxiTremorCatalog_BinaryClassifier_M5L50_2025BB");
load(classifierModelFileName);
experimentalFlag = true;

%%
tStart = dateshift(max(featureTime),'start','day');
tEnd = dateshift(datetime('now')+hours(5),'start','day');
deleteI = featureTime >= tStart;

allAmps(:,deleteI)= [];
allFeatures(deleteI,:)= [];
featureTime(deleteI)= [];
nGood(deleteI)= [];
yesScore(deleteI)= [];

%%
[allAmps_,allFeatures_,featureTime_,nGood_,yesScore_] = ...
    runUpdateCotopaxiTremorCatalogBinaryClassifier(tStart,tEnd,...
    cGentleBoostEnsemble2,goodFeaturesI,experimentalFlag);

%%
allAmps = [allAmps allAmps_];
allFeatures = [allFeatures; allFeatures_];
featureTime = [featureTime; featureTime_];
nGood = [nGood; nGood_];
yesScore = [yesScore; yesScore_];

%%
saveFlag = true;
if saveFlag
    %clearvars -except allFeatures allAmps featureTime yesScore nGood saveFile
    save(saveFileName,'allFeatures','allAmps','featureTime','yesScore','nGood','-v7.3');
end