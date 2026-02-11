clear; close all;

%
cd ~/research/now/sangay/
load sangayRandomForest.mat;
clear features;
AF = load('~/research/now/sangay/AllFeatures.mat');
timeAll = AF.tFeatures;
features = AF.features;
refTime = dn2dt(min([min(datenum(timeAll)) min(datenum(timeMain))]));

%
lia = ismembertol(seconds(timeAll - refTime),seconds(timeMain - refTime),3,'DataScale',1);
features = features(lia,:);
lia = ismembertol(seconds(timeMain - refTime),seconds(timeAll - refTime),3,'DataScale',1);
timeMain = timeMain(lia);
labelMain = labelMain(lia);
nccMain = nccMain(lia);
istrain = istrain(lia);
neffMaster = neffMaster(lia);
clear lia AF timeAll refTime

%%  legend (9 of each)
%z2p [3:11]
%peak2rms [12:20]
%kurtosis [21:29]
%skewness [30:38]
%rms [39:47]

%% experimental
exerimentalFlag = false;
if exerimentalFlag
    load('sangayGentleBoostCompact');
end

%% PUYO BULB TAIS TAMH BMAS BPAT PKYU PORT BRUN
dists = [66700,63289,102391,71016,57654,56096,91971,76817,65576];

sensitivities = [3.141950e+08 4.872110e+08 3.141950e+08 3.141950e+08...
    4.872110e+08 4.872110e+08 3.141950e+08 3.141950e+08 4.872110e+08]; % ./dists; %dists adds no information

predGood = true(119,1); %predGood(39:47) = false; %features(:,3:11) = features(:,39:47);

%%
neff = features(:,2);
tGood = neff >= 2;
NumLearningCycles = 200;
classificationAlgorithm = 'GentleBoost'; %'AdaBoostM1'; %'RobustBoost'; %'GentleBoost'; %others might include 'AdaBoostM1'

%%
features = features(tGood,predGood);
timeMain = timeMain(tGood);
nccMain = nccMain(tGood);
labelMain = labelMain(tGood);
istrain = istrain(tGood);
neff = neff(tGood);

%%
% features = [features(:,[1:11 84:119]) Msangay Err_sangay];
% log10(1./(abs(getRatios(featureSingle(:,3:11)',false))'))];

%% 552/674, t = templateTree('MaxNumSplits',10,'Surrogate','all','NumVariablesToSample','all');
% [Msangay,Err_sangay,mlv] = sangayMagnitudeCalculation(timeMaster,features(:,3:11));
% %Err_sangay = mean(abs(getRatios(mlv',false)),1,'omitnan')';
%
% features(:,3:11) = mlv;
% features(:,39:47) = log10(features(:,39:47)./sensitivities);   	% scale rms
% featureSingle = features(:,1:47);
% features(:,1:47) = [];
% features = [featureSingle features ...
%     log10(getRatios(featureSingle(:,3:11)',true)')...         	% mlv
%     (getRatios(featureSingle(:,12:20)',false)')...              % peak2rms (false is best)
%     (getRatios(featureSingle(:,21:29)',false)')...              % kurtosis
%     (getRatios(featureSingle(:,30:38)',false)')...              % skewness
%     (getRatios(featureSingle(:,39:47)',false)')...            	% rms (false is best)
%     Msangay Err_sangay];

%%
% [Msangay,Err_sangay,mlv,wa] = sangayMagnitudeCalculation(timeMaster,features(:,3:11));
% 
% features(:,3:11) = mlv;
% features(:,39:47) = log10(wa);
% featureSingle = features(:,1:47);
% features(:,1:47) = [];
% features = [featureSingle features ...
%     (getRatios(featureSingle(:,3:11)',false)')...               % mlv
%     (getRatios(featureSingle(:,12:20)',false)')...              % peak2rms (false is best)
%     (getRatios(featureSingle(:,21:29)',false)')...              % kurtosis
%     (getRatios(featureSingle(:,30:38)',false)')...              % skewness
%     (getRatios(featureSingle(:,39:47)',false)')...              % rms (false is best)
%     Msangay Err_sangay ...
%     nanmedian(abs(getRatios(wa',false)))' ...
%     nanvar(getRatios(wa',false))'];
% %(nanvar(getRatios(wa',false)))'

%%
[Msangay,Err_sangay,mlv,wa] = sangayMagnitudeCalculation(timeMain,features(:,3:11));
features(:,3:11) = mlv;
features(:,39:47) = wa; %log10(features(:,39:47)./sensitivities);   	% scale rms
featureSingle = features(:,1:47);
features(:,1:47) = [];
features = [featureSingle features ...
    (getRatios(mlv',false)')...                                 % mlv
    (getRatios(featureSingle(:,12:20)',false)')...              % peak2rms (false is best)
    (getRatios(featureSingle(:,21:29)',false)')...              % kurtosis
    (getRatios(featureSingle(:,30:38)',false)')...              % skewness
    (getRatios(featureSingle(:,39:47)',false)')...            	% rms (false is best)
    Msangay Err_sangay ...
    nanmedian(log10(getRatios(wa',false).^2))'];

%%

%??, 50 learners, ~??%, 156 features, 50 maxSplits, predictionSelector:interaction-curvature

%109k, 50 learners, ~95.6%, 156 features, 50 maxSplits
%featureSingle = [features(:,1:46) Msangay Err_sangay];
%features(:,1:46) = [];
%features = [features featureSingle (getRatios(featureSingle(:,2:10)',false)')];

%106k, 200 learners, ~96.1%, 120 features, 20 maxSplits
%featureSingle = [features(:,1:46) Msangay Err_sangay];
%features(:,1:46) = [];
%features = [features featureSingle];

%114k, 200 learners, ~95.5%, 48 features, 20 maxSplits
%features = [features(:,1:46) Msangay Err_sangay];

%111155, 200 learners, ~96.3%, 300 features, 20 maxSplits
%featureSingle = [features(:,1:46) Msangay Err_sangay];
%features(:,1:46) = [];
%features = [features featureSingle];
%    (getRatios(featureSingle(:,2:10)',false)')...   % mlv
%    (getRatios(featureSingle(:,11:19)')')...        % peak2rms
%    (getRatios(featureSingle(:,20:28)')')...        % kurtosis
%    (getRatios(featureSingle(:,29:37)',false)')...  % skewness
%    (getRatios(featureSingle(:,38:46)')')];         % rms

%(getRatios(features(:,2:10)')') (getRatios(features(:,38:46)')')];
%features = [features log10(getRatios(features(:,2:10)')') log10(getRatios(features(:,38:46)')')];
disp(size(features));

%% experimental
if exerimentalFlag
    modelFits = predict(cGentleBoostEnsemble,features);
    confusedIndex = find(modelFits ~= labelMain & istrain);
    features(confusedIndex,:) = [];
    %labels(confusedIndex) = [];
    istrain(confusedIndex) = [];
    labelMain(confusedIndex) = [];
    timeMain(confusedIndex) = [];
    nccMain(confusedIndex) = [];
end

%%
featuresTrain = features(istrain,:);
labelsTrain = labelMain(istrain);

featuresTest = features(~istrain,:);
labelsTest = labelMain(~istrain);

%%
%t = templateTree('Surrogate','all','NumVariablesToSample','all');
%t = templateTree('MaxNumSplits',20,'Surrogate','all','NumVariablesToSample','all');

t = templateTree('MaxNumSplits',10,'Surrogate','all','NumVariablesToSample','all');
%t = templateTree('MaxNumSplits',3,'Surrogate','all','NumVariablesToSample','all','PredictorSelection','curvature');
%t = templateTree('MinLeafSize',20,'Surrogate','all','NumVariablesToSample','all','PredictorSelection','interaction-curvature');
%t = templateTree('Surrogate','all','NumVariablesToSample','all','PredictorSelection','interaction-curvature');

tic;
gentleBoostEnsemble = fitcensemble(featuresTrain,labelsTrain,...
    'Type','classfication',...
    'Method',classificationAlgorithm,...
    'NumLearningCycles',NumLearningCycles,...
    'Learners',t,...
    'nprint',25,...
    'ScoreTransform','doublelogit');
toc;

%%
figure;
tic;
plot(loss(gentleBoostEnsemble,featuresTrain,labelsTrain,'mode','cumulative'));

figure(1);
hold on;
tic;
plot(loss(gentleBoostEnsemble,featuresTest,labelsTest,'mode','cumulative'));
toc;

grid on;
xlabel('Number of trees');
ylabel('Test classification error'); zoom on;

%%
close all;
disp('applying trained classifaction ensemble to testing dataset...')
Yfit = predict(gentleBoostEnsemble,featuresTest);

disp('results of final classification:');
testTrue = labelsTest == true;
testFalse = labelsTest == false;
yfittrue = Yfit == true;
yfitfalse = Yfit == false;

truetrue = testTrue & yfittrue;
falsefalse = testFalse & yfitfalse;
falsetrue = testFalse & yfittrue;
truefalse = testTrue & yfitfalse;

fprintf('true positive rate: <strong>%f</strong> (%d/%d)\n',sum(truetrue)/sum(testTrue),sum(truetrue),sum(testTrue));
fprintf('true negative rate: <strong>%f</strong> (%d/%d)\n',sum(falsefalse)/sum(testFalse),sum(falsefalse),sum(testFalse));
fprintf('false positive rate: <strong>%f</strong> (%d/%d)\n',sum(falsetrue)/sum(testFalse),sum(falsetrue),sum(testFalse));
fprintf('false negative rate: <strong>%f</strong> (%d/%d)\n',sum(truefalse)/sum(testTrue),sum(truefalse),sum(testTrue));

%%
tt = timeMain(~istrain);
mm = Msangay(~istrain);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
axx(1) = subplot(211);
plot(axx(1),tt(truetrue),mm(truetrue),'.'); hold on; zoom on; grid on; hold on;
plot(axx(1),tt(falsetrue),mm(falsetrue),'.'); hold on; zoom on; grid on;
%title('events i keep and should have');
legend('good keep','bad keep','location','northwest');

axx(2) = subplot(212);
plot(axx(2),tt(falsefalse),mm(falsefalse),'.'); hold on; zoom on; grid on; hold on;
plot(axx(2),tt(truefalse),mm(truefalse),'.'); hold on; zoom on; grid on;
%title('events i discard and should have');
legend('good discard','bad discard','location','northwest');

%axx(3) = subplot(211);
%plot(tt(falsetrue),mm(falsetrue),'.'); hold on; zoom on; grid on;
%title(['events i keep but should have discarded, N = ',num2str(sum(falsetrue))]);


%axx(4) = subplot(414);
%plot(tt(truefalse),mm(truefalse),'.'); hold on; zoom on; grid on;
%title(['events i discard but should have kept, N = ',num2str(sum(truefalse))]);

linkaxes(axx,'x');

%%
if ~exerimentalFlag
    cGentleBoostEnsemble = compact(gentleBoostEnsemble);
    tic;
    [yfit_,score] = predict(cGentleBoostEnsemble,featuresTest(end,:));
    toc;
    save('sangayGentleBoostCompact','cGentleBoostEnsemble');
    !\ls -lhrt sangayGentleBoostCompact.mat
end

%%
imp = predictorImportance(cGentleBoostEnsemble);
figure;
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
zoom on;
figure();
plot(1:length(imp),cumsum(sort(imp,'descend'))/sum(imp),'.');
zoom on; grid on;
