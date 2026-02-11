%CotopaxiGenerateClassificationModel
clear; clear all; close all;

%
tic;
cd ~/research/now/cotopaxi/tremor_binary_classification;
load("CotopaxiFeatures_v1");

%%
trainI = (1:2:maxFeatures)';
validateI = (2:2:maxFeatures)';

%
featuresTrain = [featuresFalse(trainI,:); ...
    featuresTrue(trainI,:)];

featuresValidate = [featuresFalse(validateI,:); ...
    featuresTrue(validateI,:)];

%
labelsTrain = [false(length(trainI),1); true(length(trainI),1)];
labelsValidate = [false(length(validateI),1); true(length(validateI),1)];

%
timeTrain = [timeFalse(trainI); timeTrue(trainI)];
timeValidate = [timeFalse(validateI); timeTrue(validateI)];

%
[~,sortI] = sort(timeTrain);
timeTrain = timeTrain(sortI);
labelsTrain = labelsTrain(sortI);
featuresTrain = featuresTrain(sortI,:);

%
[~,sortI] = sort(timeValidate);
timeValidate = timeValidate(sortI);
labelsValidate = labelsValidate(sortI);
featuresValidate = featuresValidate(sortI,:);
toc;

%%
t = templateTree('MaxNumSplits',5,'Surrogate','all','NumVariablesToSample','all');
NumLearningCycles = 400;
classificationAlgorithm = 'GentleBoost';

tic;
gentleBoostEnsemble = fitcensemble(featuresTrain,labelsTrain,...
    'Type','classfication',...
    'Method',classificationAlgorithm,...
    'NumLearningCycles',NumLearningCycles,...
    'Learners',t,...
    'nprint',5,...
    'ScoreTransform','doublelogit');
toc;

%%
close all;
figure();
tic;
plot(loss(gentleBoostEnsemble,featuresTrain,labelsTrain,'mode','cumulative'));
hold on;
toc;
plot(loss(gentleBoostEnsemble,featuresValidate,labelsValidate,'mode','cumulative'));
toc;

grid on;
xlabel('Number of trees');
ylabel('Test classification error'); zoom on;

%
disp('applying trained classifaction ensemble to testing dataset...');
tic;
Yfit = predict(gentleBoostEnsemble,featuresValidate);
toc;

disp('results of final classification:');
validateTrue = labelsValidate == true;
validateFalse = labelsValidate == false;
yfittrue = Yfit == true;
yfitfalse = Yfit == false;

truetrue = validateTrue & yfittrue;
falsefalse = validateFalse & yfitfalse;
falsetrue = validateFalse & yfittrue;
truefalse = validateTrue & yfitfalse;

fprintf('true positive rate: <strong>%f</strong> (%d/%d)\n',sum(truetrue)/sum(validateTrue),sum(truetrue),sum(validateTrue));
fprintf('true negative rate: <strong>%f</strong> (%d/%d)\n',sum(falsefalse)/sum(validateFalse),sum(falsefalse),sum(validateFalse));
fprintf('false positive rate: <strong>%f</strong> (%d/%d)\n',sum(falsetrue)/sum(validateFalse),sum(falsetrue),sum(validateFalse));
fprintf('false negative rate: <strong>%f</strong> (%d/%d)\n',sum(truefalse)/sum(validateTrue),sum(truefalse),sum(validateTrue));

%%
cGentleBoostEnsemble = compact(gentleBoostEnsemble);
tic;
[yfit_,score] = predict(cGentleBoostEnsemble,featuresValidate(end,:));
toc;

saveFlag = ~true;
if saveFlag
    clearvars -except labelsTrain labelsValidate featuresTrain featuresValidate cGentleBoostEnsemble kstnms newFxx
    save('CotopaxiGentleBoostEnsembleCompact_5MaxNumSplits400LearningCycles');
end

%%
% % kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];
% % lfc = 1/5;
% % hfc = 10;
% % newFs = 25;
% % nF = 101;
% % newFxx = logspace(log10(lfc),log10(hfc),nF)';
% 
% imp = predictorImportance(cGentleBoostEnsemble);
% figure;
% bar(imp);
% title('Predictor Importance Estimates');
% ylabel('Estimates');
% xlabel('Predictors');
% zoom on;
% figure();
% plot(1:length(imp),cumsum(sort(imp,'descend'))/sum(imp),'.');
% zoom on; grid on;
% 
% 
% close all;
% imp2 = reshape(imp*1e3,101,55);
% clear ax;
% nCol = 0;
% goodFeaturesI = false(101,55);
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(kstnms)
%     kstnm1 = kstnms(i);
%     for ii = i+1:length(kstnms)
%         kstnm2 = kstnms(ii);
%         comboStr = sprintf("%s-%s",kstnm1,kstnm2);
%         nCol = nCol+1;
%         
%         imp2_ = imp2(:,nCol);
%         [imp3,sortI] = sort(imp2_,'descend');
% 
%         % old numbers: 1.2, 10 (670)
%         % new numbers: 1.0, 3
%         thresh1 = 3;
%         thresh5 = imp3(10);
% 
%         ax(nCol) = subplot(5,11,nCol);
%         iI = imp2_ >= thresh1 | imp2_ >= thresh5;
%         loglog(newFxx(~iI),imp2_(~iI),'.'); zoom on; grid on;
%         hold on;
%         minGoodThresh = min(imp2_(iI));
%         loglog(newFxx(iI),imp2_(iI),'.'); zoom on; grid on;
%         ylim([0.3 2]); xlim([0.2 10]); 
%         plot([0.2 10],thresh1*[1 1],'-','linewidth',2);
%         title(comboStr,'FontSize',11);
%         legStrTmp = sprintf("%d,%3.2f",sum(iI),minGoodThresh);
%         text(0.3,1.7,legStrTmp,'FontSize',12);
%         goodFeaturesI(:,nCol) = iI;
%     end
% end
% 
% linkaxes(ax,'xy');
% goodFeaturesI = goodFeaturesI(:);
% sgtitle(sprintf("Total Number Good Features: %d, Threshold for Extras: %3.2f",sum(goodFeaturesI),thresh1));
% disp(sum(goodFeaturesI));
