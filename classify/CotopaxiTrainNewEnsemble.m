function [newEnsemble,goodFeaturesI,newGentleBoostEnsemble] = CotopaxiTrainNewEnsemble(labelsTrain,labelsValidate,featuresTrain,featuresValidate,...
    mainEnsemble,param1,param2,MaxNumSplits,NumLearningCycles,classificationAlgorithm)

if nargin < 8; MaxNumSplits = 10; end
if nargin < 9; NumLearningCycles = 25; end
if nargin < 10; classificationAlgorithm = 'GentleBoost'; end

%%
kstnms = ["CO1V";"BREF";"BVC2";"BTAM";"BMOR";"BNAS";"VC1";"NAS2";"SUCR";"SLOR";"TAMB"];
lfc = 1/5;
hfc = 10;
nF = 101;
newFxx = logspace(log10(lfc),log10(hfc),nF)';

imp = predictorImportance(mainEnsemble);
figure;
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
zoom on;
figure();
plot(1:length(imp),cumsum(sort(imp,'descend'))/sum(imp),'.');
zoom on; grid on;

%
imp2 = reshape(imp*1e3,101,55);
clear ax;
nCol = 0;
goodFeaturesI = false(101,55);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(kstnms)
    kstnm1 = kstnms(i);
    for ii = i+1:length(kstnms)
        kstnm2 = kstnms(ii);
        comboStr = sprintf("%s-%s",kstnm1,kstnm2);
        nCol = nCol+1;

        imp2_ = imp2(:,nCol);
        [imp3,~] = sort(imp2_,'descend');

        %old numbers: 1.2, 10 (670)
        %new numbers: 1.0, 3
        thresh1 = param2;
        thresh5 = imp3(param1);

        ax(nCol) = subplot(5,11,nCol);
        iI = imp2_ >= thresh1 | imp2_ >= thresh5;
        loglog(newFxx(~iI),imp2_(~iI),'.'); zoom on; grid on;
        hold on;
        minGoodThresh = min(imp2_(iI));
        loglog(newFxx(iI),imp2_(iI),'.'); zoom on; grid on;
        ylim([0.3 2]); xlim([0.2 10]);
        plot([0.2 10],thresh1*[1 1],'-','linewidth',2);
        title(comboStr,'FontSize',11);
        legStrTmp = sprintf("%d,%3.2f",sum(iI),minGoodThresh);
        text(0.3,1.7,legStrTmp,'FontSize',12);
        goodFeaturesI(:,nCol) = iI;
    end
end

linkaxes(ax,'xy');
goodFeaturesI = goodFeaturesI(:);
sgtitle(fig,sprintf("Total Number Good Features: %d, Threshold for Extras: %3.2f",sum(goodFeaturesI),thresh1));
fprintf('number of good features: <strong>%d</strong>\n',sum(goodFeaturesI));

%%
t = templateTree('MaxNumSplits',MaxNumSplits,'Surrogate','all','NumVariablesToSample','all');
tic;
newGentleBoostEnsemble = fitcensemble(featuresTrain(:,goodFeaturesI),labelsTrain,...
    'Type','classfication',...
    'Method',classificationAlgorithm,...
    'NumLearningCycles',NumLearningCycles,...           
    'Learners',t,...
    'nprint',5,...
    'ScoreTransform','doublelogit');
toc;

%%
figure(); hold on;
tic;
plot(loss(newGentleBoostEnsemble,featuresTrain(:,goodFeaturesI),labelsTrain,'mode','cumulative'));
toc;
plot(loss(newGentleBoostEnsemble,featuresValidate(:,goodFeaturesI),labelsValidate,'mode','cumulative'));
toc;
grid on;
xlabel('Number of trees');
ylabel('Test classification error'); zoom on;

%
fprintf('applying trained classifaction ensemble to testing dataset...\n');
tic; 
Yfit = predict(newGentleBoostEnsemble,featuresValidate(:,goodFeaturesI)); 
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
newEnsemble = compact(newGentleBoostEnsemble);
% tic;
% [yfit_,score] = predict(newEnsemble,featuresValidate(end,goodFeaturesI));
% toc;
%

imp = predictorImportance(newEnsemble);
figure;
bar(imp);
title('Predictor Importance Estimates');
ylabel('Estimates');
xlabel('Predictors');
zoom on;
figure();
plot(1:length(imp),cumsum(sort(imp,'descend'))/sum(imp),'.');
zoom on; grid on;
