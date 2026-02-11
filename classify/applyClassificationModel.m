function [YfitAll,t,eqmag,ncc,energyMag,merr,tFeatures,NCC,Neff,Msangay] = ...
    applyClassificationModel(cModel,features,tFeatures,minNeff,minMag)
if nargin < 4
    minNeff = 3;
end

if nargin < 5
    minMag = 1;
end

%%
predGood = true(119,1);
features = features(:,predGood);

%%
[Msangay,Err_sangay,mlv,wa] = sangayMagnitudeCalculation(tFeatures,features(:,3:11));
features(:,3:11) = mlv;
features(:,39:47) = wa;
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
NCC = features(:,1);
Neff = features(:,2);

%%
tic;
fprintf('generating predictions\n');
[YfitAll,~] = predict(cModel,features);
toc;

%%
YfitAll = (YfitAll & Neff >= minNeff & Msangay <= 3.15 & Msangay >= minMag);

%%
t = tFeatures(YfitAll);
ncc = NCC(YfitAll);
merr = Err_sangay(YfitAll);
eqmag = Msangay(YfitAll);
energyMag = 10.^(1.5.*eqmag+4.8);