clear; close all;

%tic; updateAllSangayFeatures(true); toc;
cd ~/igdata; %~/research/now/sangay/
maxDisp = 150;
writeFlag = false;
load('~/masa/old/research/now/sangay/sangayGentleBoostCompact');
load ~/igdata/AllFeatures.mat;

%% orig
[YfitAll,t,eqmag,ncc,energyMag,merr,tFeatures,NCC,Neff,Msangay] = ...
    applyClassificationModel(cGentleBoostEnsemble,features,tFeatures,2,-1.8);

% [YfitAll,t,eqmag,ncc,energyMag,merr,tFeatures,NCC,Neff,Msangay] = ...
%     applyClassificationModel(cGentleBoostEnsemble,features,tFeatures,3,1.5);
neff = Neff(YfitAll);

%%
goodI = neff >= 3 & merr < 1 & t <= datetime(2022,05,01) & eqmag >= 1.9 & eqmag <= 2.5;

%%
% sangay_kstnms = ["BPAT";"BMAS";"BRTU";"BULB";"BBIL";"BRUN";"TAMH";"PORT";...
%     "PKYU";"TAIS";"PIS1";"PUYO";"PIAT";"SAGA";"PORT"];
sangay_kstnms = ["PORT";"SAGA";"PUYO"];
lK = length(sangay_kstnms);
sangay_channels = ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"];
for i = 1:lK
    kstnm_ = sangay_kstnms(i);
    S = extractWaveforms(t(goodI),120,kstnm_,sangay_channels,"EC","",...
        false,true,1,true,[0.6 1.2 false false true 10]);
    badI = sum(~isnat(pull(S,"ref")),2)<3;
    S(badI,:) = [];
    S = S(:,1:3);
    cd ~/masa/subspace_detector/sangay/
    save(kstnm_,"S");
end