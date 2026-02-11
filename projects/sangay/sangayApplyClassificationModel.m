clear; close all;
cd ~/research/now/sangay/

%%
tic;
load('sangayGentleBoostCompact');
load AllFeatures.mat;
toc;

%%
tic;
load('sangayRandomForest','timeMaster','labelMaster');
bonafideBad = timeMaster(~labelMaster);
probablyGood = timeMaster(labelMaster);

clear labels timeMaster;
toc;

%%
dists = [66700,63289,102391,71016,57654,56096,91971,76817,65576];
sensitivities = [3.141950e+08 4.872110e+08 3.141950e+08 3.141950e+08...
    4.872110e+08 4.872110e+08 3.141950e+08 3.141950e+08 4.872110e+08];%./dists;

%%
predGood = true(119,1);
features = features(:,predGood);

%%
[Msangay,Err_sangay,mlv,wa] = sangayMagnitudeCalculation(tFeatures,features(:,3:11));
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
disp(size(features));

%%
%[tabs,~,NCC,Neff,p2rms,kurt] = filterUniqueEvents('~/research/now/sangay/SangayRegionalAnalysis_v6');
tabs = tFeatures;%    tabs(1:size(features,1));
NCC = features(:,1); %NCC(1:size(features,1));
%z2p = z2p(1:size(features,1),:);
Neff = features(:,2); %Neff(1:size(features,1));

%%
% predict new classes (takes ~55 seconds)
tic; 
[YfitAll,score] = predict(cGentleBoostEnsemble,features); 
toc;

%%
dumRef = min([dateshift(min(bonafideBad),'start','day') dateshift(min(tabs),'start','day')]);
[lia,locb] = ismembertol(seconds(bonafideBad-dumRef),seconds(tabs-dumRef),1,'DataScale',1);
locb = locb(lia);
falsePositivesFromSAGA = YfitAll(locb); %in theory, these should ALL be falses

%%
dumRef = min([dateshift(min(probablyGood),'start','day') dateshift(min(tabs),'start','day')]);
[lia,locb] = ismembertol(seconds(probablyGood-dumRef),seconds(tabs-dumRef),1,'DataScale',1);
locb = locb(lia);
falseNegativeFromSAGA = ~YfitAll(locb); %in theory, these should ALL be falses

%%
YfitAll(falsePositivesFromSAGA) = false;
YfitAll(falseNegativeFromSAGA) = true;

%%
YfitAll = YfitAll & Neff >= 3; % & Msangay < 3.2; % | NCC >= 0.5 | (Neff >= 6 & NCC >= 0.25);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax_(1) = subplot(311); plot(tabs(~YfitAll),Msangay(~YfitAll),'.'); zoom on; grid on; hold on;
plot(tabs(YfitAll),Msangay(YfitAll),'o'); grid on; hold on;

ax_(2) = subplot(312); plot(tabs(~YfitAll),NCC(~YfitAll),'.');
hold on; zoom on; grid on;
plot(tabs(YfitAll),NCC(YfitAll),'o'); 

ax_(3) = subplot(313); plot(tabs(~YfitAll),Neff(~YfitAll),'.');
hold on; grid on;
plot(tabs(YfitAll),Neff(YfitAll),'o');

linkaxes(ax_,'x');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(Err_sangay(YfitAll),Msangay(YfitAll),'.');
zoom on; grid on;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs(YfitAll),1:sum(YfitAll),'.'); 
zoom on; grid on;

%%
[N,edges] = histcounts(tabs(YfitAll),dateshift(min(tabs(YfitAll | ~YfitAll)),'start','day'):hours(8):dateshift(max(tabs(YfitAll | ~YfitAll)),'end','day')+1);
N1 = N';
edges = edges(1:end-1)';
figure('units','normalized','outerposition',[0 0 1 1]);
axsp(1) = subplot(211);
stairs(edges,N1,'linewidth',2); zoom on; grid on;
hold on;

% bar chart for ALL events
[N,edges] = histcounts(tabs(YfitAll | ~YfitAll),dateshift(min(tabs(YfitAll | ~YfitAll)),'start','day'):hours(8):dateshift(max(tabs(YfitAll | ~YfitAll)),'end','day')+1);
N = N';
edges = edges(1:end-1)';
%figure('units','normalized','outerposition',[0 0 1 1]);
pp = stairs(edges,N,'linewidth',0.01); zoom on; grid on;
pp.Color(4) = 0.7;

axsp(2) = subplot(212);
fractionKept = N1./N;
plot(edges,fractionKept,'.'); grid on; zoom on;
linkaxes(axsp,'x');
fprintf('mean fraction kept: <strong>%s</strong>\n',num2str(mean(fractionKept(edges >= datetime(2020,02,01) & N1 > 3))));
fprintf('median fraction kept: <strong>%s</strong>\n',num2str(median(fractionKept(edges >= datetime(2020,02,01) & N1 > 3))));
fprintf('sum total in modern era: <strong>%s</strong>\n',num2str(sum(N1(edges >= datetime(2020,02,01)))));

%
disp(' ');
fprintf('%s\n',num2str(sum(YfitAll)));

%%
N = 7;
tgood = tabs(YfitAll);
difftime = seconds(diff(tgood));
[tcut,startIndex,endIndex,badFlag] = cutWindows(datenum(tgood),N,N-1,false);
iet = seconds(diff(dn2dt(tcut)));
tEnd = tgood(endIndex);
tStart = tgood(startIndex);
nccGood = NCC(YfitAll);

mm = Msangay(YfitAll);
medMag = median(cutWindows(mm,N,N-1,false))';
energyMag = (1.5.*mm+4.8);

medEventsPerDay = 86400./median(iet)';
meanEventsPerDay = N*86400./seconds(tEnd - tStart);

medEnergyPerEvent = (10.^(1.5.*medMag+4.8));
meanEnergyPerEvent = (10.^(1.5.*mean(cutWindows(mm,N,N-1,false))'+4.8));
medEnergyPerDay = medEnergyPerEvent.*medEventsPerDay;

%
nMinutes = 5;
smoothEnergyRate = zpkFilter(medEnergyPerDay,-inf,1/4,1,1,0); %(N./days(tgood(endIndex) - tgood(startIndex)));

t2 = datenum(dateshift(tEnd(3),'start','minute'):minutes(nMinutes):dateshift(tEnd(end-2),'end','minute'))';
smoothEnergyRate2 = interp1(datenum(tEnd),smoothEnergyRate,t2,'next');
t2 = dn2dt(t2);

%
[~,energyRatePeaksI] = findpeaks(smoothEnergyRate2,'MINPEAKDISTANCE',1440/nMinutes/2);
[~,energyRateTroughsI] = findpeaks(-smoothEnergyRate2,'MINPEAKDISTANCE',1440/nMinutes/2);

% figure 5
figure('units','normalized','outerposition',[0 0 1 1]);
ax2(1) = subplot(311);
% plot(tabs(~YfitAll),Msangay(~YfitAll),'.');
% zoom on; grid on; hold on;
% plot(tabs(YfitAll),Msangay(YfitAll),'o'); grid on;
% zoom on; grid on; hold on;
semilogy(tgood,10.^energyMag,'.'); zoom on; grid on; hold on;
semilogy(tEnd,medEnergyPerEvent,'.');
%semilogy(tEnd,meanEnergyPerEvent,'.');
title('Energy [J per event]');

ax2(2) = subplot(312);
plot(tEnd,medEventsPerDay,'.');
zoom on; grid on; hold on;
%semilogy(tEnd,meanEventsPerDay,'.');
title('Event Rate [Events per day]');

ax2(3) = subplot(313);
semilogy(tEnd,medEnergyPerDay,'.');
zoom on; grid on; hold on;
%semilogy(tEnd,meanEnergyPerEvent.*meanEventsPerDay,'.');
title('Energy Rate [J per day]');
%hold on;
%semilogy(tEnd,cumsum(smoothEnergyRate).*days(tEnd(end) - tEnd(1)),'.'); %cumulative energy

%
smoothEnergyRate3 = interp1(datenum(tEnd),cumsum(smoothEnergyRate).*sum(days(diff(tEnd))),datenum(t2)')';
smoothEnergyRate4 = zpkFilter(diff(smoothEnergyRate3)./sum(days(diff(tEnd))),-inf,1/4,1,1,0);
t3 = t2(1:end-1);

%semilogy(t3,smoothEnergyRate4,'.');
%semilogy(t2,smoothEnergyRate3,'.');
%semilogy(t3,cumsum(smoothEnergyRate4).*sum(days(diff(t3))),'o');
zoom on; grid on;

linkaxes(ax2,'x');

%%
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(tEnd,cumsum(medEnergyPerEvent),'.'); zoom on; grid on;

%%
%close all; 
%close(7); close(8); %close(9);
medEnergyPerEvent2 = (10.^(1.5.*medfiltSH(mm,N)+4.8));

tAll = [tgood; tgood+hours(6)];
iAll = [ones(size(tgood)); zeros(size(tgood))];
eAll = [zeros(size(tgood)); medEnergyPerEvent2];
[tAll,sI] = sort(tAll);
iAll = iAll(sI);
eAll2 = eAll(sI);
eAll = cumsum(eAll2);

%
posI = find(iAll == 1);
negI = find(iAll == 0);
td = tAll(posI);

cum1 = cumsum(iAll);
cum2 = cumsum(~iAll);
newSum2 = cum1 - cum2;
energyBin = eAll(negI)-eAll(posI);
rate = newSum2(posI);
accI = posI - rate;
accI(~accI) = 1;
dr = newSum2(posI) - newSum2(accI);

acc2I = (1:length(rate))' - rate;
acc2I(acc2I < 1) = 1; %acc2I should be applicable to accI as well (a copy) but i havent tested if true
de = energyBin - energyBin(acc2I);
%dr = rate - rate(acc2I);

smoothRate = zpkFilter(rate,-inf,1/20,1,1,false);

figure(); 
axr(1) = subplot(211);
plot(td,rate,'.'); zoom on; grid on;

axr(2) = subplot(212);
plot([min(td) max(td)],[0 0],'-','linewidth',2,'Color',[0.5 0.5 0.5]);
hold on;
hhacc_rate = plot(td,dr,'.');
plot(td,medfiltSH(dr,1,true),'-','linewidth',1); 
zoom on; grid on;
title('event acceleration');
linkaxes(axr,'x');

figure();
axe(1) = subplot(211);
plot(td,energyBin,'.');
zoom on; grid on;

axe(2) = subplot(212);
plot([min(td) max(td)],[0 0],'-','linewidth',2,'Color',[0.5 0.5 0.5]);
hold on;
plot(td,de,'.'); 
plot(td,medfiltSH(de,1,true),'-','linewidth',1); 
%plot(td(2:end),medfiltSH(86400*diff(energyBin)./seconds(diff(td)),9),'.'); 
zoom on; grid on;
title('energy acceleration');
linkaxes(axe,'x');

%%
figure(); 
axx(1) = subplot(311); plot(tgood,rate,'.'); grid on; 
axx(2) = subplot(312); plot(tgood,mm,'o'); grid on; 
axx(3) = subplot(313); plot(tgood,nccGood,'o'); zoom on; grid on; 

linkaxes(axx,'x'); 
xlim([datetime(2021,02,01) datetime(2021,04,06)]); 
axx(1).Title.String = 'Events in previous 6 hours'; 
axx(2).Title.String = 'Estimated Magnitude'; 
axx(3).Title.String = '$\gamma$ (Proportional to Confidence)';

%% figure 6
% imp = predictorImportance(cGentleBoostEnsemble);
% figure();
% bar(imp);
% title('Predictor Importance Estimates');
% ylabel('Estimates');
% xlabel('Predictors');
% zoom on;
% 
% % figure 7
% figure(); 
% plot(1:length(imp),cumsum(sort(imp,'descend'))/sum(imp),'.'); 
% zoom on; grid on;

% figure