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
figure('units','normalized','outerposition',[0 0 1 1]);
ss = scatter(t,ncc,3*exp(eqmag),log10(merr),'filled'); zoom on; colorbar;
ss.MarkerEdgeColor = 'w'; ss.MarkerFaceAlpha = 0.6; ss.MarkerEdgeAlpha = 0.2;
clim([-2 0]);
grid on;

%%
TTtmp = timetable(t(end-maxDisp+1:end),t2r(t(end-maxDisp+1:end),hours(24)),...
    eqmag(end-maxDisp+1:end),ncc(end-maxDisp+1:end),neff(end-maxDisp+1:end),...
    'VariableNames',{'rate','magnitude','gamma','neff'});

disp(TTtmp);
sum(YfitAll)
Nmed = 5;
minMag = 1.9;
magI = eqmag >= minMag;
rate = t2r(t,hours(24));

nDays = 02;
minN_perDay = 04;
minN = 1+nDays*minN_perDay;

highEnergy = medfiltSH(energyMag(magI),Nmed);
highEnergy = interp1(datenum(t(magI)),cumsum(highEnergy),datenum(t),'nearest','extrap');
lowEnergy = medfiltSH(energyMag(~magI),Nmed);
lowEnergy = interp1(datenum(t(~magI)),cumsum(lowEnergy),datenum(t),'nearest','extrap');

%%
rateLong = t2r(t,hours(24*nDays));
rateLowOrig = t2r(t(eqmag < minMag),hours(24*nDays));
rateLow = rateLowOrig/nDays;

rateHighOrig = t2r(t(eqmag >= minMag),hours(24*nDays));
rateHigh = rateHighOrig/nDays;

rateLow = interp1(datenum(t(eqmag < minMag)),rateLow,datenum(t),'nearest','extrap');
rateHigh = interp1(datenum(t(eqmag >= minMag)),rateHigh,datenum(t),'nearest','extrap');

goodI = rateHigh >= 1/nDays & rateLow >= 1/nDays & rateLow+rateHigh >= minN_perDay;

%%
[N,edges] = histcounts(rateLong,[unique(rateLong); inf]);
edges = edges(1:end-1)';
N = N';
edgesI = edges >= minN;
edges2 = edges(edgesI);
N2 = N(edgesI);
energyRatio = NaN(size(rateLong));
highEnergyRate = energyRatio;
lowEnergyRate = energyRatio;

for i = 1:length(edges2)
    q = edges2(i);
    qI = find(rateLong == q);
    qI2 = qI - q + 1;
    
    %%
    highEnergyRate_ = (highEnergy(qI) - highEnergy(qI2));
    highEnergyRate(qI) = highEnergyRate_/nDays;
    lowEnergyRate_ = (lowEnergy(qI) - lowEnergy(qI2));
    lowEnergyRate(qI) = lowEnergyRate_/nDays;
    energyRatio(qI) = highEnergyRate_./lowEnergyRate_;
end

%%
figure('units','normalized','outerposition',[0 0 1 1]);
axC(1) = subplot(211);
plot(t,1:length(t),'.');
zoom on; grid on;

axC(2) = subplot(212);
plot(t,cumsum(medfiltSH(energyMag,Nmed)),'.');
zoom on; grid on;
linkaxes(axC,'x');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax = subplot(211);
plot(t,rate,'.'); zoom on; grid on;
title('Daily rate');
[N,edges] = histcounts(t,dateshift(min(t),'start','day'):dateshift(max(t),'end','day')+1);
edges = edges(1:end-1)';
N = N';
hold on; stairs(edges,N,'k-'); zoom on;

ax(2) = subplot(212);
plot(t,eqmag,'.');
zoom on; grid on; title('$M_{sangay}$'); linkaxes(ax,'x');
%sum(eqmag >= 3.15)

%%
figure('units','normalized','outerposition',[0 0 1 1]);
clear ax__;
ax__(1) = subplot(311);
semilogy(ax__(1),t(goodI),rateHigh(goodI)./rateLow(goodI),'.'); zoom on; grid on;
hold on;
semilogy(ax__(1),[min(t) max(t)],[1 1],'--','linewidth',1,'Color',[0 0 0]); zoom on; grid on;
semilogy(ax__(1),[min(t) max(t)],[10 10],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on; grid on;
semilogy(ax__(1),[min(t) max(t)],[0.1 0.1],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on; grid on;
title('large rate / small rate');

%
ax__(2) = subplot(312);
goodI = goodI & isfinite(energyRatio);
semilogy(ax__(2),t(goodI),energyRatio(goodI),'.'); zoom on; grid on;
hold on;

semilogy(ax__(2),[min(t) max(t)],[1 1],'--','linewidth',1,'Color',[0 0 0]); zoom on; grid on;
semilogy(ax__(2),[min(t) max(t)],[10 10],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on; grid on;
semilogy(ax__(2),[min(t) max(t)],[0.1 0.1],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on; grid on;
title('large energy / small energy');

ax__(3) = subplot(313);
plot(ax__(3),t,rateLong/nDays,'.'); zoom on; grid on;
linkaxes(ax__,'x');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
clear ax__;
ax___(1) = subplot(211);
plot(ax___(1),t,highEnergyRate,'.'); zoom on; grid on;
title('daily energy rate, large events');

%
ax___(2) = subplot(212);
plot(ax___(2),t,lowEnergyRate,'.'); zoom on; grid on;
title('daily energy rate, small events');
linkaxes(ax___,'xy');

%%
figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t,(highEnergyRate+lowEnergyRate)./(rateHigh+rateLow),'.'); zoom on; grid on;
hold on;
semilogy(t,medfiltSH((highEnergyRate+lowEnergyRate)./(rateHigh+rateLow),Nmed),'-','linewidth',2); zoom on; grid on;

%%
if writeFlag
    return;
end
nowTime = dn2dt(now)+hours(5);
[yyyy,mm,dd] = datevec(nowTime);
dayStr = num2str(dd);
if dd < 10
    dayStr = ['0',dayStr];
end
fileName = strcat("SangayRegionalCatalog_UPDATED");
outFile = '~/test.txt';
tOrig = t;
t = datenum(t) - 693960;
m = eqmag; NN = Neff(YfitAll);
formatSpec = '%f %f %u %f %f %u';

str = compose(formatSpec,t,m,NN,ncc,energyMag,rate);
str = string(str);
fileID = fopen(outFile,'w');
fprintf(fileID,'%s\n',str);
fclose(fileID);

%%
figure('units','normalized','outerposition',[0 0 1/2 1]);
plot(tOrig,cumsum(medfiltSH(energyMag,Nmed)),'.');
zoom on; grid on;
% 
% figure(); 
% ntest = 35; plot(tOrig(2:end),3600*(10.^medfiltSH(m(2:end),ntest))./medfiltSH(seconds(diff(tOrig)),ntest),'.'); zoom on;

%%
% % close all
% % winlen = 128; nOverlap = winlen - 8; nboot = 400;
% % dcut = cutWindows(eqmag,winlen,nOverlap,false);
% % dcut = demean(dcut);
% % dcut = dcut ./ rssq(dcut);
% % dcut = sort(dcut);
% % dcut = detrend(dcut);
% % dcutOrig = dcut ./ rssq(dcut);
% % BOOTI = randi(winlen,[winlen*nboot,1]);
% % dcut = cutWindows(eqmag,winlen,nOverlap,false);
% % tic; for i = 1:size(dcutOrig,2)
% %     dcut_ = dcut(:,i);
% %     dcut_ = dcut_(BOOTI);
% %     dcut_ = reshape(dcut_,[winlen nboot]);
% %     dcut_ = demean(dcut_);
% %     dcut_ = dcut_ ./ rssq(dcut_);
% %     dcut_ = sort(dcut_);
% %     dcut_ = detrend(dcut_);
% %     dcut_ = dcut_ ./ rssq(dcut_);
% %     dcut(:,i) = mean(dcut_,2);
% % end; toc;
% % [maxccp,plags,maxccn,nlags] = doccFreqCircShift(dcut,true);
% % [newFamilies,l_uniq_indices,Nsingletons,tree] = generateFamilies(maxccp,0.5,'weighted');
% % close all; clear tmp_stack_; for i = 1:min([35 length(newFamilies)]); fam_ = newFamilies{i}; data = dcut(:,fam_); tmp_stack_(:,i) = plot_family(data,(1:size(data,2))',20,1); end;
% % %%
% % figure(); plot(tmp_stack_,'linewidth',3); zoom on; grid on;
% % figure(); plot(l_uniq_indices,'.'); zoom on; grid on;
% % dcut = cutWindows(eqmag,winlen,winlen-1,false);
% % dcut = demean(dcut);
% % dcut = dcut ./ rssq(dcut);
% % dcut = sort(dcut);
% % dcut = detrend(dcut);
% % dcutOrig = dcut ./ rssq(dcut);
% % %%
% % t2 = tOrig(winlen:end); mag2 = eqmag(winlen:end);
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % cc = NaN(size(tmp_stack_,2),size(dcutOrig,2));
% % for i = 1:size(tmp_stack_,2)
% %     thisStack = tmp_stack_(:,i);
% %     thisStack = repmat(thisStack,[1 size(dcutOrig,2)]);
% %     cc_ = sum(thisStack.*dcutOrig)';
% %     cc(i,:) = cc_';
% %     plot(t2,cc_,'.'); hold on;
% % end
% % zoom on;
% % %%
% % [~,maxI] = max(cc);
% % uniqueIndex = unique(maxI);
% % %%
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % for i = 1:length(uniqueIndex)
% %     ui_ = uniqueIndex(i);
% %     UI = ui_ == maxI;
% %     plot(t2(UI),mag2(UI),'.'); hold on; grid on; zoom on;
% % end;
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % for i = 1:length(uniqueIndex)
% %     ui_ = uniqueIndex(i);
% %     UI = ui_ == maxI;
% %     plot(t2(UI),(1:sum(UI))','.'); hold on; grid on; zoom on;
% % end;
% % [U,S,V] = svd(dcutOrig,"econ");
% % figure('units','normalized','outerposition',[0 0 1 1]);
% % for i = 1:length(uniqueIndex)
% %     ui_ = uniqueIndex(i);
% %     UI = ui_ == maxI;
% %     fam_ = dcutOrig(:,UI);
% %     fam_ = mean(fam_,2);
% %     plot((fam_),'linewidth',3); zoom on; hold on; grid on;
% %     disp(sum(UI));
% % end