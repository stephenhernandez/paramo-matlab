clear; close all; clc;

cd ~/research/now/reventador/
% [tabs,z2p,NCC,Neff,templateIndex,scaledCC,kurt] = ...
%     filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v5',10,10);
[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v6',10); %,10);

%%
minKurt = min(kurt(:,1:2),[],2,'omitnan');
maxKurt = max(kurt(:,1:2),[],2,'omitnan');
kurtRatio = maxKurt./minKurt;

z2p = z2p(:,1); %max(z2p,[],2);
kurt = kurt(:,1); %max(kurt,[],2);

%%
minAmp = 400;
winlen = 50;
nboot = 500;
zeroPhaseFlag = true;

%%
goodI1 = z2p <= minAmp & ~(Neff == 2 & NCC < 0.035) & ~(Neff == 1 & scaledCC < 50) & kurt < 25;
goodI2 = z2p >= minAmp & ~(Neff == 2 & NCC < 0.035) & ~(Neff == 1 & scaledCC < 50) & kurt < 25;

%%
tGood1 = tabs(goodI1);
aGood1 = z2p(goodI1);
[acut,startIndex1,endIndex1] = cutWindows(aGood1,winlen,winlen-1,false);
tcut1 = cutWindows(datenum(tGood1),winlen,winlen-1,false);
tcut1 = dn2dt(tcut1);
difft1 = seconds(diff(tcut1));

%%
disp('getting variances, please be patient');
tic;
stepsize = max([1 ceil(winlen/20)]);
varIndices = (1:stepsize:size(acut,2))';
p = [2.5 97.5];
ampBars = medboot(acut(:,varIndices),nboot,p);
rateBars = medboot(difft1(:,varIndices),nboot,p);
toc;

%%
rollingBestimate = log10(exp(1))./(nanmedian(log10(acut)) - log10(minAmp) + 0.005);
rollingBestimate = rollingBestimate';
rollingBestimate = zpkFilter(rollingBestimate,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
medianAmplitude = nanmedian(acut)';
medianAmplitude = zpkFilter(medianAmplitude,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
rollingRate1 = 86400./nanmedian(difft1);
rollingRate1 = rollingRate1';
rollingRate1 = zpkFilter(rollingRate1,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
tGood2 = tabs(goodI2);
aGood2 = z2p(goodI2);
[acut2,startIndex2,endIndex2] = cutWindows(aGood2,winlen,winlen-1,false);
tcut2 = cutWindows(datenum(tGood2),winlen,winlen-1,false);
tcut2 = dn2dt(tcut2);
difft2 = seconds(diff(tcut2));

%%
disp('getting variances, please be patient');
tic;
stepsize = max([1 ceil(winlen/20)]);
varIndices = (1:stepsize:size(acut2,2))';
p = [2.5 97.5];
ampBars2 = medboot(acut2(:,varIndices),nboot,p);
rateBars2 = medboot(difft2(:,varIndices),nboot,p);
toc;

%%
rollingBestimate2 = log10(exp(1))./(nanmedian(log10(acut2)) - log10(minAmp) + 0.005);
rollingBestimate2 = rollingBestimate2';
rollingBestimate2 = zpkFilter(rollingBestimate2,-inf,1/winlen,1,1,0);

%%
medianAmplitude2 = nanmedian(acut2)';
medianAmplitude2 = zpkFilter(medianAmplitude2,-inf,1/winlen,1,1,0);

%%
rollingRate2 = 86400./nanmedian(difft2);
rollingRate2 = rollingRate2';
rollingRate2 = zpkFilter(rollingRate2,-inf,1/winlen,1,1,0);

%%
figure('units','normalized','outerposition',[0 0 0.8 1]);
ax = subplot(211);
pp = semilogy(ax,tGood1(endIndex1),rollingRate1,'linewidth',2);
zoom on; grid on;
pp.Color(4) = 0.65;
hold on;
pp2 = semilogy(ax,tGood2(endIndex2),rollingRate2,'linewidth',2);
zoom on; grid on;
pp2.Color(4) = 0.65;
legend('little rate','big rate');

%%
ax(2) = subplot(212);
t2 = sort([tGood1(endIndex1);tGood2(endIndex2)]);
rollingRate3 = interp1(tGood1(endIndex1),rollingRate1,t2,'pchip');
rollingRate4 = interp1(tGood2(endIndex2),rollingRate2,t2,'pchip');
pp3 = semilogy(ax(2),t2,rollingRate4./rollingRate3,'linewidth',2);
zoom on; grid on;
pp3.Color(4) = 0.65;
title('big rate / little rate');
hold on;
plot([t2(1) t2(end)],[1 1],'-','linewidth',1,'color',[0.5 0.5 0.5]);
linkaxes(ax,'x');
