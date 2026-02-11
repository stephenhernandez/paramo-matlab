clear; close all; clc;
%make sure `variableNames' match the names of output variables
variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";"templateNumber";"madMain";"nUsedMain"};
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    loadRepeaterCatalog("ggp",variableNames);

%%
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    filterCatalog(tMain,ccMain,10,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain);

%%
tic;
ccMainThresh = 0.16;
madThresh = 15;
ampThresh = 1e3;
%ccI = templateNumber == 1 & madMain >= madThresh & ccMain >= ccMainThresh & nUsedMain >= 14 & tMain >= datetime(2022,09,01) & ampMain >= 1e1;
ccI = (templateNumber == 2 & madMain >= 18 & ccMain >= 0.18 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= ampThresh) | ...
    (templateNumber == 1 & madMain >= 10 & ccMain >= 0.165 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= ampThresh);

%ccI = templateNumber == 1 & madMain >= 10 & ccMain >= 0.1 & nUsedMain > 2 & tMain >= datetime(2012,07,01) & ampMain >= 1e0 & magMain >= 4;

% ccI = templateNumber == 2 & madMain >= 18 & ccMain >= 0.18 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= 5e2; % for template 2
%ccI = templateNumber == 1 & madMain >= 10 & ccMain >= 0.16 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= 5e2; % for template 1
%ccI = templateNumber == 1 & madMain >= 18 & ccMain >= 0.12 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= 5e2;
%ccI = templateNumber == 2 & madMain >= 20 & ccMain >= 0.12 & nUsedMain > 2 & tMain >= datetime(1997,01,01) & ampMain >= 5e2;
% & ((ccMain >= 0.1 & madMain >= 10 & nUsedMain > 3) | (ccMain >= 0.25 & madMain >= 11 & nUsedMain == 3) | (ccMain >= 0.35 & madMain >= 14 & nUsedMain < 3));

close all; clear ax;
% Figure 1
figure();
ax(1) = subplot(311); semilogy(tMain(ccI),ccMain(ccI),'.');
zoom on; grid on; hold on;
semilogy(tMain(~ccI),ccMain(~ccI),'.'); zoom on; grid on;

ax(2) = subplot(312);
semilogy(tMain(ccI),madMain(ccI),'.'); zoom on; grid on; hold on;
semilogy(tMain(~ccI),madMain(~ccI),'.');

ax(3) = subplot(313);
semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on; hold on;
semilogy(tMain(~ccI),ampMain(~ccI),'.');
linkaxes(ax,'x');

%%
% Figure 2
figure();
plot(tMain(ccI),1:sum(ccI),'.'); zoom on; grid on;

% Figure 3
figure();
semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on;

% Figure 4
figure();
plot(tMain(ccI),nUsedMain(ccI),'.'); zoom on; grid on;

% % Figure 5
% figure();
% plot(tMain(ccI& nUsedMain == 14),1:sum(ccI&nUsedMain==14),'.'); zoom on;

%
nDayVec = unique([1./(1:7)'; (1:12)'; 24]);
rateRatio = [];

for i = 1:length(nDayVec)
    nDays = 1./nDayVec(i);
    rate2 = t2r(tMain(ccI),days(nDays),[],false)/nDays;
    rate1 = t2r(tMain(ccI),days(nDays),[],true)/nDays;
    rateRatio = [rateRatio rate2./rate1];
end

Nfilt = 3;
w = nDayVec/sum(nDayVec);
rateRatio = convFilter(rateRatio,3,true);

% stackRate = sum(rateRatio.*w',2,"omitnan");
% [cc_,locs_,width_,prom_] = ...
%     findpeaks(stackRate,'MINPEAKDISTANCE',5,'MINPEAKHEIGHT',1.5,'MinPeakProminence',0.5);

stackRate1 = mean(rateRatio,2,"omitnan");

%
%close all;
tMain2 = tMain(ccI);
%tMain2 = tdum; stackRate1 = stackRate2;
[cc_,locs_,width_,prom_] = ...
    findpeaks(stackRate1,'MINPEAKDISTANCE',5,'MINPEAKHEIGHT',2,'MinPeakProminence',1);

%stackRate = medfiltSH(mean(rateRatio,2,"omitnan"),Nfilt,true);
%stackRate = medfiltSH(sum((rateRatio).*w',2,"omitnan"),3,true);
%stackRate = zpkFilter(stackRate,-inf,1/Nfilt,1,1,1);
%box = ones(Nfilt,1)/Nfilt; stackRate = mean([fftfilt(box,stackRate) flipud(fftfilt(box,flipud(stackRate)))],2);
%stackRate = zpkFilter(sum((rateRatio).*w,2,"omitnan"),-inf,1/3,1,1,1);

% Figure 5
figure();
axStack(1) = subplot(211);
semilogy(tMain2,stackRate1,'.'); zoom on; grid on; hold on;

% [cc2_,locs2_] = findpeaks(stackRate,'MINPEAKDISTANCE',40,'MINPEAKHEIGHT',2);
% cc_ = [cc_; cc2_];
% locs_ = [locs_; locs2_];
% peakVec = unique([cc_ locs_],'rows');
% cc_ = peakVec(:,1);
% locs_ = peakVec(:,2);

ccI2 = find(ccI);
plot(tMain2(locs_),cc_,'p');

% Figure 7
% figure();
% semilogy(tMain(ccI),rateRatio,'.'); zoom on; grid on; hold on;
% semilogy(tMain(ccI),stackRate,'ko','linewidth',2); zoom on; grid on; hold on;

% Figure 6  %use this to find when the swarm ends
figure();
semilogy(tMain2,1./stackRate1,'.');
zoom on; grid on; hold on;

% Figure 7
tMain3 = tMain2(locs_);
tMain4 = tMain3(find(tMain3 >= datetime(2008,01,01),1):end);
figure();
interswarmReposeIntervalsOrig = days(diff(tMain4));
semilogy(tMain4(2:end),interswarmReposeIntervalsOrig,'.-');
zoom on; grid on;
hold on;
semilogy(tMain4(2:end),medfiltSH(interswarmReposeIntervalsOrig,11,true),'-','linewidth',4);
toc;

figure();
ax(1) = subplot(121);
semilogy(sort(ampMain(ccI)),1-((0:sum(ccI)-1)'/sum(ccI)),'.'); zoom on; grid on;
ax(2) = subplot(122);
loglog(sort(ampMain(ccI)),1-((0:sum(ccI)-1)'/sum(ccI)),'.'); zoom on; grid on;

%%
% close all;
% tMain2 = tMain(ccI);
rateRatio = [];
dailyRate1 = [];
%nWinVec = (10:10:200)';
%nWinVec = round(logspace(log10(2),log10(200),51))';
nWinVec = (2:6)';
%nWinVec = (3:8:101)';
maxNwin = max(nWinVec);
t3 = seconds(tMain2 - tMain2(1));

for i = 1:length(nWinVec)
    nWin = nWinVec(i);
    G = [1; zeros(nWin-2,1); -1];
    t4 = fftfilt(G,t3);

    t5 = t4(nWin:end-nWin)./t4(2*nWin:end);
    rate_ = t5(maxNwin-nWin+1:end-maxNwin+nWin);
    rateRatio = [rateRatio rate_];

    %dRate = nWin.*86400./t4(nWin:end-nWin);
    dRate = nWin.*86400./t4(2*nWin:end);
    dRate_ = dRate(maxNwin-nWin+1:end-maxNwin+nWin);
    dailyRate1 = [dailyRate1 dRate_];
end
w = nWinVec./sum(nWinVec); %(1./nWinVec)/sum(1./nWinVec);
tdum = tMain2(maxNwin:end-maxNwin);

rateRatio = convFilter(rateRatio,3,true);
%stackRate = mean(rateRatio,2,"omitnan");
stackRate1 = sum((rateRatio).*w',2,"omitnan"); %weighted average
stackRate2 = median(rateRatio,2,"omitnan"); %average rate ratio (snr) using multiple windows

% figure();
% ax(1) = subplot(211);
% semilogy(tdum,1./stackRate1,'.','linewidth',2); zoom on; grid on; hold on;
% ax(2) = subplot(212);
% semilogy(tdum,1./stackRate2,'.','linewidth',2); zoom on; grid on; %legend('weighted stack','mean');
% linkaxes(ax,'x');

figure();
ax(1) = subplot(211);
semilogy(tdum,stackRate1,'.','linewidth',2); zoom on; grid on; hold on;
ax(2) = subplot(212);
semilogy(tdum,stackRate2,'.','linewidth',2); zoom on; grid on; %legend('weighted stack','mean');
linkaxes(ax,'x');

dailyRate1 = convFilter(dailyRate1,3,true);
rateWA = sum((dailyRate1).*w',2,"omitnan"); %weighted average
rateMedian = median(dailyRate1,2,"omitnan"); %average rate ratio (snr) using multiple windows

figure();
ax(1) = subplot(211);
semilogy(tdum,rateWA,'.','linewidth',2); zoom on; grid on; hold on;
ax(2) = subplot(212);
semilogy(tdum,rateMedian,'.','linewidth',2); zoom on; grid on; %legend('weighted stack','mean');
linkaxes(ax,'x');


%% alterntive method to get snr, using fixed time windows instead of fixed N
snr = [];
%winVec = [0.25; 0.5; 1; 2; 4; 8]; %fixed time windows, in days
%winVec = [0.25; 0.5; 1; 2; 3; 4]; %fixed time windows, in days
winVec = [1./(2:4)'; (1:4)']; %fixed time windows, in days
%winVec = [1/8; 1/6; 1/4; 1/3; 1/2; 1; 2; 3; 4; 6; 8]; %fixed time windows, in days
for i = 1:length(winVec)
    nDays = 1./winVec(i);
    rate = t2r(tMain2,days(nDays),[],true);
    rate2 = t2r(tMain2,days(nDays),[],~true);
    snr = [snr rate2./rate];
end

w = winVec;
w = w./sum(w);
meanSnr = sum(snr.*w',2,"omitnan"); %meanSnr = 10.^(meanSnr/20);
%meanSnr = mean(snr,2); %meanSnr = 10.^(meanSnr/20);

[cc_,locs_,width_,prom_] = ...
    findpeaks(meanSnr,'MINPEAKDISTANCE',5,'MINPEAKHEIGHT',3,'MinPeakProminence',1); %changed the prominence to 0.5, is different to fixed-N example above
%[cc_,locs_,width_,prom_] = ...
%    findpeaks(meanSnr,'MINPEAKDISTANCE',5,'MINPEAKHEIGHT',1,'MinPeakProminence',1); %changed the prominence to 0.5, is different to fixed-N example above


% Figure 5
figure(5);
axStack(2) = subplot(212);
semilogy(tMain2,meanSnr,'.'); zoom on; grid on; hold on;

ccI2 = find(ccI);
plot(tMain2(locs_),cc_,'p');
title('weighted');
linkaxes(axStack,'x');

% Figure 11
tMain3 = tMain2(locs_);
figure();
semilogy(tMain3(2:end),days(diff(tMain3)),'.-');
zoom on; grid on;
hold on;
semilogy(tMain3(2:end),medfiltSH(days(diff(tMain3)),11,true),'-','linewidth',4);
toc;

tMain4 = tMain3(find(tMain3 >= datetime(2008,01,01),1):end);
interswarmReposeIntervals = days(diff(tMain4));
figure();
semilogy(tMain4(2:end),interswarmReposeIntervals,'.-');
zoom on; grid on;
hold on;
semilogy(tMain4(2:end),medfiltSH(interswarmReposeIntervals,11,true),'-','linewidth',4);
toc;

figure();
plot(sort(interswarmReposeIntervalsOrig),1-((0:length(interswarmReposeIntervalsOrig)-1)'/length(interswarmReposeIntervalsOrig)),'.');
hold on;
plot(sort(interswarmReposeIntervals),1-((0:length(interswarmReposeIntervals)-1)'/length(interswarmReposeIntervals)),'.');
zoom on; grid on;

figure(); plot(tMain2,magMain(ccI),'.'); zoom on; grid on;
figure(); semilogy(sort(magMain(ccI)),1-((0:sum(ccI)-1)/sum(ccI)),'.'); zoom on; grid on;

% % figure();
% % semilogy(tdum,rateRatio,'.'); zoom on; grid on; hold on;
%
% % figure();
% % semilogy(tdum,rateRatio,'.'); zoom on; grid on; hold on;
% % semilogy(tdum,stackRate,'ko','linewidth',2); zoom on; grid on; hold on;
% %
% % figure(); semilogy(tdum,rateRatio,'.'); zoom on; grid on; hold on;
% %
% % figure(); semilogy(tdum,1./stackRate,'.'); zoom on; grid on; hold on;

