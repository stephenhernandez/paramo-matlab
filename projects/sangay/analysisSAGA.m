clear; close all; clc;

%%
cd ~/research/now/sangay/

%%
clear; [tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/research/now/sangay/sangaySubspaceDetectorSAGA_v2.mat',10); %,10);

t1 = tabs >= datetime(2018,11,11);
z2p(t1) = z2p(t1)/1.256780e+09;
z2p(~t1) = z2p(~t1)/3.141950e+08;
repeatThresh = 50; %50 seconds is definition of repeat

%%
z2p = z2p(:,1); %max(z2p,[],2);
kurt = kurt(:,1); %max(kurt,[],2);

%%
minAmp = 1e-6;
winlen = 50;
nboot = 1e3;
zeroPhaseFlag = 0;

%%
discr1 = z2p ./ kurt; %
%goodI = kurt >= 5 & kurt <= 30 & (NCC.*scaledCC./kurt) >= 5e-2; % & (NCC./scaledCC) >= 2e-4; % & (NCC./scaledCC) >= 4e-4;
goodI = discr1 >= 2e-7 & kurt <= 30 & kurt >= 5 & NCC >= 0.01 & z2p <= 2e-4;

%%
[N,edges] = histcounts(tabs(goodI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
N = N';
edges = edges(1:end-1)';
figure();
plot(edges,N,'.');
zoom on; grid on;

%%
figure();
semilogy(tabs(goodI),z2p(goodI),'.'); zoom on; grid on;

%%
figure();
plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;

%%
tGood = tabs(goodI);
aGood = z2p(goodI);
[acut,startIndex1,endIndex1] = cutWindows(aGood,winlen,winlen-1,false);
tcut = cutWindows(datenum(tGood),winlen,winlen-1,false);
tcut = dn2dt(tcut);
difft = seconds(diff(tcut));

%%
disp('getting variances, please be patient');
tic;
stepsize = max([1 ceil(winlen/20)]);
varIndices = (1:stepsize:size(acut,2))';
% vamp = stdcorr(log10(acut(:,1:10:size(acut,2))));
% vrate = stdcorr(difft(:,1:10:size(difft,2)));

% stdAmp = stdboot(log10(acut(:,varIndices)),nboot);
% stdRate = stdboot(difft(:,varIndices),nboot);

% stdAmp = stdcorr(log10(acut(:,varIndices)));
% stdRate = stdcorr(log10(difft(:,varIndices)));

p = [2.5 97.5];
ampBars = medboot(acut(:,varIndices),nboot,p);
rateBars = medboot(difft(:,varIndices),nboot,p);
toc;

%%
rollingBestimate = log10(exp(1))./(nanmedian(log10(acut)) - log10(minAmp) + 0.005);
rollingBestimate = rollingBestimate';
rollingBestimate = zpkFilter(rollingBestimate,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
medianAmplitude = nanmedian(acut)';
medianAmplitude = zpkFilter(medianAmplitude,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
rollingRate = 86400./nanmedian(difft);
rollingRate = rollingRate';
%rollingRate = medfiltSH(rollingRate,winlen,true);
rollingRate = zpkFilter(rollingRate,-inf,1/winlen,1,1,zeroPhaseFlag);

%%
figure();
ax = subplot(211);
plot(tGood(endIndex1),rollingBestimate,'.'); zoom on;  title('b-value');
ax(2) = subplot(212);
plot(tGood(endIndex1),rollingRate,'.'); zoom on;  title('rate'); grid on;
linkaxes(ax,'x');

%%
figure();
semilogy(tGood(endIndex1),rollingRate.*rollingBestimate,'.'); zoom on; title('rate * b-value'); grid on;

%%
figure();
ss = scatter(rollingBestimate,rollingRate,3*exp(log10(nanmedian(acut)')),datenum(tGood(endIndex1)),'filled');
zoom on; grid on; colorbar;
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log'; grid on;

%%
% ttable = tGood(endIndex1);
% figure('units','normalized','outerposition',[0 0 0.75 1]);
% r_m(1) = subplot(311);
% semilogy(ttable,rollingRate,'.'); zoom on; grid on; title('rate');
% hold on;
% semilogy(tGood(endIndex1(varIndices)),86400./rateBars,'-','Color',[0.5 0.5 0.5]);
% %
% r_m(2) = subplot(312);
% medianAmplitude = nanmedian(acut)'; % - log10(minAmp) + 0.005);
% semilogy(ttable,medianAmplitude,'.'); zoom on; grid on; title('median Amplitude');
% hold on;
% semilogy(tGood(endIndex1(varIndices)),ampBars,'-','color',[0.5 0.5 0.5]);
% %
% r_m(3) = subplot(313);
% energy_rate = rollingRate.*medianAmplitude;
% semilogy(tGood(endIndex1),energy_rate,'.'); zoom on; title('energy release rate'); grid on;
% linkaxes(r_m,'x');

%%
ttable = tGood(endIndex1);
figure('units','normalized','outerposition',[0 0 0.75 1]);
r_m(1) = subplot(311);
semilogy(ttable,rollingRate,'.'); zoom on; grid on; title('rate');
hold on;
semilogy(ttable(varIndices),86400./rateBars,'-','Color',[0.5 0.5 0.5]);
ax = gca;
ax.YScale = 'log';

%
r_m(2) = subplot(312); 
% - log10(minAmp) + 0.005);
semilogy(ttable,medianAmplitude,'.'); zoom on; grid on; title('median amplitude');
hold on;
semilogy(ttable(varIndices),ampBars,'-','color',[0.5 0.5 0.5]);
ax = gca;
ax.YScale = 'log';

%
r_m(3) = subplot(313);
semilogy(ttable,rollingRate.*medianAmplitude,'.'); zoom on; title('energy release rate'); grid on;
linkaxes(r_m,'x');

%%
figure('units','normalized','outerposition',[0 0 0.5 1]);
randI = randsample(size(acut,2),100);
for i = 1:length(randI)
    hold on; plot(sort(acut(:,randI(i))),1 - ((1:winlen)'/winlen),'.'); zoom on; grid on;
end
ax = gca;
ax.XScale = 'log';
title('amplitude distributions (CCDF)');
xlabel('[counts]');

%%
figure('units','normalized','outerposition',[0.5 0 0.5 1]);
randI = randsample(size(acut,2),100);
for i = 1:length(randI)
    hold on; plot(sort(difft(:,randI(i))),1 - ((1:winlen-1)'/(winlen-1)),'.'); zoom on; grid on;
end
ax = gca;
ax.XScale = 'log';
title('inter-event time distributions (CCDF)');
xlabel('[seconds]');

%%
figure('units','normalized','outerposition',[0 0 0.8 1]);
ax(1) = subplot(511);
semilogy(tabs(goodI),nanmedian(z2p(goodI,1),2),'.'); zoom on; grid on;
title('zero-to-peak');
ax(2) = subplot(512);
semilogy(tabs(goodI),(z2p(goodI,1)./kurt(goodI,1)),'.'); zoom on; grid on;
title('gamma');
ax(3) = subplot(513);
semilogy(tabs(goodI),nanmedian(kurt(goodI,1),2),'.'); zoom on; grid on;
title('kurtosis');
ax(4) = subplot(514);
plot(tabs(goodI),scaledCC(goodI),'.'); zoom on; grid on;
title('scaledCC');
ax(5) = subplot(515);
plot(tabs(goodI),NCC(goodI),'.'); zoom on; grid on;
title('NCC');
linkaxes(ax,'x');

%%
figure('units','normalized','outerposition',[0 0 0.8 1]);
Nsmooth = 100;
nDays = 1/2;
tsynth = datenum(ttable(1):seconds(86400/Nsmooth):ttable(end))';

ysynth = interp1(datenum(ttable),rollingRate,tsynth,'pchip');
ysynth = zpkFilter(ysynth,-inf,1/Nsmooth/nDays,1,1,zeroPhaseFlag);

ysynth2 = interp1(datenum(ttable),medianAmplitude,tsynth,'pchip');
ysynth2 = zpkFilter(ysynth2,-inf,1/Nsmooth/nDays,1,1,zeroPhaseFlag);

ysynth3 = ysynth.*ysynth2;

tsynth = dn2dt(tsynth);

load('~/products/rsam/EC.SAGA..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat','S');
gapInfo = S.gapInfo;
gI = gapInfo(:,2) >= 20;
tSAGA = getTimeVec(S);
gapStart = gapInfo(gI,1);
gapEnd = sum(gapInfo(gI,:),2)-1;
lGaps = length(gapStart);
lSAGA = S.npts;

tic;
for i = 1:lGaps
    disp(i);
    gapStart_ = gapStart(i);
    gapEnd_ = gapEnd(i);
    
    %%
    if gapStart_ < lSAGA && gapEnd_ <= lSAGA
        tsI = tsynth >= tSAGA(gapStart_) & tsynth <= tSAGA(gapEnd_);
        if sum(tsI)
            ysynth(tsI) = NaN;
            ysynth2(tsI) = NaN;
            ysynth3(tsI) = NaN;
        end
    end
end
toc;

axSmooth = subplot(311);
pp = plot(tsynth,ysynth,'linewidth',3);
pp.Color(4) = 0.75;
zoom on; grid on; ax = gca; ax.YScale = 'log';

axSmooth(2) = subplot(312);
pp = plot(tsynth,ysynth2,'linewidth',3);
pp.Color(4) = 0.75;
zoom on; grid on; ax = gca; ax.YScale = 'log';

axSmooth(3) = subplot(313);
pp = plot(tsynth,ysynth3,'linewidth',3);
pp.Color(4) = 0.75;
zoom on; grid on; ax = gca; ax.YScale = 'log';

linkaxes(axSmooth,'x');

%%
load sangay_excel_times.mat
excelI = texcel >= datetime(2013,01,01);
texcel = texcel(excelI);
aExcel = aExcel(excelI);
tCombined = [tGood; texcel];
iCombined = [ones(size(tGood)); zeros(size(texcel))];
[tCombined,sI] = sort(tCombined);
iCombined = iCombined(sI);
figure(); plot(tCombined,iCombined,'o'); zoom on;
figure(); semilogy(tCombined(1:end-1),seconds(diff(tCombined)),'.'); zoom on; grid on;
[tReduced,ii,removeIndices,keepIndices] = removeRepeatedMatches(tCombined,iCombined,repeatThresh);

%%
ampsCombined = [aGood; aExcel];
ampsCombined = ampsCombined(sI);
aReduced = ampsCombined;
figure(); tI = tReduced > datetime(2011,01,01); plot(tReduced(tI),(1:sum(tI))','.'); zoom on; grid on;
figure(); semilogy(tReduced(1:end-1),seconds(diff(tReduced)),'o'); zoom on; grid on; title('did i remove repeats?');

%%
[Nsvd,edgesSvd] = histcounts(tReduced(ii == 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
Nsvd = Nsvd';
edgesSvd = edgesSvd(1:end-1)';
figure(); plot(edgesSvd,Nsvd,'.'); zoom on;

%%
[Nexcel,edgesExcel] = histcounts(tReduced(ii ~= 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
Nexcel = Nexcel';
edgesExcel = edgesExcel(1:end-1)';
figure(17); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
disp([sum(Nsvd) sum(Nexcel)]);

%%
for i = 1:length(removeIndices)
rI = removeIndices{i};
aReduced(rI) = [];
end
figure(); 
semilogy(tReduced(ii == 1),aReduced(ii == 1),'.'); zoom on; grid on; hold on; semilogy(tReduced(ii ~= 1),aReduced(ii ~= 1),'p')

%%
[Nsvd,edgesSvd] = histcounts(tGood,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
Nsvd = Nsvd';
edgesSvd = edgesSvd(1:end-1)';
figure(); plot(edgesSvd,Nsvd,'.'); zoom on;

[Nexcel,edgesExcel] = histcounts(texcel,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
Nexcel = Nexcel';
edgesExcel = edgesExcel(1:end-1)';
figure(19); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
disp([sum(Nsvd) sum(Nexcel)]);

%%
figure(); histogram(Nexcel(Nexcel>0)./Nsvd(Nexcel>0),10000); zoom on; grid on;
figure(); plot(edgesExcel(Nexcel>0),Nexcel(Nexcel>0)./Nsvd(Nexcel>0),'.');
zoom on; grid on; 

title('Nexcel (Marcelo)/ Nsvd'); hold on; 
plot([edgesExcel(1) edgesExcel(end)],[1 1],'--','linewidth',2,'Color',[0.5 0.5 0.5]);
tCommon = [];
aRatioCommon = tCommon;

for i = 1:length(keepIndices)
kI = keepIndices{i};
rI = removeIndices{i};
tCommon = [tCommon; tCombined(kI)];
aRatioCommon = [aRatioCommon; ampsCombined(kI)./ampsCombined(rI)];
end

[tCommon,sI2] = sort(tCommon);
aRatioCommon = aRatioCommon(sI2);

%%
figure(); semilogy(tCommon,aRatioCommon,'.'); zoom on; grid on;

%% do not erase
% z2p_seismic = z2p; load('REVS_infrasoundZero2Peak_forCASC','z2p','rI');
% z2p_seismic = z2p_seismic(goodI);
% %%
% figure(); semilogy(tGood(~rI),z2p./z2p_seismic(~rI),'.'); zoom on; grid on;
% tVASR = tGood(~rI);
% [VASRcut,~,endIndex1] = cutWindows(z2p./z2p_seismic(~rI),winlen,winlen-1,false);
% medVASR = nanmedian(VASRcut);
%
% %%
% hold on;
% pp = semilogy(tVASR(endIndex1),medVASR,'linewidth',5); zoom on; grid on;
% pp.Color(4) = 0.75;

