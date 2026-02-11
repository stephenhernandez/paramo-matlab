clear; close all; %clc;

cd ~/research/now/cotopaxi/
%cd ~/subspace_detector/
%cd ~/igdata/

%[tabs,z2p,NCC,Neff,templateIndex,scaledCC,kurt] = ...
%    filterUniqueEvents('~/research/now/cotopaxi/cotopaxiSubspaceDetectorBREF_v2',30,15);

[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('cotopaxiSubspaceDetectorBREF_v8',30); %,10,10);

%%
option2 = true;
if option2
    minKurt = min(kurt,[],2,'omitnan');
    maxKurt = max(kurt,[],2,'omitnan');
    kurtRatio = maxKurt./minKurt;
    z2pOrig = z2p;
    kurtOrig = kurt;
    p2rmsOrig = p2rms;
    %z2p = median(z2pOrig(:,1:3),2,'omitnan');
    %z2p = median(z2pOrig,2,'omitnan');
    %z2p(:,7:9) = []; z2p = max(z2p(:,1:9),[],2,'omitnan');
    z2p(:,7:9) = []; z2p = median(z2p(:,1:9),2,'omitnan');
    kurt = median(kurtOrig,2,'omitnan');
    p2rms = median(p2rmsOrig,2,'omitnan');

    minAmp = 2e2;
    winlen = 21;
    nboot = 5e2;

    goodI = median(z2p,2,'omitnan') >= minAmp & ...
        NCC >= 0.15 & ...
        mad(p2rmsOrig,1,2)./median(p2rmsOrig,2,'omitnan') >= 0.05 & ...
        mad(kurtOrig,1,2)./median(kurtOrig,2,'omitnan') >= 0.1 & ...
        median(kurtOrig,2,'omitnan') >= 3.5 & ...
        median(p2rmsOrig,2,'omitnan') >= 3.5 & ...
        max(kurtOrig,[],2,'omitnan')./median(kurtOrig,2,'omitnan') >= 1.4 & ...
        max(p2rmsOrig,[],2) >= 4.5 & ...
        max(kurtOrig,[],2) >= 4 &...
        tabs >= datetime(2000,08,01) & ... %tabs < datetime(2023,10,01) & ...
        Neff > 9;

    % goodI = median(z2p,2,'omitnan') >= minAmp & ...
    %     NCC >= 0.15 & ...
    %     max(p2rmsOrig,[],2) >= 4 & ... %4.5 & ...
    %     max(p2rmsOrig,[],2) <= 9 & ...
    %     tabs >= datetime(2000,12,01) & ...
    %     max(kurtOrig,[],2) >= 4 &...
    %     max(kurtOrig,[],2) <= 20;

    tGood = tabs(goodI);
    aGood = z2p(goodI);
else
    minKurt = min(kurt(:,1:2),[],2,'omitnan');
    maxKurt = max(kurt(:,1:2),[],2,'omitnan');
    kurtRatio = maxKurt./minKurt;

    z2pOrig = z2p;
    kurtOrig = kurt;
    p2rmsOrig = p2rms;
    z2p = median(z2pOrig,2,'omitnan');
    kurt = median(kurtOrig,2,'omitnan');
    p2rms = median(p2rmsOrig,2,'omitnan');

    minAmp = 3e2;
    winlen = 51;
    nboot = 5e2;

    % goodI = median(z2pOrig,2,'omitnan') >= minAmp & ...
    %     NCC >= 0.15 & ... %0.20 & ...
    %     max(p2rmsOrig,[],2) >= 4.5 & ...
    %     max(p2rmsOrig,[],2) <= 9 & ...
    %     tabs >= datetime(2000,12,01) & ...
    %     max(kurtOrig,[],2) >= 4 &...
    %     max(kurtOrig,[],2) <= 20;

    goodI = z2p >= minAmp & ...
        NCC >= 0.15 & ...
        tabs >= datetime(2000,12,01) & ...
        mad(kurtOrig,1,2)./median(kurtOrig,2,'omitnan') >= 0.2;

    tGood = tabs(goodI);
    aGood = z2p(goodI);
end

%%
% ccPredicted = 10.^(b(1) + b(2).*log10(scaledCC));
% error = abs(NCC - ccPredicted);
% eI2 = error <= 0.15;

% minAmp = 3e2;
% winlen = 51;
% nboot = 5e2;
%
% %goodI = z2p >= minAmp minAmp & scaledCC >= 8 & kurt >= 10 & kurt < 50 & NCC >= 0.01 & tabs >= datetime(2011,01,01); % & discr2 < 1 & discr2 > 0.01; %& p2rms >= 5 &p2rms <= 12.5; % & kurt <= 70 & NCC >= 0.35; % & NCC >= 0.4 & kurt >= 5 & kurt <= 25; % discr2 >= 0.03 & discr2 <= 1; winlen = 100;
% %goodI = z2p >= minAmp & NCC >= 0.3; % & scaledCC >= 15 & kurtRatio <= 5 & kurt >= 10 & tabs >= datetime(2011,01,01); % & discr2 < 1 & discr2 > 0.01; %& p2rms >= 5 &p2rms <= 12.5; % & kurt <= 70 & NCC >= 0.35; % & NCC >= 0.4 & kurt >= 5 & kurt <= 25; % discr2 >= 0.03 & discr2 <= 1; winlen = 100;
% %goodI = z2p >= minAmp & NCC >= 0.2 & isfinite(z2p(:,1)) & Neff > 0 & z2p <= 1e5; % & tabs >= datetime(2020,12,01);
% goodI = nanmedian(z2pOrig(:,1:3),2) >= minAmp & ...
%     NCC >= 0.20 & ...
%     max(p2rmsOrig,[],2) >= 4.5 & ...
%     max(p2rmsOrig,[],2) <= 9 & ...
%     tabs >= datetime(2000,12,01) & ...
%     max(kurtOrig,[],2) >= 4 &...
%     max(kurtOrig,[],2) <= 20;
% z2p = nanmedian(z2pOrig(:,1:3),2);
% tGood = tabs(goodI);
% aGood = z2p(goodI);

%%
[N,edges] = histcounts(tabs(goodI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
N = N';
edges = edges(1:end-1)';

figure();
plot(edges,N,'.');
zoom on; grid on;

%
figure();
semilogy(tabs(goodI),z2p(goodI),'.'); zoom on; grid on;

%
figure();
plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;

%tGood = tsheet(refI);
%aGood = amp(refI);minAmp = 4e2; winlen = 25; nboot = 2e2;
[acut,startIndex1,endIndex1] = cutWindows(aGood,winlen,winlen-1,false);
tcut = cutWindows(datenum(tGood),winlen,winlen-1,false);
tcut = dn2dt(tcut);
difft = seconds(diff(tcut));
disp('getting variances, please be patient');
tic;
stepsize = max([1 ceil(winlen/20)]);
varIndices = (1:stepsize:size(acut,2))';
p = [2.5 97.5];
ampBars = medboot(acut(:,varIndices),nboot,p);
rateBars = medboot(difft(:,varIndices),nboot,p);
toc;

rollingBestimate = log10(exp(1))./(nanmedian(log10(acut)) - log10(minAmp) + 0.005);
rollingBestimate = rollingBestimate';
rollingRate = 86400./nanmedian(difft);
rollingRate = rollingRate';

%
close all;

figure(); semilogy(tGood(endIndex1),rollingRate,'.'); zoom on;  title('rate'); grid on;
%figure(); histogram(rollingRate(tGood(endIndex1)>=datetime(2019,01,01)),logspace(floor(min(log10(rollingRate(tGood(endIndex1)>=datetime(2019,01,01))))),ceil(max(log10(rollingRate(tGood(endIndex1)>=datetime(2019,01,01))))),151)); zoom on; grid on; ax = gca; ax.XScale = 'log';

tEnd = tGood(endIndex1);
[tEnd,sI] = sort(tEnd);
rollingRate = rollingRate(sI);

tEnd = tEnd(end:-winlen:1);
rateGood = (rollingRate(end:-winlen:1));

%figure(); histogram(rateGood(tEnd>=datetime(2019,01,01)),logspace(floor(min(log10(rateGood(tEnd>=datetime(2019,01,01))))),ceil(max(log10(rateGood(tEnd>=datetime(2019,01,01))))),151)); zoom on; grid on; ax = gca; ax.XScale = 'log';
figure(); stem(tEnd,(rateGood),'o'); zoom on; grid on;

%%
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
rollingRate = 86400./nanmedian(difft);
rollingRate = rollingRate';

%%
figure();
plot(tGood(endIndex1),rollingBestimate,'.'); zoom on;  title('b-value');

figure();
plot(tGood(endIndex1),rollingRate,'.'); zoom on;  title('rate'); grid on;

%%
figure();
semilogy(tGood(endIndex1),rollingRate.*rollingBestimate,'.'); zoom on; title('rate * b-value'); grid on;

%%
figure();
ss = scatter(rollingBestimate,rollingRate,3*exp(log10(nanmean(acut)')),datenum(tGood(endIndex1)),'filled');
zoom on; grid on; colorbar;
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log'; grid on;

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
medianAmplitude = nanmedian(acut)'; % - log10(minAmp) + 0.005);
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
ax.YScale = 'log';
title('amplitude distributions (CCDF)');
xlabel('[counts]');

%
figure('units','normalized','outerposition',[0.5 0 0.5 1]);
randI = randsample(size(acut,2),100);
for i = 1:length(randI)
    hold on; plot(sort(difft(:,randI(i))),1 - ((1:winlen-1)'/(winlen-1)),'.'); zoom on; grid on;
end
ax = gca;
ax.XScale = 'log';
title('inter-event time distributions (CCDF)');
xlabel('[seconds]');

%close all;
nDays = 3;
[rate,~,medianMagsFixedTimeWin] = ...
    t2r(tGood,days(nDays),log10(aGood)-1);

axL1 = linkedPlot(tGood,...
    rate/nDays,...
    medianMagsFixedTimeWin,...
    (10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays);
grid(axL1,'on')
axL1(1).YScale = 'log';
axL1(3).YScale = 'log'; 

axL2 = linkedPlot(tGood,...
    aGood,...
    NCC(goodI)); 
axL2(1).YScale = 'log';

[N,edges] = histcounts(tGood,(dateshift(min(tGood),'start','day'):dateshift(max(tGood),'end','day'))');
N = N';

edges = edges(1:end-1);
figure(); stairs(edges,N,'linewidth',2); zoom on;
grid on;
ylabel('Numero Diario');
title('Sismicidad Cotopaxi');

plotAllFlag = false;
if plotAllFlag
    figure();
    ax5(1) = subplot(511);
    ax5(1) = subplot(512);
    ax5(1) = subplot(513);
    ax5(1) = subplot(514);
    ax5(1) = subplot(515);
    linkaxes(ax5,'x')
end

% %% do not erase
%
% % pI = true(size(kurt)); %(kurt >= 9.5 & kurt < 40) | (p2rms >=7 & p2rms < 12); %pksOrig >= 0.2 | (p2rms >= 7 & p2rms < 12);
% % figure('units','normalized','outerposition',[0 0 0.8 1]);
% % ax(1) = subplot(511);
% % semilogy(tabs(goodI),nanmedian(z2p(goodI,1),2),'.'); zoom on; grid on;
% % title('zero-to-peak');
% % ax(2) = subplot(512);
% % semilogy(tabs(goodI),(z2p(goodI,1)./kurt(goodI,1)),'.'); zoom on; grid on;
% % title('gamma');
% % ax(3) = subplot(513);
% % semilogy(tabs(goodI),nanmedian(kurt(goodI,1),2),'.'); zoom on; grid on;
% % title('kurtosis');
% % ax(4) = subplot(514);
% % plot(tabs(goodI),scaledCC(goodI),'.'); zoom on; grid on;
% % title('scaledCC');
% % ax(5) = subplot(515);
% % plot(tabs(goodI),NCC(goodI),'.'); zoom on; grid on;
% % title('NCC');
% % linkaxes(ax,'x')
% % %%
% % tsynth = (datenum(dateshift(ttable(1),'end','day')):1/48:datenum(dateshift(ttable(end),'start','day')))';
% % ysynth = interp1(datenum(ttable),(1:length(ttable))',tsynth);
% % figure();
% % pp = plot(dn2dt(tsynth(1:end-1)),zpkFilter(diff(ysynth)./diff(tsynth),-inf,1/480,1,1,1),'linewidth',4);
% % pp.Color(4) = 0.5; zoom on; grid on; ax = gca; ax.YScale = 'log';
% % %%
% % load cotopaxi_excel_times.mat
% % excelI = texcel >= datetime(2011,01,01);
% % texcel = texcel(excelI);
% % aExcel = aExcel(excelI);
% % tCombined = [tGood; texcel];
% % iCombined = [ones(size(tGood)); zeros(size(texcel))];
% % [tCombined,sI] = sort(tCombined);
% % iCombined = iCombined(sI);
% % figure(); plot(tCombined,iCombined,'o'); zoom on;
% % figure(); semilogy(tCombined(1:end-1),seconds(diff(tCombined)),'.'); zoom on; grid on;
% % [tReduced,ii,removeIndices,keepIndices] = removeRepeatedMatches(tCombined,iCombined,30);
% % %%
% % ampsCombined = [aGood; aExcel];
% % ampsCombined = ampsCombined(sI);
% % aReduced = ampsCombined;
% % figure(); tI = tReduced > datetime(2011,01,01); plot(tReduced(tI),(1:sum(tI))','.'); zoom on; grid on;
% % figure(); semilogy(tReduced(1:end-1),seconds(diff(tReduced)),'o'); zoom on; grid on; title('did i remove repeats?');
% % %%
% % [Nsvd,edgesSvd] = histcounts(tReduced(ii == 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
% % Nsvd = Nsvd';
% % edgesSvd = edgesSvd(1:end-1)';
% % figure(); plot(edgesSvd,Nsvd,'.'); zoom on;
% % %%
% % [Nexcel,edgesExcel] = histcounts(tReduced(ii ~= 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
% % Nexcel = Nexcel';
% % edgesExcel = edgesExcel(1:end-1)';
% % figure(17); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
% % disp([sum(Nsvd) sum(Nexcel)]);
% % %%
% % for i = 1:length(removeIndices)
% % rI = removeIndices{i};
% % aReduced(rI) = [];
% % end
% % figure(); semilogy(tReduced(ii == 1),aReduced(ii == 1),'.'); zoom on; grid on; hold on; semilogy(tReduced(ii ~= 1),aReduced(ii ~= 1),'p')
% % S = loadWaveforms(datetime(2020,05,29),1,"BREF","BHZ");
% % %%
% % [Nsvd,edgesSvd] = histcounts(tGood,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
% % Nsvd = Nsvd';
% % edgesSvd = edgesSvd(1:end-1)';
% % figure(); plot(edgesSvd,Nsvd,'.'); zoom on;
% % [Nexcel,edgesExcel] = histcounts(texcel,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
% % Nexcel = Nexcel';
% % edgesExcel = edgesExcel(1:end-1)';
% % figure(19); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
% % disp([sum(Nsvd) sum(Nexcel)]);
% % %%
% % figure(); histogram(Nexcel(Nexcel>0)./Nsvd(Nexcel>0),10000); zoom on; grid on;
% % figure(); plot(edgesExcel(Nexcel>0),Nexcel(Nexcel>0)./Nsvd(Nexcel>0),'.');
% % zoom on; grid on; title('Nexcel (vero)/ Nsvd'); hold on; plot([edgesExcel(1) edgesExcel(end)],[1 1],'--','linewidth',2,'Color',[0.5 0.5 0.5]);
% % tCommon = [];
% % aRatioCommon = tCommon;
% % for i = 1:length(keepIndices)
% % kI = keepIndices{i};
% % rI = removeIndices{i};
% % tCommon = [tCommon; tCombined(kI)];
% % aRatioCommon = [aRatioCommon; ampsCombined(kI)./ampsCombined(rI)];
% % end
% % [tCommon,sI2] = sort(tCommon);
% % aRatioCommon = aRatioCommon(sI2);
% % %%
% % figure(); semilogy(tCommon,aRatioCommon,'.'); zoom on; grid on;
% % sum(ii)
% % sum(ii) - length(tCommon)
% % length(tCommon)
% % sum(~ii)
% % length(ii)
%
% %%
% % z2p_seismic = z2p; load('REVS_infrasoundZero2Peak_forCASC','z2p','rI');
% % z2p_seismic = z2p_seismic(goodI);
% % %%
% % figure(); semilogy(tGood(~rI),z2p./z2p_seismic(~rI),'.'); zoom on; grid on;
% % tVASR = tGood(~rI);
% % [VASRcut,~,endIndex1] = cutWindows(z2p./z2p_seismic(~rI),winlen,winlen-1,false);
% % medVASR = nanmedian(VASRcut);
% %
% % %%
% % hold on;
% % pp = semilogy(tVASR(endIndex1),medVASR,'linewidth',5); zoom on; grid on;
% % pp.Color(4) = 0.75;
%
