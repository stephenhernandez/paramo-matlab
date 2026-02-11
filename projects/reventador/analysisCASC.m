clear; close all; clc;

%
%cd ~/research/now/reventador/
% load('CASC_SubspaceDetector_Results_3');
% load('CASC_SubspaceDetector_Results_2','b');

%
%[tabs,z2p,NCC,Neff,templateIndex,scaledCC,kurt] = ...
%    filterUniqueEvents('~/research/now/reventador/CASC_SubspaceDetector_Results_v4',30,20);
%tPotential,z2p,NCC,Neff,p2rms,kurt

%[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
%    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v6',10); %,10);
% [tabs,z2p,NCC,Neff,p2rms,kurt] = ...
%     filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v7',10); %,10);
%[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
%    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v10',10); %,10);
[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/masa/old/research/now/reventador/ReventadorSubspaceDetectorResults_v10',10); %,10);
tabs = tabs + seconds(10);
rateSpectrumFlag = true;

%%
minKurt = min(kurt,[],2,'omitnan');
maxKurt = max(kurt,[],2,'omitnan');
kurtRatio = maxKurt./minKurt;

%
refEllipse = referenceEllipsoid('wgs84');
[stla,stlo] = metaDataFromStationList(["CASC";"BONI";"ANTS";"ANTG"]);
d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse);
sensitivities = [3.141950e+08; 2.01494e9; 5.03735e8];

%
z2pOrig = z2p;
z2p(:,1:3) = d_(1)*d_(1)*z2p(:,1:3)/sensitivities(1);

boniI = tabs <= datetime(2020,01,339);
z2p(boniI,4:6) = d_(2)*d_(2)*z2p(boniI,4:6)/sensitivities(2);
z2p(~boniI,4:6) = d_(2)*d_(2)*z2p(~boniI,4:6)/sensitivities(3);

antsI = tabs <= datetime(2014,08,14);
z2p(antsI,7:9) = d_(3)*d_(3)*z2p(antsI,7:9)/sensitivities(2);
z2p(~antsI,7:9) = d_(3)*d_(3)*z2p(~antsI,7:9)/sensitivities(3);

z2p(:,10:12) = d_(4)*d_(4)*z2p(:,10:12)/sensitivities(1);

%
z2p = median([median(z2p(:,1:3),2,'omitnan') ...
    median(z2p(:,4:6),2,'omitnan') ...
    median(z2p(:,7:9),2,'omitnan') ...
    median(z2p(:,10:12),2,'omitnan')],2,'omitnan');
%z2p = median(z2p,2,'omitnan');
%kurt = max(kurt,[],2);

%
%ccPredicted = b(1) + b(2).*log10(discr);
%error = abs(NCC - ccPredicted);
%eI2 = error <= 0.35;

%
minAmp = 300; %2e3;  %1e-2;
maxAmp = 1e5; %1e4;
winlen = 101; %75;
nboot = 6e2;
zeroPhaseFlag = false;

%
%discr1 = z2p.*kurt;
%discr2 = z2p./kurt;
%goodI = (eI2 & z2p >= minAmp) & kurt >= 10 & scaledCC >= 4 & discr2 >= 1 & kurt <= 25; %p2rms >= 4 & p2rms <= 10 & NCC >= 0.4 & kurt >= 5 & kurt <= 25; % discr2 >= 0.03 & discr2 <= 1; winlen = 100;
%goodI = z2p >= minAmp & kurt <= 10 & discr2 <= 100 & NCC >= 0.01; %; % & kurt >= 5  & kurtRatio <= 3 & kurt >= 1
%goodI = z2p >= minAmp & NCC >= 0.02;
%goodI = z2p >= minAmp & kurt <= 20;
%goodI = z2p >= minAmp & NCC >= 0.04 & kurt <= 25 & kurt >= 8 & tabs >= datetime(2013,01,01);
%goodI = (z2p >= minAmp & kurt <= 25 & kurt >= 5 & tabs >= datetime(2013,01,01)) & ~(Neff > 1 & NCC < 0.03);
%goodI = (z2p >= minAmp & kurt <= 25 & kurt >= 5 & tabs >= datetime(2013,01,01)) & ~(Neff < 2 & scaledCC < 50) & ~(Neff > 1 & NCC < 0.03);
%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 1 & scaledCC < 50) & kurt < 25;
%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 1) & kurt < 25;
%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 0) & kurt < 25;

% goodI = z2p >= minAmp & z2p <= maxAmp & NCC >= 0.6 & Neff >= 12 & tabs >= datetime(2018,01,01);
goodI = tabs >= datetime(2016,01,01) & z2p >= minAmp & z2p <= maxAmp & NCC >= 0.35 & ...
    median(p2rms,2,'omitnan') >= 3 & max(p2rms,[],2,'omitnan') <= 7 & ...
    min(p2rms,[],2,'omitnan') >= 2.5 & median(kurt,2,'omitnan') >= 3 & ...
    max(kurt,[],2,'omitnan') <= 20 & min(kurt,[],2,'omitnan') >= 2.5;

%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.4) & ~(Neff == 0) & kurt < 20 & tabs >= datetime(2016,07,01) & tabs <= datetime(2020,01,01);
%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.4) & ~(Neff == 0) & kurt < 20 & tabs <= datetime(2014,08,01);

% figure 1
[N,edges] = histcounts(tabs(goodI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
N = N';
edges = edges(1:end-1)';
figure();
stairs(edges,N,'.-');
zoom on; grid on;

%
tGood = tabs(goodI);
aGood = z2p(goodI);

figure('units','normalized','outerposition',[0 0 0.75 1]);
ax2(1) = subplot(211);
semilogy(tGood,aGood,'.'); zoom on; grid on;
ax2(2) = subplot(212);
plot(tGood,(1:length(tGood))','.'); zoom on; grid on;
linkaxes(ax2,'x');

%%
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
close all; figure(); semilogy(tGood(endIndex1),rollingRate,'.'); zoom on;  title('rate'); grid on;
figure(); histogram(rollingRate(tGood(endIndex1)>=datetime(2019,01,01)),logspace(floor(min(log10(rollingRate(tGood(endIndex1)>=datetime(2019,01,01))))),ceil(max(log10(rollingRate(tGood(endIndex1)>=datetime(2019,01,01))))),151)); zoom on; grid on; ax = gca; ax.XScale = 'log';
tEnd = tGood(endIndex1);
[tEnd,sI] = sort(tEnd(end:-winlen:1));
rateGood = (rollingRate(end:-winlen:1));
rateGood = rateGood(sI);
figure(); histogram(rateGood(tEnd>=datetime(2019,01,01)),logspace(floor(min(log10(rateGood(tEnd>=datetime(2019,01,01))))),ceil(max(log10(rateGood(tEnd>=datetime(2019,01,01))))),151)); zoom on; grid on; ax = gca; ax.XScale = 'log';
figure(); stem(tEnd,(rateGood),'.'); zoom on; grid on;

%% figure 2
% figure('units','normalized','outerposition',[0 0 0.75 1]);
% ax2(1) = subplot(211);
% semilogy(tGood,aGood,'.'); zoom on; grid on;
% ax2(2) = subplot(212);
% plot(tGood,(1:length(tGood))','.'); zoom on; grid on;
% linkaxes(ax2,'x');
%hold on;
%semilogy(tabs(~goodI),z2p(~goodI),'.'); zoom on; grid on;

%%
% figure();
% plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;

%%
% tStart = [datetime(2014,04,01); ...
%     datetime(2016,11,10);...
%     datetime(2017,09,23);...
%     datetime(2018,04,02);...
%     datetime(2018,07,11);...
%     datetime(2018,09,19);...
%     datetime(2018,11,24)];
%
% tEnd = [datetime(2014,07,23);
%     datetime(2017,09,22);...
%     datetime(2018,04,01);...
%     datetime(2018,07,10);...
%     datetime(2018,09,18);...
%     datetime(2018,11,23);...
%     datetime(2019,03,27)];

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
semilogy(tGood(endIndex1),rollingRate./rollingBestimate,'.'); zoom on; title('rate / b-value'); grid on;

%%
% figure();
% ss = scatter(rollingBestimate,rollingRate,3*exp(log10(nanmedian(acut)')),datenum(tGood(endIndex1)),'filled');
% zoom on; grid on; colorbar;
% ax = gca;
% ax.YScale = 'log';
% ax.XScale = 'log'; grid on;

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
% figure('units','normalized','outerposition',[0 0 0.8 1]);
% ax(1) = subplot(611);
% semilogy(tabs(goodI),nanmedian(z2p(goodI,1),2),'.'); zoom on; grid on;
% title('zero-to-peak');
% ax(2) = subplot(612);
% semilogy(tabs(goodI),(z2p(goodI,1)./kurt(goodI,1)),'.'); zoom on; grid on;
% title('gamma');
% ax(3) = subplot(613);
% semilogy(tabs(goodI),nanmedian(kurt(goodI,1),2),'.'); zoom on; grid on;
% title('kurtosis');
% ax(4) = subplot(614);
% plot(tabs(goodI),scaledCC(goodI),'.'); zoom on; grid on;
% title('scaledCC');
% ax(5) = subplot(615);
% plot(tabs(goodI),NCC(goodI),'.'); zoom on; grid on;
% title('NCC');
% ax(6) = subplot(616);
% plot(tabs(goodI),kurtRatio(goodI),'.'); zoom on; grid on;
% title('kurtRatio');
% linkaxes(ax,'x');

%%
% tsynth = (datenum(dateshift(ttable(1),'end','day')):1/48:datenum(dateshift(ttable(end),'start','day')))';
% ysynth = interp1(datenum(ttable),(1:length(ttable))',tsynth);
% figure();
% pp = plot(dn2dt(tsynth(1:end-1)),zpkFilter(diff(ysynth)./diff(tsynth),-inf,1/480,1,1,1),'linewidth',4);
% pp.Color(4) = 0.5; zoom on; grid on; ax = gca; ax.YScale = 'log';

%tsynth = datenum(dateshift(ttable(1),'end','day'):seconds(1728):dateshift(ttable(end),'start','day'))';
Nsmooth = 100;
%rollingRate = medfiltSH(rollingRate,10,true);
tsynth = datenum(ttable(1):seconds(86400/Nsmooth):ttable(end))';
ysynth = interp1(datenum(ttable),rollingRate,tsynth,'pchip');
%ysynth = medfiltSH(ysynth,10,true);


%%
if rateSpectrumFlag
    t = tGood;
    amps = aGood;

    clearvars -except t amps rateSpectrumFlag;
    tI = t >= datetime(2010,01,01); % & amps >= 1e3;
    t = t(tI); 
    amps = amps(tI);

    hourlyRate = t2r(t,hours(1));
    nDays = 1; dailyRate = t2r(t,days(nDays))/nDays;
    tHour = (dateshift(min(t),'start','hour'):hours(1):dateshift(max(t),'end','hour'))';
    tDay = (dateshift(min(t),'start','day'):days(nDays):dateshift(max(t),'end','day'))';
    rHour2 = interp1(datenum(t),hourlyRate,datenum(tHour));
    rDay2 = interp1(datenum(t),dailyRate,datenum(tDay));
    rDay2(~isfinite(rDay2)) = 0;
    rHour2(~isfinite(rHour2)) = 0;

    close all; 
    figure(); plot(t,[hourlyRate dailyRate],'.'); zoom on;
    figure(); plot(tHour,rHour2,'.'); zoom on; grid on; hold on; plot(tDay,rDay2,'.');
    figure(); ax(1) = subplot(211); semilogy(t,amps,'.'); zoom on; grid on; 
    ax(2) = subplot(212); plot(t,dailyRate,'.'); 
    linkaxes(ax,'x');

    [pxx,fxx] = pwelch(detrend(rDay2),blackmanharris(256),[],[],1/nDays);
    [pxx2,fxx2] = pwelch(detrend(24*rHour2),blackmanharris(4096),[],[],24);

    figure(); p = loglog(1./fxx,pxx,'-'); zoom on; grid on;
    figure(4); hold on; p = loglog(1./fxx2,pxx2,'-'); zoom on; grid on;
    xlabel('days per cycle');

    [pxx,fxx] = pmtm(detrend(rDay2),4,[],1/nDays);
    [pxx2,fxx2] = pmtm(detrend(24*rHour2),4,[],24);

    figure(); p = loglog(1./fxx,pxx,'-'); zoom on; grid on;
    hold on; p = loglog(1./fxx2,pxx2,'-'); zoom on; grid on;
    xlabel('days per cycle');

    txx = 1./fxx(2:end);
    newTxx = (log10(min(txx)):0.01:log10(max(txx)))';
    pxxInterp = interp1(log10(txx),log10(pxx(2:end)),newTxx);

    figure(); 
    loglog(10.^newTxx,10.^pxxInterp,'-','linewidth',2); zoom on; grid on;
    b_ = robustfit(newTxx,pxxInterp);
    b_ = flipud(b_);
    disp(b_);
    yfit = polyval(b_,newTxx);
    hold on; loglog(10.^newTxx,10.^yfit,'-','linewidth',2); zoom on; grid on;
    %
    ylabel('$PSD$')
    xlabel('$Days/Cycle$')
    legend('PSD of rate data',sprintf('best fit slope=%f',b_(1)),'location','northwest');
end

%%
close all;
nDays = 1;
axL1 = linkedPlot(tGood,aGood,t2r(tGood,days(nDays))/nDays);
axL1(1).YScale = 'log';
axL1(1).YLabel.String = "Amplitud Explosiones [cuentas]";
axL1(2).YLabel.String = "Eventos / Día";
axis tight;

nDays = 14; [rate,~,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = ...
t2r(tGood,days(nDays),log10(aGood));
%axL1 = linkedPlot(tGood,aGood,medianMagsFixedTimeWin,1.5*medianMagsFixedTimeWin.*rate/nDays,rate/nDays); axL1(1).YScale = 'log';
axL1 = linkedPlot(tGood,aGood,medianMagsFixedTimeWin,rate/nDays,(10.^(1.5*medianMagsFixedTimeWin)).*rate/nDays); axL1(1).YScale = 'log'; axL1(end).YScale = 'log';

%% run this, do not erase, temporarily commentd out on 12 apr 2022
% figure('units','normalized','outerposition',[0 0 0.8 1]);
% figure(11); hold on;
% ysynth = zpkFilter(ysynth,-inf,1/Nsmooth/2,1,1,zeroPhaseFlag);
% tsynth = dn2dt(tsynth);
%
% load('~/products/rsam/EC.CASC..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat','S');
% gapInfo = S.gapInfo;
% gI = gapInfo(:,2) >= 20;
% tCASC = getTimeVec(S);
% gapStart = gapInfo(gI,1);
% gapEnd = sum(gapInfo(gI,:),2)-1;
% lGaps = length(gapStart);
% lCASC = S.npts;
% 
% %%
% tic;
% for i = 1:lGaps
%     disp(i);
%     gapStart_ = gapStart(i);
%     gapEnd_ = gapEnd(i);
%     
%     %%
%     if gapStart_ < lCASC && gapEnd_ <= lCASC
%         tsI = tsynth >= tCASC(gapStart_) & tsynth <= tCASC(gapEnd_);
%         if sum(tsI)
%             ysynth(tsI) = NaN;
%         end
%     end
% end
% toc;
% 
% %%
% pp = plot(tsynth,ysynth,'linewidth',3);
% pp.Color(4) = 0.75;
% zoom on; grid on; ax = gca; ax.YScale = 'log';
% 
% %%
% load reventador_excel_times.mat
% excelI = texcel >= datetime(2011,01,01);
% texcel = texcel(excelI);
% aExcel = aExcel(excelI);
% tCombined = [tGood; texcel];
% iCombined = [ones(size(tGood)); zeros(size(texcel))];
% [tCombined,sI] = sort(tCombined);
% iCombined = iCombined(sI);
% figure(); plot(tCombined,iCombined,'o'); zoom on;
% figure(); semilogy(tCombined(1:end-1),seconds(diff(tCombined)),'.'); zoom on; grid on;
% [tReduced,ii,removeIndices,keepIndices] = removeRepeatedMatches(tCombined,iCombined,30);
% 
% %%
% ampsCombined = [aGood; aExcel];
% ampsCombined = ampsCombined(sI);
% aReduced = ampsCombined;
% figure(); tI = tReduced > datetime(2011,01,01); plot(tReduced(tI),(1:sum(tI))','.'); zoom on; grid on;
% figure(); semilogy(tReduced(1:end-1),seconds(diff(tReduced)),'o'); zoom on; grid on; title('did i remove repeats?');
% 
% %%
% [Nsvd,edgesSvd] = histcounts(tReduced(ii == 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
% Nsvd = Nsvd';
% edgesSvd = edgesSvd(1:end-1)';
% figure(); plot(edgesSvd,Nsvd,'.'); zoom on;
% 
% %%
% [Nexcel,edgesExcel] = histcounts(tReduced(ii ~= 1),dateshift(min(tReduced),'start','day'):dateshift(max(tReduced),'end','day'));
% Nexcel = Nexcel';
% edgesExcel = edgesExcel(1:end-1)';
% figure(17); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
% disp([sum(Nsvd) sum(Nexcel)]);
% 
% %%
% for i = 1:length(removeIndices)
% rI = removeIndices{i};
% aReduced(rI) = [];
% end
% figure(); 
% semilogy(tReduced(ii == 1),aReduced(ii == 1),'.'); zoom on; grid on; hold on; semilogy(tReduced(ii ~= 1),aReduced(ii ~= 1),'p')
% 
% %%
% [Nsvd,edgesSvd] = histcounts(tGood,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
% Nsvd = Nsvd';
% edgesSvd = edgesSvd(1:end-1)';
% figure(); plot(edgesSvd,Nsvd,'.'); zoom on;
% 
% [Nexcel,edgesExcel] = histcounts(texcel,dateshift(min(texcel),'start','day'):dateshift(max(texcel),'end','day'));
% Nexcel = Nexcel';
% edgesExcel = edgesExcel(1:end-1)';
% figure(19); hold on; plot(edgesExcel,Nexcel,'s'); zoom on;
% disp([sum(Nsvd) sum(Nexcel)]);
% 
% %%
% figure(); histogram(Nexcel(Nexcel>0)./Nsvd(Nexcel>0),10000); zoom on; grid on;
% figure(); plot(edgesExcel(Nexcel>0),Nexcel(Nexcel>0)./Nsvd(Nexcel>0),'.');
% zoom on; grid on; title('Nexcel (vero)/ Nsvd'); hold on; plot([edgesExcel(1) edgesExcel(end)],[1 1],'--','linewidth',2,'Color',[0.5 0.5 0.5]);
% tCommon = [];
% aRatioCommon = tCommon;
% 
% for i = 1:length(keepIndices)
% kI = keepIndices{i};
% rI = removeIndices{i};
% tCommon = [tCommon; tCombined(kI)];
% aRatioCommon = [aRatioCommon; ampsCombined(kI)./ampsCombined(rI)];
% end
% 
% [tCommon,sI2] = sort(tCommon);
% aRatioCommon = aRatioCommon(sI2);
% figure(); semilogy(tCommon,aRatioCommon,'.'); zoom on; grid on;
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

