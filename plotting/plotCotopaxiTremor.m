clear; close all;
matDir = "~/masa/old/research/now/cotopaxi/";
cd(matDir);
load(fullfile(matDir,"CotopaxiTremorLocationResults_v4"));

Mcorr = McorrOrig;
mcorrI = ~isfinite(Mcorr);
nBad = sum(mcorrI,2);
magLogThresh = log10(400);
magErrorThresh = 5e-3;
smoothN = 5;
minN = 5;
snrThresh = 8;

%% figure 1
mI = magErr <= magErrorThresh & M >= magLogThresh;
timetable(unique(dateshift(t(mI),'start','day')),groupcounts(dateshift(t(mI),'start','day')))
[~,maxMI] = max(abs(Mcorr-M),[],2,"omitnan");
for i = 1:size(Mcorr,1)
    if nBad(i) < 2
        Mcorr(i,maxMI(i)) = NaN;
    end
end
mcorrI = ~isfinite(Mcorr);
nBad = sum(mcorrI,2);
M = median(Mcorr,2,"omitnan");
magErr2 = mad(Mcorr,1,2);

mI2 = magErr2 <= magErrorThresh & M >= magLogThresh;
timetable(unique(dateshift(t(mI2),'start','day')),groupcounts(dateshift(t(mI2),'start','day')))
[~,maxMI] = max(abs(Mcorr-M),[],2,"omitnan");
for i = 1:size(Mcorr,1)
    if nBad(i) < 2
        Mcorr(i,maxMI(i)) = NaN;
    end
end

mcorrI = ~isfinite(Mcorr);
nBad = sum(mcorrI,2);
M = median(Mcorr,2,"omitnan");
magErr2 = mad(Mcorr,1,2);
magErr3 = std(Mcorr,0,2,"omitnan");
snr = median(Mcorr,2,"omitnan")./range(Mcorr,2);

% figure(1);
% hold on;
% semilogy(t,magErr2,'.'); zoom on; grid on;

maxMonth = datetime(2023,08,01);
mag2smooth = medfiltSH(magErr2,smoothN,true);
mag3smooth = medfiltSH(magErr3,smoothN,true);
metric1 = mag2smooth.*mag3smooth;
mI2 = M >= magLogThresh & ...
    7-nBad>= minN & ...
    metric1 <= magErrorThresh & ...
    mag3smooth <= 0.1 & ...
    mag2smooth <= 5e-2 & ...
    magErr3 <= 0.10 & ...
    snr >= snrThresh & t <= maxMonth;

timetable(unique(dateshift(t(mI2),'start','day')),groupcounts(dateshift(t(mI2),'start','day')))
t2 = t(mI2);
tOrig = t;
M2 = 10.^M(mI2);
uniqueDays = unique(dateshift(t2,'start','day'));

lDays = length(uniqueDays); clear sumEnergy;
for i = 1:lDays
    tI = t2 >= uniqueDays(i) & t2 < uniqueDays(i)+1;
    sumEnergy(i,1) = sum((M2(tI)*1e-9).^2);
end

%% figure 2
figure('units','normalized','outerposition',[0 0 1.5/2 1]);
loglog(groupcounts(dateshift(t(mI2),'start','day')),sumEnergy,'o','LineWidth',3);
zoom on; grid on; xlabel("Minutes of Tremor"); ylabel("$\propto$ Energy");

%% figure 3

%stairs([uniqueDays; max(uniqueDays)+1],[groupcounts(dateshift(t(mI2),'start','day'))/60; 0],'-','linewidth',3);
t3 = (dateshift(min(t2),'start','day'):dateshift(max(t2),'end','day'))';
y3 = zeros(size(t3));
[lia,locb] = ismember(uniqueDays,t3);
y3(locb(lia)) = groupcounts(dateshift(t(mI2),'start','day'))/60;
figure('units','normalized','outerposition',[0 0 3/4 1]);
stairs(t3,y3,'-','linewidth',2); zoom on; grid on;
zoom on; grid on; ylabel("Hours of Tremor Per Day"); ylim([0 24]); %axis tight;

%% figure 4
y4 = zeros(size(t3));
[lia,locb] = ismember(uniqueDays,t3);
y4(locb(lia)) = 1e9*sqrt(sumEnergy./groupcounts(dateshift(t2,'start','day')));
figure('units','normalized','outerposition',[0 0 3/4 1]);
stairs(t3,y4,'-','linewidth',2); zoom on; grid on;
zoom on; grid on; title("RSAM Amplitude at 1 km.");

%% figure 5
nDays = 1;
[rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = ...
    t2r(t2,days(nDays),log10(M2));

figure('units','normalized','outerposition',[0 0 3/4 1]);
plot(t2,rate/nDays/60,'.'); zoom on; grid on;
title("Hours of Tremor Per Day");

%% figure 6
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(t2,10.^meanMagsFixedTimeWin,'.'); zoom on; grid on;
ylim([50 1e4]);
title("Mean RSAM Amplitude at 1 km. per 24 hours");

%% figure 7
figure('units','normalized','outerposition',[0 0 3/4 1]);
loglog(groupcounts(dateshift(t(mI2),'start','day')),1e9*sqrt(sumEnergy./groupcounts(dateshift(t(mI2),'start','day'))),'o','linewidth',3);
zoom on; grid on;
xlabel("Minutes of Tremor");
ylabel("Average RSAM at 1 km");

%% figure 8
figure('units','normalized','outerposition',[0 0 3/4 1]);
plot(groupcounts(dateshift(t(mI2),'start','day')),1e9*sqrt(sumEnergy./groupcounts(dateshift(t(mI2),'start','day'))),'o','linewidth',3);
zoom on; grid on; xlabel("Minutes of Tremor");
ylabel("Average RSAM at 1 km");

%% figure 9
% figure('units','normalized','outerposition',[0 0 1.5/2 1]);
% semilogy(t2,1:length(t2),'.');
% zoom on; grid on;
% ylabel("Cumulative Minutes of Tremor");

figure('units','normalized','outerposition',[0 0 1.5/2 1]);
plot((1:length(t2))',t2,'.');
zoom on; grid on;
xlabel("Cycle Number");

%%
gaps = minutes(diff(t2));
minGap = 10;
minDur = 5;

gapI = gaps >= minGap;
episodeStarts = [t2(1)-minutes(1); t2(find(gapI)+1)-minutes(1)];
episodeEnds = [episodeStarts(2:end) - minutes(gaps(gapI))+minutes(1); t2(end)];
quiescence = episodeStarts(2:end) - episodeEnds(1:end-1);
episodeDurations = episodeEnds - episodeStarts;

durI = episodeDurations > minutes(minDur);
eStarts = episodeStarts(durI);
eEnds = episodeEnds(durI);
eDurs = episodeDurations(durI);
reposeTimes = eStarts(2:end) - eEnds(1:end-1);
eStarts = eStarts(1:end-1);
eEnds = eEnds(1:end-1);
eDurs = eDurs(1:end-1);
%timetable(eStarts,eEnds,eDurs,reposeTimes)

%% figure 10
% figure('units','normalized','outerposition',[0 0 3/4 1]);
% plot(episodeStarts,1:length(episodeStarts),'.');
% zoom on; grid on; hold on; plot(episodeEnds,1:length(episodeEnds),'.');
% ylabel("Cycle Number");
% legend("Tremor Start","Tremor End","Location","Best");

figure('units','normalized','outerposition',[0 0 3/4 1]);
plot((1:length(episodeStarts))',episodeStarts,'.');
zoom on; grid on; hold on; plot((1:length(episodeEnds)),episodeEnds,'.');
xlabel("Cycle Number");
legend("Tremor Start","Tremor End","Location","Best");

%% figure 11
figure('units','normalized','outerposition',[0 0 3/4 1]);
plot(sort(minutes(eDurs)),1:length(eDurs),'.'); zoom on; grid on; hold on;
plot(sort(minutes(reposeTimes+eDurs)),1:length(reposeTimes),'.'); zoom on; grid on;
plot(sort(minutes(reposeTimes)),1:length(reposeTimes),'.'); zoom on; grid on;
legend("Tremor Episode Durations","Tremor Cycle Durations","Repose Durations","Location","NorthWest");
ax = gca;
ax.XScale = 'log';

%% figure 12
figure('units','normalized','outerposition',[0 0 3/4 1]);
plot(eEnds,minutes(reposeTimes)./minutes(reposeTimes+eDurs),'.');
zoom on; grid on; title("Repose Duration / Total Cycle Duration");
hold on;
plot(eEnds,medfiltSH(minutes(reposeTimes)./minutes(reposeTimes+eDurs),smoothN,true),'.');

%% figure 13
figure('units','normalized','outerposition',[0 0 3/4 1]);
clear ax; ax(1) = subplot(211); hold on;
for i = 1:length(eStarts)
    aa = area([eStarts(i) eEnds(i)],minutes(reposeTimes(i)+eDurs(i))*[1 1]);
    aa.FaceAlpha = 0.3;
    ax(1).ColorOrderIndex = 1;
    aa.EdgeColor = 'none';
end
ax(1).ColorOrderIndex = 2;
for i = 1:length(eStarts)
    aa = area([eEnds(i) eEnds(i)+reposeTimes(i)],minutes(reposeTimes(i)+eDurs(i))*[1 1]);
    aa.FaceAlpha = 0.3;
    ax(1).ColorOrderIndex = 2;
    aa.EdgeColor = 'none';
end

stairs(eStarts,minutes(reposeTimes+eDurs),'k','linewidth',2); grid on;
ax(2) = subplot(212); semilogy(t2,M2,'.'); zoom on; grid on;
linkaxes(ax,'x');

ax(2).ColorOrderIndex = 2; hold(ax(2),'on');
for i = 1:length(eStarts)
    aa = area(ax(2),[eEnds(i) eEnds(i)+reposeTimes(i)],ceil(max(M2))*[1 1]);
    aa.FaceAlpha = 0.3;
    ax(2).ColorOrderIndex = 2;
    aa.EdgeColor = 'none';
end
ax(1).YScale = 'log'; ylim(ax(1),[minGap ceil(max(minutes(reposeTimes+eDurs)))]);

%% figure 14
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(t2,cumsum(M2.^2),'.'); zoom on; grid on;
ylabel("$\propto$ Cumulative Energy");

%% figure 15
figure('units','normalized','outerposition',[0 0 3/4 1]);
snrI = snr >= snrThresh;
semilogy(t(snrI),metric1(snrI),'.','linewidth',2);
zoom on; grid on; hold on;
semilogy(t(snrI),[mag2smooth(snrI) mag3smooth(snrI)],'.','linewidth',2);

%% figure 16
figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t(mI2),10.^Mcorr(mI2,:),'.'); zoom on; grid on;
hold on; t(~mI2) = NaT; M(~mI2) = NaN;
semilogy(t,10.^M,'k.-','linewidth',3,'markersize',14);
title("RSAM at 1 km. from source");
ylabel("[$nm \cdot s^{-1}$]");

%% figure 17
figure('units','normalized','outerposition',[0 0 3/4 1]);
SS = scatter(minutes(reposeTimes(1:end-1)),minutes(eDurs(2:end)),2*exp(log10(minutes(reposeTimes(1:end-1)))),datenum(eEnds(2:end)),'filled');
colorbar; zoom on; grid on;
ax = gca; ax.YScale = 'log'; ax.XScale = 'log';
xlabel("Current Repose Time [min.]");
ylabel("Next Episode Duration");

%% figure 18
figure('units','normalized','outerposition',[0 0 3/4 1]);
BB = bar(uniqueDays+0.5,sumEnergy,1); zoom on; grid on;
ylabel("$\propto$ Energy"); title("Daily Energy");
BB.EdgeColor = 'k'; BB.EdgeAlpha = 0.8; BB.FaceAlpha = 0.5;

%% figure 19
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(eEnds,cumsum(minutes(eDurs)),'.'); zoom on; grid on; hold on;
semilogy(eStarts+eDurs+reposeTimes,cumsum(minutes(reposeTimes)),'.');
semilogy(eStarts+eDurs+reposeTimes,cumsum(minutes(eDurs+reposeTimes)),'.');
legend('Episode Durations','Repose Times','Total Cycle Durations','location','best');
ylabel("Cumulative Minutes");

%%
TremorStartTime = eStarts;
TremorEndTime = eEnds;
TremorDuration = eDurs;
ReposeTime = reposeTimes;
TotalCycleDuration = diff(episodeStarts(durI));
timetable(TremorStartTime,TremorEndTime,TremorDuration,ReposeTime,TotalCycleDuration)

%%
maxAmp = NaN(size(uniqueDays));
meanAmp = maxAmp;
medianAmp = maxAmp;
stdAmp = maxAmp;
nMinutes = maxAmp;
sumEnergy = maxAmp;
for i = 1:length(maxAmp)
    tI = t2 >= uniqueDays(i) & t2 < uniqueDays(i)+1;
    nMinutes(i) = sum(tI);
    M2_ = M2(tI);
    sumEnergy(i) = sum((1e-9*M2_).^2);
    maxAmp(i) = max(M2_);
    meanAmp(i) = mean(M2_);
    medianAmp(i) = median(M2_);
    stdAmp(i) = std(M2_);
end
T = table(uniqueDays,nMinutes,sumEnergy,maxAmp,meanAmp,medianAmp,stdAmp);
writetable(T,"~/Desktop/Cotopaxi_2022_2023_dailytremor.csv");

close all;
axL = linkedPlot(uniqueDays,...
    [nMinutes,sumEnergy,maxAmp,meanAmp,medianAmp,stdAmp],'.-');

titles = ["Number of Minutos";"$\propto$ Energy";"Max. Amp.";"Mean Amp.";"Median Amp.";"Amp. Standard Deviation"];
for i = 1:length(axL)
    axL(i).Title.String = titles(i);
end