%%
clear;
close all; clc;
%cd ~/research/now/reventador/

%load ~/research/now/reventador/ReveMagnitudes_v1.mat;
%load ~/igdata/ReveMagnitudes_v2.mat;
[t,M1,magErr] = getReventadorMagnitudes();
medMag = M1; %medfiltSH(M1,7,true); % - 2;
goodI = medMag >= 0.3 & magErr <= 1/2 & t >= datetime(2016,01,01) & t <= datetime(2019,05,01);
tGoodOrig = t(goodI);
aGoodOrig = medfiltSH(medMag(goodI),1,true); %10.^(3+1.44*medMag(goodI));
[tGood,aGood] = filterCatalog(tGoodOrig,aGoodOrig,30);

tStart = [datetime(2016,11,10);...
    datetime(2017,09,23);...
    datetime(2018,04,02);...
    datetime(2018,07,11);...
    datetime(2018,09,19);...
    datetime(2018,11,24)];

tEnd = [datetime(2017,09,22);...
    datetime(2018,04,01);...
    datetime(2018,07,10);...
    datetime(2018,09,18);...
    datetime(2018,11,23);...
    datetime(2019,03,27)];

%%
% load ~/igdata/reventadorExcelData.mat;
% refTime = min([min(tGood) min(tsheet)]); 
% [lia,locb] = ismembertol(seconds(tsheet-refTime),seconds(tGood-refTime),10,'DataScale',1);
% locb = locb(lia);
% [lia2,locb2] = ismember((1:length(tGood))',locb);
% allRegionalI = (1:length(tGood))';
% badRegionalI = allRegionalI(~lia2);
% figure();
% ax(1) = subplot(311);
% semilogy(tsheet(lia),amp(lia),'.'); zoom on; grid on; hold on;
% semilogy(tsheet(~lia),amp(~lia),'o'); zoom on; grid on;
% 
% ax(2)= subplot(312);
% semilogy(tGood(locb),max(aGood(locb),[],2),'.');
% 
% ax(3)= subplot(313);
% semilogy(tGood(badRegionalI),max(aGood(badRegionalI),[],2),'.');
% linkaxes(ax,'x');

%%
% figure();
% plot(tGood(locb),(1:length(locb))','.'); zoom on; grid on;
% hold on;
% plot(tGood(badRegionalI),(1:length(badRegionalI))','.');

%%
nDays = 14;
%[rate,meanMag] = t2r(tGood,days(nDays),aGood);
[rate,meanMag,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = t2r(tGood,days(nDays),aGood);
figure();
axx(1) = subplot(311);
semilogy(tGood,rate/nDays,'.'); grid on;

axx(2) = subplot(312);
plot(tGood,meanMag,'.'); grid on;

energySmoothed = (10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays;

axx(3) = subplot(313); semilogy(tGood,energySmoothed,'.'); zoom on; grid on;

linkaxes(axx,'x');

%%
juvenileContent = [91.45;...
    70.10;...
    81.33;...
    85.04;...
    86.48;...
    77.53];

frp = [10.46039578;...
    9.905761317;...
    4.290847458;...
    4.194;...
    4.989756098;...
    6.20296875];

columnHeights = [1274.625;...
    1266.066869;...
    1378.042254;...
    1445.755814;...
    1604.569106;...
    1251.266094];

totalDense = [96.30;...
    79.93;...
    80.93;...
    78.50;...
    75.64;...
    97.71];

totalVesicular = [3.70;...
    20.07;...
    19.07;...
    21.50;...
    24.36;...
    2.29];

convexityXsolidity = [0.78;...
0.78;...
0.79;...
0.79;...
0.76;...
0.78];

lstarts = length(tStart);
figure('units','normalized','outerposition',[0 0 3/4 1]);
axF = gobjects(lstarts,1);
medianAmps = NaN(lstarts,1);
averageRate = medianAmps;
averageEnergy = medianAmps;
durations = NaT(lstarts,1) - NaT(1);
totEnergy = medianAmps;
for i = 1:lstarts
    ia = tGood >= tStart(i) & tGood <= tEnd(i);
    axF(i) = subplot(lstarts,1,i);
    thisT = tGood(ia);
    thisA = aGood(ia);
    thisE = energySmoothed(ia);
    plot(thisT,thisA,'.'); zoom on; grid on;
    %meanAmp_ = mean(thisA);

    medianAmp_ = median(thisA,"omitnan");
    sumE = sum(thisE);
    medianAmps(i) = medianAmp_;
    dur_ = (tEnd(i) - tStart(i));
    durations(i) = dur_;
    avgRate_ = sum(ia)/days(dur_);
    averageRate(i) = avgRate_;
    totEnergy(i) = sumE;
    avgEnergy_ = median(thisE,"omitnan"); %sumE/days(dur_);
    averageEnergy(i) = avgEnergy_;

    titleStr = sprintf('Median Amplitude: %f; Duration: %d; Avg. Rate: %f; Avg. Energy: %f;',...
        medianAmp_,days(dur_),avgRate_,avgEnergy_);
    title(titleStr);
    axis(axF(i),'tight');
end
linkaxes(axF,'y');

%%
symbolSize = 50;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax3(1) = subplot(411);
plot(tEnd,medianAmps,'.','MarkerSize',symbolSize);
%ylim([2000 9000]);
zoom on; grid on;

ax3(2) = subplot(412);
semilogy(tEnd,averageRate,'.','MarkerSize',symbolSize);
ylim([10 100]);
zoom on; grid on;

ax3(3) = subplot(413);
semilogy(tEnd,averageEnergy,'.','MarkerSize',symbolSize);
zoom on; grid on;

ax3(4) = subplot(414);
semilogy(tEnd,juvenileContent,'.','MarkerSize',symbolSize);
zoom on; grid on;

linkaxes(ax3,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax4(1) = subplot(311);
SS(1) = scatter(juvenileContent,medianAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(1).YScale = 'lin';
%ylim([2000 9000]);
zoom on; grid on;
xlabel('Average Juvenile Content');

ax4(2) = subplot(312);
SS(2) = scatter(juvenileContent,averageRate,symbolSize,datenum(tEnd),'filled');
SS(2).MarkerEdgeColor = 'k';
SS(2).MarkerEdgeAlpha = 0.5;
SS(2).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(2).YScale = 'log';
ylim([10 100]);
zoom on; grid on;
xlabel('Average Juvenile Content');

ax4(3) = subplot(313);
SS(3) = scatter(juvenileContent,averageEnergy,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(juvenileContent,log10(averageEnergy)));
yq = polyval(b_,juvenileContent);
ll = semilogy(ax4(3),juvenileContent,10.^yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Juvenile Content');
linkaxes(ax4,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax4(1) = subplot(311);
SS(1) = scatter(columnHeights,medianAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(1).YScale = 'lin';
%ylim([2000 9000]);
zoom on; grid on;
xlabel('Average Column Heights');

ax4(2) = subplot(312);
SS(2) = scatter(columnHeights,averageRate,symbolSize,datenum(tEnd),'filled');
SS(2).MarkerEdgeColor = 'k';
SS(2).MarkerEdgeAlpha = 0.5;
SS(2).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(2).YScale = 'log';
ylim([10 100]);
zoom on; grid on;
xlabel('Average Column Heights');

ax4(3) = subplot(313);
SS(3) = scatter(columnHeights,averageEnergy,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
ax4(3).YScale = 'log';
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(columnHeights,log10(averageEnergy)));
yq = polyval(b_,columnHeights);
ll = semilogy(ax4(3),columnHeights,10.^yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Column Heights');
linkaxes(ax4,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax5(1) = subplot(311);
SS(1) = scatter(frp,medianAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax5(1).YScale = 'lin';
%ylim([2000 9000]);
zoom on; grid on;
xlabel('Average Fire Radiative Power');

ax5(2) = subplot(312);
SS(2) = scatter(frp,averageRate,symbolSize,datenum(tEnd),'filled');
SS(2).MarkerEdgeColor = 'k';
SS(2).MarkerEdgeAlpha = 0.5;
SS(2).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax5(2).YScale = 'log';
ylim([10 100]);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(frp,log10(averageRate)));
yq = polyval(b_,frp);
ll = semilogy(ax5(2),frp,10.^yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Fire Radiative Power');

ax5(3) = subplot(313);
SS(3) = scatter(frp,averageEnergy,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(frp,log10(averageEnergy)));
yq = polyval(b_,frp);
ll = semilogy(ax5(3),frp,10.^yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Fire Radiative Power');
linkaxes(ax5,'x');


%%
for i = 1:lstarts
    figure('units','normalized','outerposition',[0 0 3/4 1]);
    ia = tGood >= tStart(i) & tGood <= tEnd(i);

    thisT = tGood(ia);
    thisA = aGood(ia);

    axIET(1) = subplot(211);
    SS_ = scatter(thisA(1:end-1),seconds(diff(thisT)),5*exp((thisA(1:end-1))),datenum(thisT(1:end-1)),'filled');
    axIET(1).YScale = 'log';
    axIET(1).XScale = 'linear';
    Cbar = colorbar; colormap turbo;
    Cbar.TickLabels = datestr(Cbar.Ticks);
    zoom on; grid on;
    xlabel('Amplitud of Current Event');
    ylabel('Waiting Time to Next Event');
    %xlim([2e3 8e9]);
    ylim([1e1 1e5]);
    SS_.MarkerEdgeColor = 'k';
    SS_.MarkerEdgeAlpha = 0.5;
    SS_.MarkerFaceAlpha = 0.5;

    hold on;
    b_ = flipud(robustfit((thisA(1:end-1)),log10(seconds(diff(thisT)))));
    yq = polyval(b_,[floor((min(thisA))) ceil((max(thisA)))]);
    ll = semilogy([floor((min(thisA))) ceil((max(thisA)))],10.^yq,'linewidth',3);
    ll.Color(4) = 0.5;
    title([datestr(tStart(i)),'-',datestr(tEnd(i))])

    axIET(2) = subplot(212);
    SS_ = scatter(seconds(diff(thisT)),thisA(2:end),5*exp((thisA(2:end))),datenum(thisT(2:end)),'filled');
    axIET(2).YScale = 'linear';
    axIET(2).XScale = 'log';
    zoom on; grid on;
    xlabel('Time Since Last Event');
    ylabel('Amplitude of Next Event');
    xlim([1e1 1e5]);
    %ylim([2e3 8e9]);
    SS_.MarkerEdgeColor = 'k';
    SS_.MarkerEdgeAlpha = 0.5;
    SS_.MarkerFaceAlpha = 0.5;
    Cbar = colorbar; colormap turbo;
    Cbar.TickLabels = datestr(Cbar.Ticks);

    hold on;
    b_ = flipud(robustfit(log10(seconds(diff(thisT))),(thisA(2:end))));
    yq = polyval(b_,[floor(log10(min(seconds(diff(thisT))))) ceil(log10(max(seconds(diff(thisT)))))]);
    ll = semilogx(10.^[floor(log10(min(seconds(diff(thisT))))) ceil(log10(max(seconds(diff(thisT)))))],yq,'linewidth',3);
    ll.Color(4) = 0.5;
end
