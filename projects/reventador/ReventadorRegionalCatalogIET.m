clear; close all; clc;

cd ~/research/now/reventador/

[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v10',10); %,10);

%%
tabs = tabs + seconds(10);

%
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

% casc
z2p(:,1:3) = d_(1)*d_(1)*z2p(:,1:3)/sensitivities(1);

% boni
boniI = tabs <= datetime(2020,01,339);
z2p(boniI,4:6) = d_(2)*d_(2)*z2p(boniI,4:6)/sensitivities(2);
z2p(~boniI,4:6) = d_(2)*d_(2)*z2p(~boniI,4:6)/sensitivities(3);

% ants
antsI = tabs <= datetime(2014,08,14);
z2p(antsI,7:9) = d_(3)*d_(3)*z2p(antsI,7:9)/sensitivities(2);
z2p(~antsI,7:9) = d_(3)*d_(3)*z2p(~antsI,7:9)/sensitivities(3);

% antg
z2p(:,10:12) = d_(4)*d_(4)*z2p(:,10:12)/sensitivities(1);

%%
z2p = median([median(z2p(:,1:3),2,'omitnan') median(z2p(:,4:6),2,'omitnan') median(z2p(:,7:9),2,'omitnan') median(z2p(:,10:12),2,'omitnan')],2,'omitnan');

%%
minAmp = 300; %1e-2;
maxAmp = 1e5;
winlen = 51;
nboot = 5e2;
zeroPhaseFlag = false;

%%
goodI = z2p >= minAmp & z2p <= maxAmp & NCC >= 0.35 & ...
    median(p2rms,2,'omitnan') >= 3 & max(p2rms,[],2,'omitnan') <= 7 & ...
    min(p2rms,[],2,'omitnan') >= 2.5 & median(kurt,2,'omitnan') >= 3 & ...
    max(kurt,[],2,'omitnan') <= 20 & min(kurt,[],2,'omitnan') >= 2.5;

%%
tGood = tabs(goodI);
aGood = z2p(goodI);

%%
load reventadorExcelData.mat;
refTime = min([min(tGood) min(tsheet)]); 
[lia,locb] = ismembertol(seconds(tsheet-refTime),seconds(tGood-refTime),10,'DataScale',1);
locb = locb(lia);
[lia2,locb2] = ismember((1:length(tGood))',locb);
allRegionalI = (1:length(tGood))';
badRegionalI = allRegionalI(~lia2);

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
nDays = 3;
[rate,meanamp] = t2r(tGood,days(nDays),aGood);
figure();
axx(1) = subplot(311);
semilogy(tGood,rate/nDays,'.'); grid on;

axx(2) = subplot(312);
semilogy(tGood,meanamp,'.'); grid on;

axx(3) = subplot(313); semilogy(tGood,rate.*meanamp/nDays,'.'); zoom on; grid on;

linkaxes(axx,'x');

%%
%close all;
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
meanAmps = NaN(lstarts,1);
averageRate = meanAmps;
durations = NaT(lstarts,1) - NaT(1);
for i = 1:lstarts
    ia = tGood >= tStart(i) & tGood <= tEnd(i);
    axF(i) = subplot(lstarts,1,i);
    thisT = tGood(ia);
    thisA = aGood(ia);
    semilogy(thisT,thisA,'o');
    %meanAmp_ = mean(thisA);
    meanAmp_ = nanmedian(thisA);

    meanAmps(i) = meanAmp_;
    durations(i) = tEnd(i) - tStart(i);
    averageRate(i) = sum(ia)/days(durations(i));
    title(strcat('Mean Amplitude: ',num2str(meanAmp_),...
        ', Duration: ',num2str(days(durations(i))),...
        ', Average Rate: ',num2str(averageRate(i))));
    axis(axF(i),'tight');
end
linkaxes(axF,'y');

%%
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax3(1) = subplot(411);
semilogy(tEnd,meanAmps,'.');
ylim([400 4000]);
zoom on; grid on;

ax3(2) = subplot(412);
semilogy(tEnd,averageRate,'.');
ylim([10 100]);
zoom on; grid on;

ax3(3) = subplot(413);
semilogy(tEnd,averageRate.*meanAmps,'.');
zoom on; grid on;

ax3(4) = subplot(414);
semilogy(tEnd,juvenileContent,'.');
zoom on; grid on;

linkaxes(ax3,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax4(1) = subplot(311);
SS(1) = scatter(juvenileContent,meanAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(1).YScale = 'log';
ylim([400 4000]);
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
SS(3) = scatter(juvenileContent,averageRate.*meanAmps,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(juvenileContent,(averageRate.*meanAmps)));
yq = polyval(b_,juvenileContent);
ll = semilogy(ax4(3),juvenileContent,yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Juvenile Content');
linkaxes(ax4,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax4(1) = subplot(311);
SS(1) = scatter(columnHeights,meanAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax4(1).YScale = 'log';
ylim([400 4000]);
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
SS(3) = scatter(columnHeights,averageRate.*meanAmps,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(columnHeights,(averageRate.*meanAmps)));
yq = polyval(b_,columnHeights);
ll = semilogy(ax4(3),columnHeights,yq,'linewidth',3);
ll.Color(4) = 0.5;
xlabel('Average Column Heights');
linkaxes(ax4,'x');

%%
symbolSize = 150;
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax5(1) = subplot(311);
SS(1) = scatter(frp,meanAmps,symbolSize,datenum(tEnd),'filled');
SS(1).MarkerEdgeColor = 'k';
SS(1).MarkerEdgeAlpha = 0.5;
SS(1).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
ax5(1).YScale = 'log';
ylim([400 4000]);
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
xlabel('Average Fire Radiative Power');

ax5(3) = subplot(313);
SS(3) = scatter(frp,averageRate.*meanAmps,symbolSize,datenum(tEnd),'filled');
SS(3).MarkerEdgeColor = 'k';
SS(3).MarkerEdgeAlpha = 0.5;
SS(3).MarkerFaceAlpha = 0.5;
Cbar = colorbar; colormap turbo;
Cbar.TickLabels = datestr(Cbar.Ticks);
zoom on; grid on;
hold on;
b_ = flipud(robustfit(frp,(averageRate.*meanAmps)));
yq = polyval(b_,frp);
ll = semilogy(ax5(3),frp,yq,'linewidth',3);
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
    SS_ = scatter(thisA(1:end-1),seconds(diff(thisT)),5*exp(log10(thisA(1:end-1))),datenum(thisT(1:end-1)),'filled');
    axIET(1).YScale = 'log';
    axIET(1).XScale = 'log';
    Cbar = colorbar; colormap turbo;
    Cbar.TickLabels = datestr(Cbar.Ticks);
    zoom on; grid on;
    xlabel('Amplitud of Current Event');
    ylabel('Waiting Time to Next Event');
    xlim([1e2 1e5]);
    ylim([1e1 1e5]);
    SS_.MarkerEdgeColor = 'k';
    SS_.MarkerEdgeAlpha = 0.5;
    SS_.MarkerFaceAlpha = 0.5;

    hold on;
    b_ = flipud(robustfit(log10(thisA(1:end-1)),log10(seconds(diff(thisT)))));
    yq = polyval(b_,[floor(log10(min(thisA))) ceil(log10(max(thisA)))]);
    ll = loglog(10.^[floor(log10(min(thisA))) ceil(log10(max(thisA)))],10.^yq,'linewidth',3);
    ll.Color(4) = 0.5;

    axIET(2) = subplot(212);
    SS_ = scatter(seconds(diff(thisT)),thisA(2:end),5*exp(log10(thisA(2:end))),datenum(thisT(2:end)),'filled');
    axIET(2).YScale = 'log';
    axIET(2).XScale = 'log';
    zoom on; grid on;
    xlabel('Time Since Last Event');
    ylabel('Amplitude of Next Event');
    xlim([1e1 1e5]);
    ylim([1e2 1e5]);
    SS_.MarkerEdgeColor = 'k';
    SS_.MarkerEdgeAlpha = 0.5;
    SS_.MarkerFaceAlpha = 0.5;
    Cbar = colorbar; colormap turbo;
    Cbar.TickLabels = datestr(Cbar.Ticks);

    hold on;
    b_ = flipud(robustfit(log10(seconds(diff(thisT))),log10(thisA(2:end))));
    yq = polyval(b_,[floor(log10(min(seconds(diff(thisT))))) ceil(log10(max(seconds(diff(thisT)))))]);
    ll = loglog(10.^[floor(log10(min(seconds(diff(thisT))))) ceil(log10(max(seconds(diff(thisT)))))],10.^yq,'linewidth',3);
    ll.Color(4) = 0.5;
end


