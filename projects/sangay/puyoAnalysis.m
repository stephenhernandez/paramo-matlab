clear; close all; clc;
cd ~/research/now/sangay/
load puyo_match_filter_analysis.mat

%%
minPKS = 0.6;
tabs = tabs(pksOrig >= minPKS);
maxAmpRMS = maxAmpRMS(pksOrig >= minPKS);

%%
load modis.txt
tmodis = datetime(modis(:,3)+2000,modis(:,2),modis(:,1),modis(:,4),modis(:,5),zeros(length(modis),1));
rf = modis(:,6); %radiant flux
crf = modis(:,7); %cumulative rf
cv = modis(:,8); %cumulative volume
tI = tabs >= min(tmodis);% & tabs <= max(tmodis);
cvInterpolated = interp1(datenum(tmodis)-datenum(min(tmodis)),cv,datenum(tabs(tI))-datenum(min(tmodis)),'linear','extrap');
%figure(); plot(tabs(tI),cvInterpolated,'o'); zoom on;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs(1:end-1),medfiltSH(86400*medfiltSH(1./seconds(diff(tabs)),31,false),31),'.'); zoom on;
ax = gca; ax.YScale = 'log';
ylabel('estimated instantaneous event rate [\#/day]');

figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs,maxAmpRMS,'.'); zoom on;
ax = gca; ax.YScale = 'log';

figure('units','normalized','outerposition',[0 0 1 1]);
plot(tabs,1:length(tabs),'o'); zoom on;

%%
puyoCatalogStart = datetime(2012,06,29);
puyoCatalogEnd = datetime(2019,11,29);
episodes = [puyoCatalogStart datetime(2013,05,29);...   % tail end of continuous phase
    datetime(2015,01,02) datetime(2015,04,25);...       % 2015
    datetime(2016,03,05) datetime(2016,07,31);...       % 2016
    datetime(2017,07,21) datetime(2017,10,29);...       % 2017
    datetime(2018,08,11) datetime(2018,12,05);...       % 2018
    datetime(2019,05,01) puyoCatalogEnd];               % 2019

nEpisodes = size(episodes,1);
episodeDurations = days(episodes(:,2) - episodes(:,1));
pseudoMag = (1)*log10(maxAmpRMS);
pseudoEnergy = mw2m0(pseudoMag);
cumulativeSeismicEnergyPerEpisode = NaN(nEpisodes,1);
totalRadiantFlux = cumulativeSeismicEnergyPerEpisode;
nEventsPerEpisode = cumulativeSeismicEnergyPerEpisode;
lStr = ["pre-2015 (incomplete)" "2015" "2016" "2017" "2018" "2019 (on-going)"];

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,3);
ax(3) = subplot(2,2,[2 4]);
for i = 1:nEpisodes
    tI = tabs >= episodes(i,1) & tabs <= episodes(i,2);
    nEventsPerEpisode(i) = sum(tI);
    cumulativeSeismicEnergyPerEpisode(i) = sum(pseudoEnergy(tI));
    
    %disp([nEventsPerEpisode(i) episodeDurations(i) round(nEventsPerEpisode(i)/episodeDurations(i)) ...
    %    log10(cumulativeSeismicEnergyPerEpisode(i)) 1e-6*round(max(cvInterpolated(tI) - min(cvInterpolated(tI))))]);
    
    tdum = datenum(tabs(tI));
    tdum = tdum - min(tdum);
    tNorm = tdum / max(tdum);
    hold(ax(1),'on');
    hold(ax(2),'on');
    hold(ax(3),'on');
    plot(ax(1),tNorm,cumsum(pseudoEnergy(tI)),'.','linewidth',4);
    semilogy(ax(2),tNorm,1e-6*(cvInterpolated(tI)-min(cvInterpolated(tI))),'.','linewidth',4);
    plot(ax(3),tdum,1:sum(tI),'.','linewidth',4);
    
    tI = tmodis >= episodes(i,1) & tmodis <= episodes(i,2);
    totalRadiantFlux(i) = max(crf(tI)-min(crf(tI)));
    
    disp([nEventsPerEpisode(i) episodeDurations(i) round(nEventsPerEpisode(i)/episodeDurations(i)) ...
        log10(cumulativeSeismicEnergyPerEpisode(i)) ...
        1e-6*round(max(cvInterpolated(tI) - min(cvInterpolated(tI)))) ...
        totalRadiantFlux(i)]);
end
legend(ax(1),lStr,'location','northwest'); %,'NumColumns',2);
legend(ax(2),lStr,'location','northwest'); %,'NumColumns',2);
legend(ax(3),lStr,'location','southeast'); %,'NumColumns',2);
xlabel(ax(1),'normalized time');
ylabel(ax(1),'cumulative pseudo energy');
xlabel(ax(2),'normalized time');
ylabel(ax(2),'cumulative volume [$Mm^3$]');
xlabel(ax(3),'days since eruption start');
ylabel(ax(3),'cumulative events (PUYO)');
zoom on;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax(1) = subplot(211);
ax(2) = subplot(212);
for i = 1:nEpisodes
    tI = tabs >= episodes(i,1) & tabs <= episodes(i,2);
    tdum = datenum(tabs(tI));
    tdum = tdum - min(tdum);
    hold(ax(1),'on');
    hold(ax(2),'on');
    S = plot(ax(1),tdum,pseudoMag(tI),'.','linewidth',4);
    [f,xi] = ksdensity(pseudoMag(tI),(-0.4:0.01:1.4)');
    plot(ax(2),xi,f,'-','linewidth',4);
end
legend(ax(1),lStr,'location','southeast'); %,'NumColumns',2);
legend(ax(2),lStr,'location','northwest'); %,'NumColumns',2);
xlabel(ax(1),'days since eruption start'); %);
ylabel(ax(1),'pseudo magnitude');
title(ax(1),'magnitude time evolution');
title(ax(2),'magnitude probability distributions');
ylabel(ax(2),'$\propto$ P');
xlabel(ax(2),'pseudo magnitude');
zoom on;

%%
lStr = ["pre-2015 (incomplete)" "2015" "2016" "2017" "2018" "2019 (til 06 Nov.)"];
figure('units','normalized','outerposition',[0 0 1 1]);
ax = gca;
for i = 1:nEpisodes
    tI = tmodis >= episodes(i,1) & tmodis <= episodes(i,2);
    tdum = datenum(tmodis(tI));
    tdum = tdum - min(tdum);
    hold(ax(1),'on');
    S = plot(ax(1),tdum,crf(tI)-min(crf(tI)),'.-','linewidth',2);
end
legend(ax(1),lStr,'location','northeast'); %,'NumColumns',2);
xlabel(ax(1),'days since eruption start');
ylabel(ax(1),'cumulative radiant flux [MW]');
zoom on;


%%
close all;
lStr = ["pre-2015 (incomplete)" "2015" "2016" "2017" "2018" "2019"];
figure('units','normalized','outerposition',[0 0 1 1]);
ax = gca;
for i = 2:nEpisodes
    figure('units','normalized','outerposition',[0 0 1 1]);
    ax = gca;
    tI = tabs >= episodes(i,1) & tabs <= episodes(i,2);
    histogram(tabs(tI),dateshift(min(tabs(tI)),'start','day'):dateshift(max(tabs(tI)),'end','day'));
    %legend(ax(1),lStr(i),'location','northeast'); %,'NumColumns',2);
    %xlabel(ax(1),'days since eruption start');
    ylabel(ax(1),'Eventos / Dia');
    zoom on;
end
