close all;
dataMain = [];
load ~/igdata/ec_boundaries.mat;
cd ~/research/now/plata/laplataregiontemplatesfamiliesfromgrowclust20112016/;
files = dir('group_*.txt');
lFiles = length(files);
for i = 1:lFiles
    tmpFile = files(i).name;
    if length(tmpFile) == 18
        famNumber = str2double(tmpFile(13:14));
    elseif length(tmpFile) == 17
        famNumber = str2double(tmpFile(13));
    end
    groupNumber = str2double(tmpFile(8));
    data = importdata(tmpFile);
    data = [data repmat(groupNumber,size(data,1),1) repmat(famNumber,size(data,1),1)];
    dataMain = [dataMain; data];
end

tMain = datetime(dataMain(:,1:6));
%sI = tMain >= datetime(2012,01,01) & tMain <= datetime(2013,01,01);
%dataMain = dataMain(sI,:);

tMain = datetime(dataMain(:,1:6));
[tMain,sI] = sort(tMain);
dataMain = dataMain(sI,:);
lon = dataMain(:,7);
lat = dataMain(:,8);
depth = dataMain(:,9);
groupNumber = dataMain(:,10);
familyNumber = dataMain(:,11);

badI = seconds(diff(tMain)) == 0;
dataMain(badI,:) = [];
tMain = datetime(dataMain(:,1:6));
lon = dataMain(:,7);
lat = dataMain(:,8);
depth = dataMain(:,9);
groupNumber = dataMain(:,10);
familyNumber = dataMain(:,11);
boundaryBox = [min(lon)-0.1 max(lon)+0.1 min(lat)-0.1 max(lat)+0.1];
minLon = boundaryBox(1);
maxLon = boundaryBox(2);
minLat = boundaryBox(3);
maxLat = boundaryBox(4);

%%
for i = 5
    gI = groupNumber == i;
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(8,1,[1 2]);
    plot(tMain(gI),1:sum(gI),'.'); zoom on; title(strcat("Group: ",num2str(i)));
    grid on; ax = gca; ax.LineWidth = 1.5;
    ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

    grid on
    uniqFams = unique(familyNumber(gI));
    nFams = 0;
    legStr = [];
    for j = 1:length(uniqFams)
        nFams = nFams + 1;
        fI = gI & familyNumber == uniqFams(j);
        if nFams <= 7
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'.','linewidth',1); zoom on;
        elseif nFams <= 14
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'p','linewidth',1); zoom on;
        elseif nFams <= 21
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'v','linewidth',1); zoom on;
        elseif nFams <= 28
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'d','linewidth',1); zoom on;
        elseif nFams <= 35
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'*','linewidth',1); zoom on;
        end

        %title(strcat("Group: ",num2str(i)));
        hold on; grid on; axis equal;

        legStr = [legStr; ...
            strcat("Family ",num2str(uniqFams(j)),"(",num2str(sum(fI)),")")];
        axis(boundaryBox);
    end

    hold on; plot(lonEC,latEC,'k-','linewidth',2);
    legend(legStr,'Location','Best');
    geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);
    grid on; ax = gca; ax.LineWidth = 1.5;
    ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
end

%%
%jica = importdata("~/jica_catalog_v1.txt");

cd ~/research/now/plata/
dataMainOrig = load('plata_catalog_10JAN2023');
dataMainOrig = dataMainOrig.dataMainOrig;
tJICA = datetime(dataMainOrig(:,1:6));
[tJICA,sI] = sort(tJICA);
dataMain = dataMainOrig(sI,:);

[tJICA,~,removeIndices,~] = removeRepeatedMatches(tJICA,dataMain(:,11),30);
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    dataMain(rI,:) = [];
end
cc1 = dataMain(:,9);
cc2 = dataMain(:,11);
nUsed = dataMain(:,13);
Iall = dataMain(:,10);
ccI = cc1 >= 0.08 & cc1 <= 1 & cc2 >= 13 & nUsed >= 9;% & t >= datetime(2020,01,01);
tOrig = tJICA;
tJICA = tJICA(ccI);
cc1Good = cc1(ccI);
cc2Good = cc2(ccI);
Igood = Iall(ccI);

% tJICA = datetime(jica(:,1:6));
figure(); plot(tJICA,1:length(tJICA),'.'); zoom on;

gI = groupNumber == 5;

t = tMain(gI);
y = (1:sum(gI))';
lon5 = lon(gI);
lat5 = lat(gI);
dep5 = depth(gI);
fam5 = familyNumber(gI);

commonI = tJICA >= min(t) & tJICA <= max(t);
t2 = tJICA(commonI);
y2 = (1:length(t2))';

[t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(t,t2,45,false);

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t2,y2,'.'); zoom on; grid on; hold on; plot(t2(t2_commonI),y2(t2_commonI),'o');
title("JICA Catalog during period of Alex's Catalog");
legend("All Jica events in this time period",...
    "JICA events common to Alex's Catalog","Location","NorthWest");
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

figure('units','normalized','outerposition',[0 0 1 1]);
plot(tMain,1:length(tMain),'.'); zoom on;
title("Alex's Catalog");
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

[unique(fam5(t1_commonI)) groupcounts(fam5(t1_commonI))]
figure(1);
subplot(8,1,[3 4 5 6 7 8]); plot(lon5(t1_commonI),lat5(t1_commonI),'ko','linewidth',2); zoom on;

load ~/research/now/plata/ispt_waveforms_from_alex_families.mat;
A = A_2;
A2 = A(gI);
clear A_2;

amps = (pull(A,'depmax') - pull(A,'depmin'))*0.5*1e9;
amps2 = amps(gI);

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(tMain,amps,'.'); zoom on; grid on;
hold on; semilogy(t,amps2,'v');
ylim([1e0 1e6]);
title("Amplitude at ISPT");
legend('All Groups','Group 5','location','northwest');
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t,amps2,'.'); zoom on; grid on; hold on; semilogy(t(t1_commonI),amps2(t1_commonI),'d');
title("Group 5 amplitudes at ISPT");
legend("All Group 5 Events","Group 5 Events present in JICA Catalog");
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

%%
phasenet = importdata("Catalog_Phasenet_24575_EXPECTATION_EDT_OT_WT_2.txt");
tPhase = datetime(phasenet(:,1:6));
latPhase = phasenet(:,7);
lonPhase = phasenet(:,8);
depthPhase = phasenet(:,9);
magPhase = phasenet(:,10);

pI = lonPhase >= minLon & lonPhase <= maxLon & latPhase >= minLat & latPhase <= maxLat;
tPhase = tPhase(pI);
latPhase = latPhase(pI);
lonPhase = lonPhase(pI);
depthPhase = depthPhase(pI);
magPhase = magPhase(pI);

figure('units','normalized','outerposition',[0 0 1 1]);
plot(lonPhase,latPhase,'.'); zoom on; grid on; hold on; plot(lonEC,latEC,'k-','linewidth',2); axis equal;
axis([-83 -77 -4 2]);
title("PhaseNet Catalog");
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

phaseMagI = magPhase > 0;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tPhase(phaseMagI),magPhase(phaseMagI),'.'); zoom on; grid on;
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
title("PhaseNet Magnitudes");

commonI2 = tJICA >= min(tPhase) & tJICA <= max(tPhase);
t3 = tJICA(commonI2);
Igood = Igood(commonI2);

[t1_commonI2,t2_commonI2,just_t1I2,just_t2I2] = synchronizeCatalog(tPhase,t3,45,true,false);
hold on;
plot(t3(t2_commonI2),magPhase(t1_commonI2),'ko','linewidth',2); zoom on; grid on;
ylim([0 6]);

for i = 1:length(just_t2I2)
    disp(i)
    i_ = just_t2I2(i);
    t_ = t3(i_);
    plot([t_ t_],[0 6],'m-','linewidth',1);
end

figure(7); hold on;
plot(lonPhase(t1_commonI2),latPhase(t1_commonI2),'ko','linewidth',2); zoom on; grid on;
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
t3 = t3(t2_commonI2);
Igood = Igood(t2_commonI2);
hold on;
famIndex = NaN(size(t3));
legStr = [];
legStr2 = legStr;
IgoodUniq = unique(Igood);

lonPhase2 = lonPhase(t1_commonI2);
latPhase2 = latPhase(t1_commonI2);

figure('units','normalized','outerposition',[0 0 1 1]);
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(IgoodUniq)
    figure(9);
    hold on;
    nI = Igood == IgoodUniq(i);
    famIndex(nI) = i;
    if i <= 7
        plot(t3(nI),1:sum(nI),'.');
    elseif i <= 14
        plot(t3(nI),1:sum(nI),'p');
    else
        plot(t3(nI),1:sum(nI),'v');
    end
    legStr = [legStr; ...
        strcat("Template ",num2str(IgoodUniq(i)),"(",num2str(sum(nI)),")")];

    figure(10);
    hold on;
    if i <= 7
        plot(lonPhase2(nI),latPhase2(nI),'.','linewidth',1); zoom on;
    elseif i <= 14
        plot(lonPhase2(nI),latPhase2(nI),'p','linewidth',1); zoom on;
    else
        plot(lonPhase2(nI),latPhase2(nI),'v','linewidth',1); zoom on;
    end

    hold on; grid on; axis equal;
    legStr2 = [legStr2; ...
        strcat("Template ",num2str(IgoodUniq(i)),"(",num2str(sum(nI)),")")];
    axis(boundaryBox);
end

figure(9);
legend(legStr,'Location','NorthWest');
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

figure(10);
hold on; plot(lonEC,latEC,'k-','linewidth',2);
legend(legStr2,'Location','Best');
geoshow('~/igdata/fallasactualizado/fallas2008completas.shp','Color','r','linewidth',0.5);
grid on; ax = gca; ax.LineWidth = 1.5;
ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
