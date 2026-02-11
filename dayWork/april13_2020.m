%cd ~/data/2020/EC/PUYO/HHZ.D/
clear; close all; clc;
tic;
if isunix
    rawDataDir = '~/rawdata/';
else
    rawDataDir = '~/data/';
end

tStart = datetime(2012,01,01);
tEnd = dn2dt(floor(now));

%tStart = datetime(2015,01,260);
%tEnd = datetime(2015,01,279);

%tStart = datetime(2020,04,06);
%tEnd = datetime(2020,04,13);

diffFlag = false;
[~,tabs,pks,templateIndex,maxAmpRMS,pksOrig] = ...
    singleStationRepeaterSearch(0.5,'~/research/now/sangay/puyo_template_2020.mat',4,1,60,tStart,tEnd,1e6,diffFlag,[],[],rawDataDir);
toc;

if isunix
    cd ~/research/now/sangay
    save('doubleTemplatePUYOAnalysisAtSangay');
end

%%
clear; close all;
cd ~/research/now/sangay/
load doubleTemplatePUYOAnalysisAtSangay.mat


%%
cd ~/research/now/sangay/
clear; close all; clc; 
load ~/research/now/sangay/SangayRegionalAnalysis_v2.mat

%%
% lastT = max(tabs);
% tNew = load('~/research/now/sangay/sangay_elevenSensors_update','tabs'); tNew = tNew.tabs;
% mar = load('~/research/now/sangay/sangay_elevenSensors_update','maxAmpRMS'); mar = mar.maxAmpRMS;
% tINew = load('~/research/now/sangay/sangay_elevenSensors_update','templateIndex'); tINew = tINew.templateIndex;
% stds = load('~/research/now/sangay/sangay_elevenSensors_update','saveStds'); stds = stds.saveStds;
% pksONew = load('~/research/now/sangay/sangay_elevenSensors_update','pksOrig'); pksONew = pksONew.pksOrig;
% pksNew = load('~/research/now/sangay/sangay_elevenSensors_update','pks'); pksNew = pksNew.pks;
% neffNew = load('~/research/now/sangay/sangay_elevenSensors_update','Neff'); neffNew = neffNew.Neff;
% 
% tGood = tNew > lastT;
% maxAmpRMS = [maxAmpRMS; mar(tGood,:)];
% tabs = [tabs; tNew(tGood)];
% templateIndex = [templateIndex; tINew(tGood)];
% 
% pksOrig = [pksOrig; pksONew(tGood)];
% pks = [pks; pksNew(tGood)];
% Neff = [Neff; neffNew(tGood)];
% 
% clear tNew mar tINew stds pksONew pksNew neffNew lastT;
% save('~/research/now/sangay/SangayRegionalAnalysis_v1');

%
% % cd ~/research/now/sangay/;
% % templateFileName = '~/research/now/sangay/threeSangayTemplatesElevenSensorsPWS';
% % loopStart = dateshift(max(tabs),'start','day');
% % loopEnd = dateshift(dn2dt(now)+hours(5),'start','day');
% % [indiv_events,tabs,pks,templateIndex,maxAmpRMS,pksOrig,saveStds,Neff] = singleStationRepeaterSearch(7,templateFileName,3,1,120,loopStart,loopEnd,5e7,false,false);
% % save('sangay_elevenSensors_update');

[tabs,sI] = sort(tabs);
maxAmpRMS = maxAmpRMS(sI,:);
Neff = Neff(sI);
pks = pks(sI);
templateIndex = templateIndex(sI);
pksOrig = pksOrig(sI);

sI = tabs >= datetime(2000,01,01);
tabs = tabs(sI);
maxAmpRMS = maxAmpRMS(sI,:);
Neff = Neff(sI);
pks = pks(sI);
templateIndex = templateIndex(sI);
pksOrig = pksOrig(sI);

[tabs,pksOrig,removeIndices] = removeRepeatedMatches(tabs,pksOrig,60);
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    maxAmpRMS(rI,:) = [];
    Neff(rI) = [];
    templateIndex(rI) = [];
    pks(rI) = [];
end

close all;
figure('units','normalized','outerposition',[0 0 1 1]); plot(sort(tabs),1:length(tabs),'o'); zoom on; grid on;
figure('units','normalized','outerposition',[0 0 1 1]); hold on; for i = 1:length(unique(templateIndex)); plot(tabs(templateIndex == i),1:length(tabs(templateIndex == i)),'.'); end
zoom on; grid on;
figure('units','normalized','outerposition',[0 0 1 1]); hold on; for i = 1:length(unique(templateIndex)); plot(tabs(templateIndex == i),pksOrig(templateIndex == i),'.'); end
zoom on; grid on;
figure('units','normalized','outerposition',[0 0 1 1]); hold on; for i = 1:length(unique(templateIndex)); pp(i) = semilogy(tabs(templateIndex == i),maxAmpRMS(templateIndex == i,9),'.'); end
zoom on; grid on;
ax = gca;
ax.YScale = 'log';
pp(2).Marker = '.';

[tSort,sI] = sort(tabs);
maxAmpRMS = maxAmpRMS(sI,:);
pks = pks(sI);
templateIndex = templateIndex(sI);
pksOrig = pksOrig(sI);

diffSeconds = seconds(diff(tSort));
secThresh = 60;
badI = false(size(sI));
badI2 = find(diffSeconds <= secThresh);
badI3 = badI2+1;
badPKS1 = pksOrig(badI2);
badPKS2 = pksOrig(badI3);
secondGood = badPKS2 >= badPKS1;
badI2(secondGood) = badI3(secondGood);

badI(badI2) = true;

dumbIndex = (1:length(diffSeconds))';
figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(dumbIndex,seconds(diff(tSort)),'.'); zoom on; grid on;
hold on;
plot(dumbIndex(badI2),diffSeconds(badI2),'.');

badI = badI | Neff == 0 | (Neff == 1 & pksOrig <= 13) | ...
    (Neff == 2 & pksOrig <= 12.5) | ...
    (Neff == 3 & pksOrig <= 12) | ...
    (Neff == 4 & pksOrig <= 11.5) | ...
    (Neff == 5 & pksOrig <= 11) | ...
    (Neff == 6 & pksOrig <= 10.5) | ...
    (Neff == 7 & pksOrig <= 10) | ...
    (Neff == 8 & pksOrig <= 9.5) | ...
    (Neff == 9 & pksOrig <= 9) | ...
    (Neff == 10 & pksOrig <= 8.5) |...
    (Neff >= 11 & pksOrig <= 8);

figure('units','normalized','outerposition',[0 0 1 1]);
ss_ = scatter(tSort(~badI),1:length(tSort(~badI)),[],Neff(~badI),'filled'); zoom on; grid on; colorbar;
ss_.MarkerFaceAlpha = 0.4;
ss_.MarkerEdgeColor = 'k';
ss_.MarkerEdgeAlpha = 0.1;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for i = 1:length(unique(templateIndex))
    plot(tabs(templateIndex == i & ~badI),1:length(tabs(templateIndex == i & ~badI)),'.');
end
zoom on; grid on;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
for i = 1:length(unique(templateIndex))
    semilogy(tabs(templateIndex == i & ~badI),maxAmpRMS(templateIndex == i & ~badI,:),'.');
end
zoom on; grid on;
ax = gca;
ax.YScale = 'log';

figure('units','normalized','outerposition',[0 0 1 1]); hold on; for i = 1:length(unique(templateIndex)); plot(tabs(templateIndex == i & ~badI),pksOrig(templateIndex == i & ~badI),'.'); end; zoom on; grid on;

t_ = sort(tabs(~badI));

%
figure('units','normalized','outerposition',[0 0 1 1]);
Nmed = 151;
%plot(t_(1:end-1),1./seconds(diff(t_)),'.'); hold on;
plot(t_(1:end-1),3600./medfiltSH(seconds(diff(t_)),Nmed,true),'linewidth',2); zoom on;

%
sensBBIL = 4.872110e+08;
sensBPAT = sensBBIL;
sensBRTU = 2.748990e+08;
sensBULB = sensBBIL;
sensCHSH = 3.141950e+08;
sensCOHC1 = 2.014940e+09;
sensCOHC2 = 5.037350e+08;
tCOHC1 = datetime(2016,01,076);
sensJSCH1 = sensCHSH;
sensJSCH2 = 1.256780e+09;
tJSCH1 = datetime(2018,01,207);
sensPORT = sensCHSH;
sensPUYO = sensCHSH;
sensTAIS = sensCHSH;
sensTAMH = sensCHSH;

sensitivities = [sensBBIL; sensBPAT; sensBRTU; sensBULB; sensCHSH; sensCOHC1; ...
    sensJSCH1; sensPORT; sensPUYO; sensTAIS; sensTAMH];
kstnm = ["BBIL","BPAT","BRTU","BULB","CHSH","COHC","JSCH","PORT","PUYO","TAIS","TAMH"];
[stla,stlo,stel] = metaDataFromStationList(kstnm);
refEllipse = referenceEllipsoid('wgs84');
r = distance(stla,stlo,-2.0051,-78.3415,refEllipse)*1e-3;
t_ = tabs(~badI);
mags = maxAmpRMS(~badI,:);
for i = 1:length(kstnm)
    kstnm_ = kstnm(i);
    if strcmp(kstnm_,"COHC")
        tI = t_ <= tCOHC1;
        mags(tI,i) = mags(tI,i)/sensCOHC1;
        mags(~tI,i) = mags(~tI,i)/sensCOHC2;
        mags_ = log10(mags(:,i));
        %mags_ = mags_ + 2*log10(r(i)) + 0.00769*r(i) + 3;
        mags_ = mags_ + 1.11*log10(r(i)) + 0.00189*r(i) + 4.5;
        mags(:,i) = mags_;
    elseif strcmp(kstnm_,"JSCH")
        tI = t_ <= tJSCH1;
        mags(tI,i) = mags(tI,i)/sensJSCH1;
        mags(~tI,i) = mags(~tI,i)/sensJSCH2;
        mags_ = log10(mags(:,i));
        %mags_ = mags_ + 2*log10(r(i)) + 0.00769*r(i) + 3;
        mags_ = mags_ + 1.11*log10(r(i)) + 0.00189*r(i) + 4.5;
        mags(:,i) = mags_;
    else
        disp(i);
        mags(:,i) = mags(:,i)/sensitivities(i);
        mags_ = log10(mags(:,i));
        %mags_ = mags_ + 2*log10(r(i)) + 0.00769*r(i) + 3;
        mags_ = mags_ + 1.11*log10(r(i)) + 0.00189*r(i) + 4.5;
        mags(:,i) = mags_;
    end
end

mags(~isfinite(mags)) = NaN;
mags(mags < 0) = NaN;

mest = nanmedian(mags,2);
merr = mad(mags,1,2);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_,mest,'o'); zoom on;
yyaxis right
plot(t_,merr,'.'); zoom on;

log10energy = 1.5*mest + 4.8; %convert to joules
energy = 10.^log10energy;
cumEnergy = cumsum(energy);

figure('units','normalized','outerposition',[0 0 1 1]);
SmoothCumEnergy = cumsum(medfiltSH(energy,Nmed,true));
plot(t_,SmoothCumEnergy,'.'); zoom on;

% hold on;
% plot(t_,cumsum(medfiltSH(energy,201,true)),'o'); zoom on;

%close all; 
figure('units','normalized','outerposition',[0 0 1 1]);
histogram(tabs(~badI),dateshift(min(tabs(~badI)),'start','day'):hours(12):dateshift(max(tabs(~badI)),'end','day')); zoom on;
yyaxis right
plot(tabs(~badI),Neff(~badI),'.'); zoom on;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(211); 
plot(t_,mest,'.'); zoom on; 
subplot(212);
plot(t_,mest,'.'); zoom on; xlim([datetime(2020,06,01) dateshift(dn2dt(now)+hours(5),'end','day')])

%hold on; plot(t_,medfiltSH(mest,201,true),'linewidth',6);

%            2015    2016    2017   2018  2019
totEnergy = [1.2e10 3.1e10 1.5e10 1.4e10 2.5e10];
% need to determine whether SAGA data shows similar trends

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(211)
histogram(tabs(~badI),dateshift(min(tabs(~badI)),'start','day'):hours(24):dateshift(max(tabs(~badI)),'end','day')); 
zoom on;
ylabel('Numero Diario');
xlim([datetime(2019,05,01) datetime(2020,06,14)]);
grid on;
hold on; 
plot([datetime(2020,06,01) datetime(2020,06,01)],[0 250],'k--');
plot([datetime(2020,06,14) dateshift(dn2dt(now)+hours(5),'end','day')],[0 250],'k--');
subplot(212)
histogram(tabs(~badI),dateshift(min(tabs(~badI)),'start','day'):hours(24):dateshift(max(tabs(~badI)),'end','day')); 
zoom on;
xlim([datetime(2020,06,01) dateshift(dn2dt(now)+hours(5),'end','day')]);
ylabel('Numero Diario');
ylim([0 200]);
grid on;
disp(['Last Point: ',datestr(max(t_),30)])

%%
%ulba_template_2020
%cd ~/data/2020/EC/BULB/BHZ.D/
%clear; close all; clc; [~,tabs,pks,templateIndex,maxAmpRMS,pksOrig] = ...
%singleStationRepeaterSearch(0.6,'~/research/now/sangay/ulba_template_2020.mat',1,1,60,datetime(2020,01,01),dn2dt(floor(now)),1e4,0);

%%
plotFlag = false;
if plotFlag
    maxAmpRMSNEW = maxAmpRMS; templateIndexNew = templateIndex;
    pksNew = pks; pksOrigNew = pksOrig;
    tabsNew = tabs;
    load ~/research/now/sangay/puyo_match_filter_analysis.mat
    tabs(end)
    figure(); plot(tabsNew,1:length(tabsNew),'o'); zoom on;
    tI = tabs >= datetime(2020,01,10);
    pksOrig(1:10)
    figure(); plot(tabs,pksOrig,'.'); zoom on;
    figure(); plot(maxAmpRMS,pksOrig,'.'); zoom on;
    clearvars -except tabs pks templateIndex maxAmpRMS pksOrig tabsNew pksNew templateIndexNew maxAmpRMSNEW pksOrigNew
    tI = tabs >= datetime(2020,01,10);
    tabs(tI) = []; maxAmpRMS(tI) = []; pks(tI) = []; templateIndex(tI) = [];
    pksOrig(tI) = [];
    tabs = [tabs; tabsNew]; maxAmpRMS = [maxAmpRMS; maxAmpRMSNEW];
    pks = [pks; pksNew]; templateIndex = [templateIndex; templateIndexNew]; pksOrig = [pksOrig; pksOrigNew];
    figure(); histogram(tabs,dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
    figure(); histogram(tabs,dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day')); zoom on;
    figure(3); hold on; plot(tabs,pksOrig,'o'); zoom on;
    figure(); plot(tabs,maxAmpRMS,'o'); zoom on;
    figure(); semilogy(tabs,cumsum(maxAmpRMS),'.'); zoom on;
    N = 101;
    rollingMed = medfiltSH(maxAmpRMS,N);
    figure(); semilogy(tabs,rollingMed,'.'); zoom on;
    tt = datenum(tabs);
    N = 51;
    iet = tt(N:end) - tt(1:end-N+1); % in days
    figure(); semilogy(tabs(N:end),rollingMed(N:end)./iet,'.'); zoom on;
end
