clear; close all; clc;

%load ~/research/now/reventador/ReventadorSubspaceDetectorResults_v10.mat
load ~/igdata/ReventadorSubspaceDetectorResults_v10.mat

[stla,stlo] = metaDataFromStationList(["CASC";"BONI";"ANTS";"ANTG"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse)*1e-3;

CASC = 1e9*median(z2p(:,1:3),2,"omitnan")/3.141950e+08;

bI = tabs >= datetime(2020,01,338);
BONI = 1e9*median(z2p(:,4:6),2,"omitnan");
BONI(bI) = BONI(bI)/5.037350e+08;
BONI(~bI) = BONI(~bI)/2.014940e+09;

bI = tabs >= datetime(2014,01,226);
ANTS = 1e9*median(z2p(:,7:9),2,"omitnan");
ANTS(bI) = ANTS(bI)/5.037350e+08;
ANTS(~bI) = ANTS(~bI)/3.141950e+08;

bI = tabs >= datetime(2017,01,105);
ANTG = 1e9*median(z2p(:,10:12),2,"omitnan");
ANTG(bI) = ANTG(bI)/3.141950e+08;
ANTG(~bI) = ANTG(~bI)/3.141950e+08;

amps = [CASC BONI ANTS ANTG];
fourI = sum(isfinite(amps),2) == 4;

fourI = fourI & CASC >= 100;
fourI = fourI & BONI >= 40;
fourI = fourI & ANTS >= 80;
fourI = fourI & ANTG >= 40;

G_ = full(Gvdcc(4));
G_ = G_(1:end-1,:);
G_ = [getDD(log10(d_)) getDD(d_) G_];

findI = find(fourI);
findI = sort(findI(randsample(length(findI),5e3)));

%
mbest_ = [];
rdum = logspace(log10(5),log10(700),501)';

G = [];
d = [];
for i = 1:length(findI)
    amps_ = [CASC(findI(i)) BONI(findI(i)) ANTS(findI(i)) ANTG(findI(i))];
    G = [G; G_];
    d = [d; -getDD(log10(amps_'))];
    disp(findI(i))
end

G = [G; 0 0 1 1 1 1];
d = [d; 0];

mbest = pinv(G)*d;
mbest_ = [mbest_ mbest];
format long g
disp(mbest);

att1 = mbest_(1)*log10(rdum) + mbest_(2)*rdum;  %hernandez
att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer

figure('units','normalized','outerposition',[0 0 1 1]);
semilogx(rdum,att1,'.'); zoom on; grid on; hold on;
semilogx(rdum,att2,'.');
ylim([0 6]);

gamma17 = -interp1(rdum,att1,17) + 2;
gamma100 = -interp1(rdum,att1,100) + 3;

attGamma17 = att1 + gamma17;
hold on;
semilogx(rdum,attGamma17,'.','color',[0.5 0.5 0.5]); grid on;

attGamma100 = att1 + gamma100;
hold on;
semilogx(rdum,attGamma100,'.','color',[0.3 0.3 0.3]); grid on;

legend('new','uhrhammer','fix17','fix100','location','northwest');

%%
McorrOrig = [CASC BONI ANTS ANTG];
Mcorr = McorrOrig; %pre-allocate
for i = 1:4
    Mcorr(:,i) = -4.3+log10(Mcorr(:,i)) + (mbest(1)*log10(d_(i)) + mbest(2)*d_(i) + gamma100 + mbest(i+2));
end
fprintf('mean mad: %g\n',mean(mad(Mcorr,1,2),'omitnan'));

figure();
M1 = median(Mcorr,2,"omitnan");
t = tabs;
mI = M1 >= 0; %-0.6;
plot(t(mI),M1(mI),'.'); zoom on;

nDays = 14;
[rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = t2r(t(mI),days(nDays),M1(mI));
figure(); plot(t(mI),rate/nDays,'.'); zoom on; grid on;
figure(); plot(t(mI),meanMagsFixedTimeWin,'.'); zoom on; grid on;
figure(); semilogy(t(mI),sumEnergyFixedTimeWin/nDays,'.'); zoom on; grid on;
figure(); plot(t(mI),medianMagsFixedTimeWin,'.'); zoom on; grid on;
figure(); semilogy(t(mI),(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); zoom on; grid on;


%%
clear; close all; clc;
T2 = load('~/igdata/ReveMagnitudes_v2.mat');
tabs = T2.t;
t = tabs;
magErr = T2.magErr;
M1 = T2.M1; %median(Mcorr,2,"omitnan") - 2; %medMag - 2;
mI = M1 >= 0.2 & magErr <= 0.5& t >= datetime(2021,01,01);

nDays = 7;
[rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = t2r(t(mI),days(nDays),M1(mI));
%figure(); plot(t(mI),rate/nDays,'.'); zoom on; grid on;
%figure(); plot(t(mI),meanMagsFixedTimeWin,'.'); zoom on; grid on;
%figure(); semilogy(t(mI),sumEnergyFixedTimeWin/nDays,'.'); zoom on; grid on;
%figure(); plot(t(mI),medianMagsFixedTimeWin,'.'); zoom on; grid on;
%figure(); semilogy(t(mI),(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); zoom on; grid on;

%cd ~/research/now/reventador/
cd ~/igdata
T = readtable("data2.xlsx");
texcel = T.Fecha;
% M_x_DeT_VENTO
% M_x_DeFN1
% M_x_DeFN2
% M_x_DeFN3
% M_x_DeFN4
% M_x_DeFlancoSurEste
% M_x_DeFlancoNorEste
% M_x_DeFlancoEste
% M_x_DeSOBRE
% M_x_DeBAJO
% MAXIMAGEN
% DailyAlertsTotal
% DailyAlertsFlancoNorte
% MinDespejadosPorD_a
% MinDespejado_07_17_
% x_d_aNoObservado

vento = T.M_x_DeT_VENTO;
fn1 = T.M_x_DeFN1;
fn2 = T.M_x_DeFN2;
fn3 = T.M_x_DeFN3;
fn4 = T.M_x_DeFN4;
Tsobre = T.M_x_DeSOBRE;
Tbajo = T.M_x_DeBAJO;
minDespejados = T.MinDespejadosPorD_a;

uniqueDays = unique(dateshift(texcel,'start','day'));
lDays = length(uniqueDays);
maxvento = NaN(lDays,1); medvento = NaN(lDays,1);
maxfn1 = NaN(lDays,1); medfn1 = NaN(lDays,1);
maxfn2 = NaN(lDays,1); medfn2 = NaN(lDays,1);
maxfn3 = NaN(lDays,1); medfn3 = NaN(lDays,1);
maxfn4 = NaN(lDays,1); medfn4 = NaN(lDays,1);
maxTsobre = NaN(lDays,1);
maxTbajo = NaN(lDays,1);
%
nPerDay = NaN(lDays,1);
for i = 1:lDays
    day_ = uniqueDays(i);
    tI = texcel >= day_ & texcel < day_+1;
    nPerDay(i) = sum(nPerDay);
    vento_ = vento(tI);
    fn1_ = fn1(tI);
    fn2_ = fn2(tI);
    fn3_ = fn3(tI);
    fn4_ = fn4(tI);
    Tsobre_ = Tsobre(tI);
    Tbajo_ = Tbajo(tI);

    bajoI = Tbajo_ > 0;
    if sum(bajoI)
        maxTbajo(i) = max(Tbajo_(bajoI));
    end

    sobreI = Tsobre_ > 0;
    if sum(sobreI)
        maxTsobre(i) = max(Tsobre_(sobreI));
    end

    fn1I = fn1_ > 0;
    if sum(fn1I)
        maxfn1(i) = max(fn1_(fn1I));
        medfn1(i) = median(fn1_(fn1I),"omitnan");
    end

    fn1I = fn1_ > 0;
    if sum(fn1I)
        maxfn1(i) = max(fn1_(fn1I));
        medfn1(i) = median(fn1_(fn1I),"omitnan");
    end
    fn2I = fn2_ > 0;
    if sum(fn2I)
        maxfn2(i) = max(fn2_(fn2I));
        medfn2(i) = median(fn2_(fn2I),"omitnan");
    end
    fn3I = fn3_ > 0;
    if sum(fn3I)
        maxfn3(i) = max(fn3_(fn3I));
        medfn3(i) = median(fn3_(fn3I),"omitnan");
    end
    fn4I = fn4_ > 0;
    if sum(fn4I)
        maxfn4(i) = max(fn4_(fn4I));
        medfn4(i) = median(fn4_(fn4I),"omitnan");
    end
    ventoI = vento_ > 0;
    if sum(ventoI)
        maxvento(i) = max(vento_(ventoI));
        medvento(i) = median(vento_(ventoI),"omitnan");
    end
end

dailyEnergy14Days = (10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays;
dailyEnergy14DaysInterp = interp1(datenum(t(mI)),dailyEnergy14Days,datenum(uniqueDays));
rateInterp = interp1(datenum(t(mI)),rate/nDays,datenum(uniqueDays));
medianInterp = interp1(datenum(t(mI)),medianMagsFixedTimeWin,datenum(uniqueDays));

zerophaseFlag = true;

% figure 1
figure('units','normalized','outerposition',[0 0 2/4 1]);
ax(1) = subplot(311);
semilogy(uniqueDays,rateInterp,'-','linewidth',4); zoom on; grid on;
title("Daily Explosion Rate");

ax(2) = subplot(312); plot(uniqueDays,medfiltSH(maxvento,1,zerophaseFlag),'.'); zoom on; grid on;
hold on; plot(uniqueDays,medfiltSH(maxvento,nDays,zerophaseFlag),'-','linewidth',4);
title("Max Lateral Vent Temperature");

ax(3) = subplot(313);
plot(uniqueDays,[medfiltSH(maxfn1,1,true) medfiltSH(maxfn2,1,true) medfiltSH(maxfn3,1,true) medfiltSH(maxfn4,1,true)],'.','color',[0.5 0.5 0.5]);
hold on;
ax(3).ColorOrderIndex = 1;
hh = plot(uniqueDays,[medfiltSH(maxfn1,nDays,true) medfiltSH(maxfn2,nDays,true) medfiltSH(maxfn3,nDays,true) medfiltSH(maxfn4,nDays,true)],'-','linewidth',5);
grid on;
title("Max North Flank Temperatures");
linkaxes(ax,'x');
legend(hh,"Pulse 1","Pulse 2","Pulse 3","Pulse 4","Location","Best");

% figure 2
figure('units','normalized','outerposition',[0 0 2/4 1]);
ax(1) = subplot(311);
pp = scatter(t(mI),M1(mI),10,[0.5 0.5 0.5],'filled'); hold on; pp.MarkerFaceAlpha = 0.5;
ax(1).ColorOrderIndex = 1;
plot(uniqueDays,medianInterp,'-','linewidth',4); zoom on; grid on;
title("Median Pseudo-Magnitude, 14-Day Moving Window");


ax(2) = subplot(312); 
vI = vento > 0;
plot(texcel(vI),vento(vI),'.'); zoom on; grid on; hold on;
pp1 = plot(texcel(vI),medfiltSH(vento(vI),nDays,zerophaseFlag),'-','linewidth',4);
pp1.Color(4) = 0.8;
pp2 = plot(texcel(vI),movmax(vento(vI),3),'-','linewidth',3);
pp2.Color(4) = 0.5;
title("Max Lateral Vent Temperatures");
legend([pp1 pp2],strcat(num2str(nDays),'-Day Moving Median'),'3-Day Moving Maximum',...
    'location','best');

ax(3) = subplot(313);
plot(uniqueDays,[medfiltSH(maxfn1,1,zerophaseFlag) medfiltSH(maxfn2,1,zerophaseFlag) ...
    medfiltSH(maxfn3,1,zerophaseFlag) medfiltSH(maxfn4,1,zerophaseFlag)],'.','color',[0.5 0.5 0.5]);
hold on;
ax(3).ColorOrderIndex = 1;
hh = plot(uniqueDays,[medfiltSH(maxfn1,nDays,zerophaseFlag) medfiltSH(maxfn2,nDays,zerophaseFlag) ...
    medfiltSH(maxfn3,nDays,zerophaseFlag) medfiltSH(maxfn4,nDays,zerophaseFlag)],'-','linewidth',5);
grid on;
title("Max North Flank Temperatures");
linkaxes(ax,'x');
legend(hh,"Pulse 1","Pulse 2","Pulse 3","Pulse 4","Location","Best");
%axis tight;
xlim([datetime(2021,05,01) datetime(2023,01,10)]);

% figure 3
figure('units','normalized','outerposition',[0 0 2/4 1]);
ax(1) = subplot(311);
semilogy(t(mI),(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); zoom on; grid on;
title("$\propto$ Energy");

ax(2) = subplot(312);
plot(uniqueDays,medfiltSH(maxTsobre,1,zerophaseFlag),'.','color',[0.5 0.5 0.5]);
hold on;
ax(2).ColorOrderIndex = 1;
plot(uniqueDays,medfiltSH(maxTsobre,nDays,zerophaseFlag),'-','linewidth',5);
grid on;
title("Max T Above Summit");

ax(3) = subplot(313);
plot(uniqueDays,medfiltSH(maxTbajo,1,zerophaseFlag),'.','color',[0.5 0.5 0.5]);
hold on;
ax(3).ColorOrderIndex = 1;
plot(uniqueDays,medfiltSH(maxTbajo,nDays,zerophaseFlag),'-','linewidth',5);
grid on;
title("Max T Below Summit");
linkaxes(ax,'x');
xlim([datetime(2021,05,01) datetime(2023,01,10)]);

% figure 4
figure('units','normalized','outerposition',[0 0 2/4 1]);
ax(1) = subplot(311);
semilogy(uniqueDays,rateInterp,'-','LineWidth',4); zoom on; grid on;
title("Daily Explosion Rate");

ax(2) = subplot(312);
plot(uniqueDays,medfiltSH(maxTsobre,1,zerophaseFlag),'.','color',[0.5 0.5 0.5]);
hold on;
ax(2).ColorOrderIndex = 1;
plot(uniqueDays,medfiltSH(maxTsobre,nDays,zerophaseFlag),'-','linewidth',5);
grid on;
title("Max T Above Summit");

ax(3) = subplot(313);
plot(uniqueDays,medfiltSH(maxTbajo,1,zerophaseFlag),'.','color',[0.5 0.5 0.5]);
hold on;
ax(3).ColorOrderIndex = 1;
plot(uniqueDays,medfiltSH(maxTbajo,nDays,zerophaseFlag),'-','linewidth',5);
grid on;
title("Max T Below Summit");
linkaxes(ax,'x');
axis tight;
%xlim([datetime(2021,05,01) datetime(2023,01,10)]);

% figure 5
figure('units','normalized','outerposition',[0 0 3/4 1]);
clear ax;
ax(1) = subplot(211);
vI = vento > 0;
plot(texcel(vI),vento(vI),'.'); zoom on; grid on; hold on;
pp1 = plot(texcel(vI),medfiltSH(vento(vI),nDays,zerophaseFlag),'-','linewidth',4);
pp1.Color(4) = 0.8;
pp2 = plot(texcel(vI),movmax(vento(vI),3),'-','linewidth',3);
pp2.Color(4) = 0.5;
legend([pp1 pp2],strcat(num2str(nDays),'-Day Moving Median'),'3-Day Moving Maximum');
title("Max Lateral Vent Temperatures");
ax(2) = subplot(212);
mI = minDespejados > 0;
bar(texcel(mI)+0.5,minDespejados(mI)/60,1); zoom on; grid on; hold on;
linkaxes(ax,'x');
ylim(ax(2),[0 24]);
emptyI = minDespejados < 0;
t2 = texcel(emptyI);
title(ax(2),"Visibility");
ylabel(ax(2),"Hours / Day");
for i = 1:length(t2)
    ax(2).ColorOrderIndex = 2;
    A = area(ax(2),[t2(i) t2(i)+1],[24 24]);
    A.FaceAlpha = 0.25; A.EdgeColor = "none";
end
linkaxes(ax,'x');
axis tight;

% figure 6
figure('units','normalized','outerposition',[0 0 3/4 1]);
ax(1) = subplot(211);
semilogy(uniqueDays,rateInterp,'-','linewidth',4); zoom on; grid on;
title(ax(1),"Daily Explosion Rate");
ax(2) = subplot(212); vI = vento > 0;
plot(texcel(vI),vento(vI),'.'); zoom on; grid on; hold on;
pp = plot(texcel(vI),movmax(vento(vI),3),'-','linewidth',3);
title(ax(2),"Max Lateral Vent Temperatures");
linkaxes(ax,'x');
pp.Color(4) = 0.8;
legend(pp,'3-Day Moving Maximum');
