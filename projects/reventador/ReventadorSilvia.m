clear; close all; clc;
%T = readtable("pulsos_all.xlsx");
cd ~/research/now/reventador/
T = readtable("data2.xlsx");
%texcel = T.FECHA;
%t = datetime(texcel,'ConvertFrom','excel');
t = T.Fecha;
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

vento = T.M_x_DeT_VENTO;%T_VENTO;
fn1 = T.M_x_DeFN1; %FN1;
fn2 = T.M_x_DeFN2; %FN2;
fn3 = T.M_x_DeFN3; %FN3;
fn4 = T.M_x_DeFN4; %FN4;

uniqueDays = unique(dateshift(t,'start','day'));
lDays = length(uniqueDays);
maxvento = NaN(lDays,1); medvento = NaN(lDays,1);
maxfn1 = NaN(lDays,1); medfn1 = NaN(lDays,1);
maxfn2 = NaN(lDays,1); medfn2 = NaN(lDays,1);
maxfn3 = NaN(lDays,1); medfn3 = NaN(lDays,1);
maxfn4 = NaN(lDays,1); medfn4 = NaN(lDays,1);
%
nPerDay = NaN(lDays,1);
for i = 1:lDays
    day_ = uniqueDays(i);
    tI = t >= day_ & t < day_+1;
    nPerDay(i) = sum(nPerDay);
    vento_ = vento(tI);
    fn1_ = fn1(tI);
    fn2_ = fn2(tI);
    fn3_ = fn3(tI);
    fn4_ = fn4(tI);

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
figure(); 
plot(uniqueDays,[maxfn1 maxfn2 maxfn3 maxfn4],'o'); zoom on; grid on; 
%hold on; plot(uniqueDays,[medfn1 medfn2 medfn3 medfn4],'.');

%
T2 = load('~/igdata/ReveMagnitudes_v2.mat');
tabs = T2.t;
t = tabs;
magErr = T2.magErr;
medMag = T2.M1; %median(Mcorr,2,"omitnan") - 2; %medMag - 2;
goodI = medMag >= 0.2 & magErr <= 0.5& t >= datetime(2021,01,01);
tGood = tabs(goodI);
aGood = medMag(goodI); %10.^(3+1.44*medMag(goodI));

% pulso 1: 17/05/2021 - 30/10/2021
% pulso 2: 11/10/2021 - 05/07/2022
% pulso 3: 21/06/2022 - 30/12/2022
% pulso 4: 07/11/2022 - 09/01/2023
tStart = [datetime(2021,05,17);...
    datetime(2021,10,11);...
    datetime(2022,06,21);...
    datetime(2022,11,07)];

tEnd = [datetime(2021,10,30);...
    datetime(2022,07,05);...
    datetime(2022,12,30);...
    datetime(2023,01,09)];

%
nDays = 3;
[rate,meanamp] = t2r(tGood,days(nDays),aGood);
figure();
axx(1) = subplot(311);
semilogy(tGood,rate/nDays,'.'); grid on;

axx(2) = subplot(312);
semilogy(tGood,meanamp,'.'); grid on;

axx(3) = subplot(313);
semilogy(tGood,rate.*meanamp/nDays,'.'); zoom on; grid on;

linkaxes(axx,'x');

%
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
durations = NaT(lstarts,1) - NaT(1);
totEnergy = medianAmps;
for i = 1:lstarts
    ia = tGood >= tStart(i) & tGood <= tEnd(i);
    axF(i) = subplot(lstarts,1,i);
    thisT = tGood(ia);
    thisA = aGood(ia);
    semilogy(thisT,thisA,'o');
    %meanAmp_ = mean(thisA);
    medianAmp_ = median(thisA,"omitnan");

    medianAmps(i) = medianAmp_;
    durations(i) = tEnd(i) - tStart(i);
    averageRate(i) = sum(ia)/days(durations(i));
    totEnergy(i) = sum(thisA);
    title(strcat('Median Amplitude: ',num2str(medianAmp_),...
        ', Duration: ',num2str(days(durations(i))),...
        ', Average Rate: ',num2str(averageRate(i))));
    axis(axF(i),'tight');
end
linkaxes(axF,'y');
