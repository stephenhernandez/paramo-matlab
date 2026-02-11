%april6_2020
cd ~/products/rsam/;
clear; close all; clc;

maxY = 200;
%load EC.FENY..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.1;
%load EC.VCH1..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.CEAZ..BHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.ALCE..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.FER1..BHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.FER2..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.AMIL..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Villamil'; floorLevel = 0.1;
%load EC.ISPT..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Isla de La Plata'; floorLevel = 1;
%load EC.PPLP..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Puerto Lopez'; floorLevel = 1;
%load EC.AOTA..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat;  ciudad1 = 'Otavalo'; floorLevel = 0.5;
%load EC.ZALD..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
load EC.BREF..BHZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.PARU..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.IESS..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.SALI..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.ZUMB..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.ELAR..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.PARU..HHZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.24MA..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.CIVI..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.AEPN..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.AV03..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.AV05..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.AV11..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.AV18..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.CRPG..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.TING..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat; ciudad1 = 'Quito'; floorLevel = 0.5; % very noisy
%load EC.ALAT..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Quito'; floorLevel = 0.5;
%load EC.ACUE..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'Cuenca'; floorLevel = 0.1;

%load EC.AAT1..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat;  ciudad1 = 'Atuntaqui';

%load EC.APJL..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat;  ciudad1 = 'Pujili'; floorLevel = 0.1;
%load EC.ZALD..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat;
%load EC.BTER..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat;
%load EC.BONI..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat;
%load EC.PULU..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat; ciudad1 = 'GGP'; floorLevel = 0.1;

T = load('EC.AOTA..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat','S'); ciudad2 = 'Otavalo';
%T = load('EC.GYKA..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Guayaquil';
%T = load('EC.AMNT..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S');
%T = load('EC.ABAB..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Babahoyo';
%T = load('EC.AOTA..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Otavalo';
%T = load('EC.AMIL..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Unknown';
%T = load('EC.BREF..BHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Cotopaxi'; floorLevel = 1;
%T = load('EC.ZALD..HNZ_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat','S'); ciudad2 = 'Quito'; floorLevel = 0.5;
%T = load('EC.APO2..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Portoviejo';
%T = load('EC.ARIO..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Riobamba';
%T = load('EC.ALAT..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Latacunga';
%T = load('EC.ALJ1..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Loja';
%T = load('EC.FENY..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat','S'); ciudad2 = 'Quito'; floorLevel = 0.5;
%T = load('EC.PULU..HHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat','S');
S(2,1) = T.S;

%%
zeroPhaseFlag = false;
Sfilt = medfiltWaveforms(S,60*1 + 1,false,zeroPhaseFlag); %90 minutes
%Sfilt = medfiltWaveforms(Sfilt,60*2 + 1,false,true);
Sfilt = resampleWaveforms(Sfilt,4/3600);
%Sfilt = medfiltWaveforms(Sfilt,9,false,true); %90 minutes

%%
t = getTimeVec(Sfilt);
t = t - hours(5);
t2 = getTimeVec(Sfilt(2));
t2 = t2 - hours(5);

%%
alphaLevel = 0.8;
alphaLevel2 = 0.4;
otherDays = (datetime(2019,12,01):dn2dt(ceil(now)))';

%%
dd1 = Sfilt(1).d;
dd2 = Sfilt(2).d;
restrictionMeasureStart = datetime(2020,03,17);

%%
% tI = t >= datetime(2020,01,24) & t <= restrictionMeasureStart;
% tI2 = t2 >= datetime(2020,01,24) & t2 <= restrictionMeasureStart;
maxmax(1) = 1; %median(dd1(tI));
maxmax(2) = 1; %median(dd2(tI2));

%%
tI = t >= otherDays(1);
dd1 = dd1(tI);
t1 = t(tI);

tI2 = t1 <= restrictionMeasureStart;
lt = sum(tI2);
oneWeek = 7*sum(t1 >= otherDays(1) & t1 < otherDays(2));
nWeeks = floor(lt/oneWeek);

dOrig = dd1;
tPred = t1;
lIncludeLockdown = length(tPred);

dd1 = dd1(1:nWeeks*oneWeek);
t1 = t1(1:nWeeks*oneWeek);
nWeeksLockDown = floor(lIncludeLockdown/oneWeek);
extra = mod(lIncludeLockdown,oneWeek);

dd1 = reshape(dd1,[oneWeek,nWeeks]);
badI = dd1 <= floorLevel;
dd1(badI) = NaN;
typicalWeek = nanmedian(dd1,2);
dispersion = 1.4826*mad(dd1,1,2); %try to trasnform to gaussian model
ll = typicalWeek - 1.96*dispersion;
ul = typicalWeek + 1.96*dispersion;

%%
DayNumber = weekday(otherDays);
saturdayI = DayNumber == 7;
sundayI = DayNumber == 1;
sundays = otherDays(sundayI);
saturdays = otherDays(saturdayI);

figure('units','normalized','outerposition',[0 0 1 1]);
sp(1) = subplot(211); P1 = plot(t,Sfilt(1).d/maxmax(1),'linewidth',3); zoom on;
hold on;
sp(2) = subplot(212);
P2 = plot(t2,Sfilt(2).d/maxmax(2),'linewidth',3); zoom on; hold on;
linkaxes(sp,'x');

ylim([0 maxY])
xlim([min(otherDays) dn2dt(ceil(now))])

for i = 1:length(saturdays)
    figure(1);
    subplot(211); Asat(i,1) = area([saturdays(i),saturdays(i)+2],maxY*[1,1]);
    Asat(i,1).LineStyle = 'none';
    Asat(i,1).FaceAlpha = alphaLevel2;
    Asat(i,1).FaceColor = Asat(1,1).FaceColor;
    
    subplot(212);
    Asat(i,2) = area([saturdays(i),saturdays(i)+2],maxY*[1,1]);
    Asat(i,2).LineStyle = 'none';
    Asat(i,2).FaceAlpha = alphaLevel2;
    Asat(i,2).FaceColor = Asat(1,2).FaceColor;
end
weekendFaceColor = Asat(1,1).FaceColor;

%%
% guayGapStarts = [datetime(2020,02,11,20,43,00); datetime(2020,02,26,21,37,00); datetime(2020,03,24,15,28,00)];
% guayGapEnds = [datetime(2020,02,12,17,51,00); datetime(2020,02,27,19,20,00); datetime(2020,03,25,13,25,00)];
% for i = 1:length(guayGapStarts)
% figure(1);
% subplot(212);
% Agap(i) = area([guayGapStarts(i),guayGapEnds(i)],maxY*[1,1]);
% Agap(i).LineStyle = 'none';
% Agap(i).FaceColor = [0.5 0.5 0.5];
% Agap(i).FaceAlpha = 0.4;
% end

%%
subplot(211)
title({['Nivel de Ruido Sismico en ',ciudad1,', Ecuador'],strcat('Estacion: ',Sfilt(1).kstnm,'; Banda de frecuencia analizada: 1 - 8 Hz.')})
ax = gca;
%ax.XTick = datetime(2020,02,01):7:dn2dt(ceil(now));
%ax.XTickLabel = datestr(ax.XTick,'mmm dd');
%ax.XTickLabelRotation = +20;
plot(ax,dn2dt(datenum(2020,03,17).*[1 1]),[0 maxY],'k-','linewidth',3);
ylim([0 maxY])
ylabel('Amplitud relativa al maximo');

subplot(212)
title({['Nivel de Ruido Sismico en ',ciudad2,', Ecuador'],strcat('Estacion: ',Sfilt(2).kstnm,'; Banda de frecuencia analizada: 1 - 8 Hz.')})
ax = gca;
%ax.XTick = datetime(2020,02,01):7:dn2dt(ceil(now));
%ax.XTickLabel = datestr(ax.XTick,'mmm dd');
%ax.XTickLabelRotation = +20;
plot(ax,dn2dt(datenum(2020,03,17).*[1 1]),[0 maxY],'k-','linewidth',3);
ylim([0 maxY])
%ylabel('Amplitud relativa al maximo');

%legend('Quito','Guayaquil')
%title('Nivel de Ruido Sismico [1 - 8 Hz.]')
P1.Color(4) = alphaLevel;
P2.Color(4) = alphaLevel;

%%
figure('units','normalized','outerposition',[0 0 1 0.75]);
p1 = plot(dateshift(t1(1:oneWeek),'start','minute'),dd1,'-','color',[0.5 0.5 0.5],'linewidth',0.1); zoom on;
hold on;
for i = 1:length(p1)
    p1(i).Color(4) = 0.2;
end
p3 = plot(dateshift(t1(1:oneWeek),'start','minute'),ll,'r-','linewidth',4); zoom on;
p4 = plot(dateshift(t1(1:oneWeek),'start','minute'),ul,'r-','linewidth',4); zoom on;
p2 = plot(dateshift(t1(1:oneWeek),'start','minute'),typicalWeek,'k-','linewidth',4); zoom on;

p2.Color(4) = alphaLevel2;
p3.Color(4) = alphaLevel2;
p4.Color(4) = alphaLevel2;
grid on;

ylim([0 maxY]);

%%
predictedTimeSeries = repmat(typicalWeek,[nWeeksLockDown,1]);
predictedTimeSeries = [predictedTimeSeries; typicalWeek(1:extra)];
extendedDispersion = repmat(dispersion,[nWeeksLockDown,1]);
extendedDispersion = [extendedDispersion; extendedDispersion(1:extra)];

%%
figure(1);
subplot(211);
hold on;
pp = plot(tPred,predictedTimeSeries,'k','linewidth',1);
pp.Color(4) = alphaLevel2;

%%
% figure();
% pnew = plot(tPred,dOrig,'linewidth',5);
% pnew.Color(4) = alphaLevel;
% hold on;
% pp = plot(tPred,predictedTimeSeries,'k','linewidth',1);
% pp.Color(4) = alphaLevel;
% zoom on;

%%
% zScore = (dOrig - predictedTimeSeries)./extendedDispersion;
% figure();
% ss_ = scatter(tPred,dOrig,[],zScore,'filled'); zoom on;
% ylim([0 maxY])
% colorbar;
% caxis([-10 10]);
% hold on;
% pPred = plot(tPred,predictedTimeSeries,'k','linewidth',0.5);
% pPred.Color(4) = alphaLevel2;

%%
% figure(); plot(tPred,zScore,'.'); zoom on; ylim([-20 20]);

%%
% figure(); ss_ = scatter(tPred,dOrig,[],(dOrig./predictedTimeSeries)-1,'filled'); zoom on;
% ylim([0 maxY])
% colorbar;
% caxis([0 2]-1);
% hold on;
% pPred = plot(tPred,predictedTimeSeries,'k','linewidth',2);
% pPred.Color(4) = alphaLevel2;

%%
figure('units','normalized','outerposition',[0 0.1 1 0.6]); zoom on;
hold on;
for i = 1:length(saturdays)
    %figure(3);
    %hold on;
    Asat(i,1) = area([saturdays(i),saturdays(i)+2],maxY*[1,1],'FaceColor',weekendFaceColor);
    Asat(i,1).LineStyle = 'none';
    Asat(i,1).FaceAlpha = alphaLevel2;
    Asat(i,1).FaceColor = Asat(1,1).FaceColor;
end

n = length(tPred);
prcntChange = (dOrig./predictedTimeSeries)-1;
p = plot(tPred,dOrig,'linewidth',3);
drawnow;
pEdge = p.Edge;
cdata = NaN(n,3);

baseCMapN = 255;
baseCMap = parula(baseCMapN);
fakeX = (1:baseCMapN)';
caxisMin = -1; caxisMax = 1;
cRange = caxisMax - caxisMin;
transformedX = baseCMapN*(prcntChange + 1)/cRange;

tic;
for i = 1:3
    cdata_ = interp1(fakeX,baseCMap(:,i),transformedX);
    cI = transformedX >= baseCMapN;
    cdata_(cI) = baseCMap(end,i);
    cI = transformedX < 1;
    cdata_(cI) = baseCMap(1,1);
    cdata(:,i) = cdata_;
end
cdata = [uint8(baseCMapN*cdata) uint8(round(alphaLevel*baseCMapN)*ones(n,1))]';
toc;

set(pEdge, 'ColorBinding','interpolated', 'ColorData',cdata);
c = colorbar;
caxis([caxisMin caxisMax]);

% figure(3);
hold on; 
pp = plot(tPred,predictedTimeSeries,'k');
pp(1).Color(4) = alphaLevel2;

plot(dn2dt(datenum(restrictionMeasureStart).*[1 1]),[0 maxY],'k-','linewidth',3);
ylim([0 maxY])
xlim([min(otherDays) dn2dt(ceil(now))])
ax = gca; ax.Box = 'on';
title({['Nivel de Ruido Sismico en ',ciudad1,', Ecuador'],...
    strcat('Estacion: ',Sfilt(1).kstnm,'; Banda de frecuencia analizada: 1 - 8 Hz.')})
ylabel('Amplitud [Cuentas]');
c.Label.Interpreter = 'latex';
c.Label.String = 'Cambio Relativo [fraccion]';

%%
% figure('units','normalized','outerposition',[0 0.1 1 0.7]);
% plot(tPred,prcntChange,'linewidth',1);
% zoom on;

%%
% figure(); ss2 = scatter(tPred,(dOrig./predictedTimeSeries)-1,[],zScore,'filled'); zoom on; ylim([-1 1]);
% ss2.MarkerFaceAlpha = 0.25;
% colorbar;
% caxis([-2 2]);
% ss2.MarkerEdgeColor = 'k';
% ss2.MarkerEdgeAlpha = 0.1;
