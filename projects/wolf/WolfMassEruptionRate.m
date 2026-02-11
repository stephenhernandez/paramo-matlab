clear; close all; 
kstnm = "FER1"; 
kcmpnm = "BHZ";

%
cd ~/research/now/wolf/

load WolfRSAM_TimeSeries.mat
tEnd = S.ref + S.e;
tEnd = dateshift(tEnd,'start','day');
TEXTWRITE = false;

%%
try
    %S = rmsGather(datetime(2022,01,01),datetime(2022,01,29),60,5,7,kstnm,kcmpnm,true,false,"EC",""); % true-false is mean of abs value
    W = rmsGather(tEnd-1,dn2dt(ceil(now + 5/24)),60,5,7,kstnm,kcmpnm,true,false,"EC",""); % true-false is mean of abs value
    S = mergeWaveforms([S; W]);
    clear W;
    save('WolfRSAM_TimeSeries.mat','S');
catch
    fprintf(2,'Could not update FER1 data\n');
    return;
end
ltcorrection = 0;

%
tLast = dn2dt(ceil(datenum(S.ref+S.e)));
fig = figure('units','normalized','outerposition',[0 0 1 1]); 
hold on; 
t = getTimeVec(S);
t = t - hours(ltcorrection);
ll = plot(t,fftfilt(ones(1,1)/1,S.d),'.','linewidth',3);
ll.Color(4) = 0.8; grid on; zoom on;
ax = gca;
ax.YScale = 'log';
title(kstnm);

% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% hold on;
% t = getTimeVec(S);
% t = t - hours(ltcorrection);
% ll = plot(t,fftfilt(ones(10,1)/10,S.d),'.','linewidth',3);
% ll.Color(4) = 0.8; grid on; zoom on;
% ax = gca;
% ax.YScale = 'log';
% title(kstnm);

figure(1);
hold on;
t = getTimeVec(S);
t = t - hours(ltcorrection);
dd = S.d; dd(~isfinite(dd)) = 0;
dsmooth = medfiltSH(dd,1440,true);
%dsmooth = zpkFilter(medfiltSH(dd,3,true),-inf,1/30,1,1,1);
ll = plot(t,dsmooth,'linewidth',8);
ll.Color(4) = 0.8; grid on; zoom on;
ax = gca;
ax.YScale = 'log';
title(kstnm);

if TEXTWRITE
    tI = t >= datetime(2021,12,01);
    t = t(tI);
    dsmooth = dsmooth(tI);
    Sf = dealHeader(S,dsmooth,1/60,t(1));
    sac2txt(Sf,'~/FER1_RSAM_5Hz7Hz_30minuteMedianFilter.txt');

    t = getTimeVec(S);
    d = double(pull(S));
    tI = t >= datetime(2021,12,01);
    t = t(tI);
    d = d(tI);
    Sf = dealHeader(S,d,1/60,t(1));
    sac2txt(Sf,'~/FER1_RSAM_5Hz7Hz_NoSmoothing.txt');
end

figure(1);
xlim([datetime(2021,12,01) tLast])

%%
t = getTimeVec(S);
d = double(pull(S));

%%
close all;
d = d.^1;
dsmooth = medfiltSH(d,30,true);
badI = ~isfinite(dsmooth) | d <= 0.2;

smoothingWindow = 60*24;
eruptionStartTime = datetime(2022,01,07,05,20,00);
noiseWinStart = datetime(2022,01,04,06,00,00);
noiseWinEnd = datetime(2022,01,07,02,20,00);

dsmooth(badI) = 0;
dsmooth2 = zpkFilter(dsmooth,-inf,1/smoothingWindow,1,2,false);

dsmooth2(badI) = NaN;

%%
tI = t >= noiseWinStart;
noiseCorection = mean(dsmooth2(tI & t <= noiseWinEnd),'omitnan');

cumulativeEruptedLava = (1/60)*1e6*112; %in m^3
lowerBound = 0.5*cumulativeEruptedLava;
upperBound = 1.5*cumulativeEruptedLava;

%%

figure();

postEruptionI = tI & t >= (eruptionStartTime + minutes(smoothingWindow)) & ~badI;
tPostEruption = t(postEruptionI);
estimationBand =  tPostEruption <= datetime(2022,01,22); % hard code

baseEruptionRateTrend = dsmooth2(postEruptionI) - noiseCorection;
baseEruptionRateTrend = baseEruptionRateTrend / sum(baseEruptionRateTrend(estimationBand),'omitnan');

plot(tPostEruption,cumulativeEruptedLava*baseEruptionRateTrend,'.'); zoom on;
grid on;
ax = gca; ax.YScale = 'log';
xlabel('UTC'); 
ylabel('[$m^{3} \cdot s^{-1}$]');

figure(1); hold on;
ll = plot(tPostEruption,lowerBound*baseEruptionRateTrend,'k-'); zoom on;
ll.Color(4) = 0.5;

figure(1); hold on;
ll = plot(tPostEruption,upperBound*baseEruptionRateTrend,'k-'); zoom on;
ll.Color(4) = 0.5;

title({'Volcan Wolf Mass Eruption Rate'; ['Smoothing Window: ',num2str(smoothingWindow),' minutes']});


