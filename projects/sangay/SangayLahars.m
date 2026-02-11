%function SangayLahars
clear; close all; clc;

cd ~/products/rsam/temporary/
%load EC.SAGA..HHZ_10Hz_30DUR_MedAmpRMS_preFiltFalse.mat; thresh = 50;
%load EC.SAGA..HHZ_10Hz_60DUR_MeanAbsAmp_preFiltFalse.mat; thresh = 50;
%load EC.SAGA..HHZ_4Sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat; thresh = 5;
%load EC.SAGA..HHZ_10Hz_60DUR_MedRMSAmp_preFiltFalse.mat; thresh = 100;
load EC.SAGA..HHZ_10Hz_60DUR_MeanAbsAmp_preFiltFalse.mat; thresh = 100; 
%load EC.SAGA..HHZ_20Hz_60DUR_MedRMSAmp_preFiltFalse.mat; thresh = 100;

%
%tStart = datetime(2019,12,01);
%Scut = cutWaveforms(S,tStart,0,datetime(2023,03,10,13,30,00)-tStart);
tStart = datetime(2013,09,01);
Scut = cutWaveforms(S,tStart,0,datetime(2023,10,01)-tStart);
S = Scut;

spm = 1;
lpNminutes = 30;
npoles = 1;
zeroPhaseFlag = false;
% %exponent = 2; ampThresh = 1e5; %1e2;
% exponent = 1; ampThresh = 1e3; %1e2;
% if zeroPhaseFlag
%     rateThresh = 1e3*0.25/lpNminutes;
% else
%     rateThresh = 1e3*0.50/lpNminutes;
% end
% minPeakDistance = lpNminutes;

lpNminutes = 30;
npoles = 1;
zeroPhaseFlag = 0;
exponent = 1;
ampThresh = 800;
rateThresh = 0.15;
minPeakDistance = lpNminutes;
dispMax = 140; %200;

%
Sf = S;
Sf = nanGapWaveforms(Sf,0);
t = getTimeVec(Sf);
d = Sf.d;
tI = t <= datetime(2018,11,11);
d(tI) = d(tI)*4;
dI = d >= thresh;
d(~dI) = NaN;
Sf1 = dealHeader(Sf,d);
Sf = Sf1; %nanGapWaveforms(Sf1,0); %  interpolateWaveforms(Sf1);
d = Sf.d;
d = medfiltSH(d,3,true);
d = d.^exponent;
if zeroPhaseFlag
    d = zpkFilter(d,-inf,1/lpNminutes,spm,npoles,zeroPhaseFlag);
else
    d = fftfilt(ones(lpNminutes,1)/lpNminutes,d);
end
Sf = dealHeader(Sf,d);

%
Sf2 = Sf;
Sf2 = dealHeader(Sf2,Sf2.d);
t2 = getTimeVec(Sf2);

S_prefilt = Sf2;
Sf2 = filterWaveforms(Sf2,-inf,1/(lpNminutes*spm*60),npoles,[],zeroPhaseFlag,false);
Sf = Sf2;
rate = differentiateWaveforms(Sf2);
rateTSOrig_ = rate.d;
rate = dealHeader(rate,rateTSOrig_);

%
close all;
Sorig = [Sf2; rate];

df3 = pull(rate);
amp = pull(Sf2);
[ratePeak,rateI] = findpeaks(df3,'MINPEAKHEIGHT',rateThresh); %,'MINPEAKDISTANCE',minPeakDistance);
[ampPeak,ampI] = findpeaks(amp,'MINPEAKHEIGHT',ampThresh); %,'MINPEAKDISTANCE',minPeakDistance);
ampI2 = rateI;

%%
for i = 1:length(rateI)-1
    r1 = rateI(i);
    r2 = rateI(i+1);
    aI = ampI>r1 & ampI < r2;
    aI = find(aI);
    ampI2_ = ampI(aI);
    amp_ = amp(ampI2_);
    if isempty(ampI2_)
        ampI2(i) = 0;
    elseif length(ampI2_) == 1
        ampI2(i) = ampI2_;
    else
        aI_ = find(amp_ == max(amp_));
        ampI2(i) = ampI2_(aI_(1));
    end
end
rateI = rateI(1:end-1);
ratePeak = ratePeak(1:end-1);
ampI2(end) = [];
badI = ampI2==0;
ampI2(badI) = [];
rateI(badI) = [];
ratePeak(badI) = [];

%%
zci = @(v,x) find(diff(sign(v)) == x);
riseI = zci(df3,2);
peakI = zci(df3,-2);

clear peakRI maxR peakI2;
nn = 0;
for i = 1:length(riseI)-1
    df3_ = df3(riseI(i):riseI(i+1));
    [max_,maxI_] = max(df3_); 
    maxI2 = peakI(i+1) - riseI(i);
    if max_ >= rateThresh
        nn = nn + 1;                        
        maxR(nn,1) = max_;                  
        peakRI(nn,1) = riseI(i) + maxI_;    
        peakI2(nn,1) = riseI(i)+maxI2;      
    end
end

%%
t1 = getTimeVec(Sf2);
t = getTimeVec(rate);
peakI = peakI2;

tPeakFlow = t1(ampI2); %peakI+1);
peakFlow = amp(ampI2); %interp1(datenum(t1),Sf2.d,datenum(tPeakFlow));
tPeakRate = t(rateI); %peakRI-1);
peakRate = ratePeak; %maxR;
FlowAmpAtPeakRate = amp(rateI); %interp1(datenum(t1),Sf2.d,datenum(tPeakRate));

%%
goodI = peakRate >= rateThresh & peakFlow >= ampThresh;
peakRate(~goodI) = [];
tPeakRate(~goodI) = [];
peakFlow(~goodI) = [];
tPeakFlow(~goodI) = [];
FlowAmpAtPeakRate(~goodI) = [];

%
close all; 
figure(); 
MarkerSize = 18;
ax(1) = subplot(211); 
pp = plot(t1,Sf2.d,'linewidth',2); zoom on; grid on; hold on; pp.Color(4) = 0.8;
plot(tPeakFlow,peakFlow,'.','markersize',MarkerSize); 
plot(tPeakRate,FlowAmpAtPeakRate,'p','markersize',MarkerSize,'linewidth',2); 

ax(2) = subplot(212); 
pp = plot(t,rate.d,'linewidth',2); zoom on; grid on; hold on; pp.Color(4) = 0.8;
plot(tPeakRate,peakRate,'.','markersize',MarkerSize); 
plot([min(tPeakRate) max(tPeakRate)],[0 0],'--','Color',[0.5 0.5 0.5],...
    'linewidth',2);
linkaxes(ax,'x')
ax(1).YScale = 'log';

figure(); 
loglog(peakRate,peakFlow,'.'); zoom on; grid on;
xlabel("Peak Rate");
ylabel("Peak Flow");

figure(); 
plot(tPeakFlow,cumsum(peakFlow),'.'); zoom on; grid on;
ylabel("Cumulative Flow$^{2}$");

figure(); 
loglog(peakRate,minutes(tPeakFlow - tPeakRate),'.'); zoom on; grid on;
xlabel("Peak Rate");
ylabel("Minutes From Peak Rate Time to Peak Flow Time");

figure(); 
loglog(peakFlow,minutes(tPeakFlow - tPeakRate),'.'); zoom on; grid on;
xlabel("Peak Flow");
ylabel("Minutes From Peak Rate Time to Peak Flow Time");

iet = minutes(diff(tPeakFlow));
figure(); 
loglog(sort((iet)),1 - (0:length(iet)-1)'/length(iet),'.'); zoom on;
grid on; 
xlabel("Inter-event Waiting Time");

figure(); 
loglog(sort(peakFlow),1 - (0:length(peakFlow)-1)'/length(iet),'.'); zoom on;
grid on; 
xlabel("Peak RSAM");

figure(); 
loglog(FlowAmpAtPeakRate,peakFlow./FlowAmpAtPeakRate,'.'); zoom on; grid on;
xlabel("Flow Amp At Peak Rate");
ylabel("Peak Flow / Flow Amp At Peak Rate");

[~,sI] = sort(peakFlow,'descend');

TTtmp = timetable(tPeakRate(sI(1:dispMax)),peakRate(sI(1:dispMax)),...
    peakFlow(sI(1:dispMax)),minutes(tPeakFlow(sI(1:dispMax)) - tPeakRate(sI(1:dispMax))),...
    'VariableNames',{'peak_rate','peak_amplitude','min. after peak rate'});
disp(TTtmp);

%%
t = getTimeVec(rate);
tlahar = t(rateI);
peakRate = ratePeak;

[peakRate,pI] = sort(peakRate,'descend');
tlahar = tlahar(pI);

pI = tlahar >= datetime(2013,09,01);
peakRate = peakRate(pI);
tlahar = tlahar(pI);

%%
maxDisp = length(peakRate); %1000;
maxDisp = min([length(tlahar) maxDisp]);

df2 = pull(Sf2);
vq = interp1(datenum(t2),df2,datenum(tlahar(1:maxDisp)));

[tsort,sortI] = sort(tlahar(1:maxDisp));
peakSort = peakRate(1:maxDisp);
peakSort = peakSort(sortI);
vqSort = vq(sortI);

%%
[peakLaharAmp,locsLahar] = findpeaks(df2,'MINPEAKHEIGHT',ampThresh,'MINPEAKDISTANCE',minPeakDistance);

%%
tIndices = [zeros(size(tlahar)); ones(size(locsLahar))];
ampMerge = [peakRate; peakLaharAmp];
tMerge = [tlahar; t2(locsLahar)];

[tMerge,sI] = sort(tMerge);
tIndices = tIndices(sI);
ampMerge = ampMerge(sI);

tIndicesDiff = diff(tIndices);
peakRateLocID = find(tIndicesDiff == 1);
peakLaharAmpLocID = peakRateLocID+1;

%%
t_peakRate = tMerge(peakRateLocID);
t_peakAmp = tMerge(peakLaharAmpLocID);
peakRate = ampMerge(peakRateLocID);
peakAmp = ampMerge(peakLaharAmpLocID);
minAfter = minutes(t_peakAmp - t_peakRate);

[peakAmp,sI] = sort(peakAmp,'descend');
t_peakRate = t_peakRate(sI);
t_peakAmp = t_peakAmp(sI);
peakRate = peakRate(sI);
minAfter = minAfter(sI);

sI = minAfter <= 8*60 & peakAmp >= ampThresh ...
    & ~(t_peakAmp >= datetime(2019,08,12) & t_peakAmp < datetime(2019,08,13)) ...
    & ~(t_peakAmp >= datetime(2019,03,20) & t_peakAmp < datetime(2019,03,21)) ...
    & ~(t_peakAmp >= datetime(2019,01,18) & t_peakAmp < datetime(2019,01,19)) ...
    & ~(t_peakAmp >= datetime(2021,08,12) & t_peakAmp < datetime(2021,08,13)) ...
    & ~(t_peakAmp >= datetime(2014,08,03) & t_peakAmp < datetime(2014,08,04)) ...
    & ~(t_peakAmp >= datetime(2016,08,01) & t_peakAmp < datetime(2016,08,02)) ...
    & ~(t_peakAmp >= datetime(2020,12,20) & t_peakAmp < datetime(2020,12,21)) ...
    & ~(t_peakAmp >= datetime(2021,03,31) & t_peakAmp < datetime(2021,04,01)) ...
    & ~(t_peakAmp >= datetime(2021,07,16) & t_peakAmp < datetime(2021,07,17)) ...
    & ~(t_peakAmp >= datetime(2018,09,07) & t_peakAmp < datetime(2018,09,08)) ...
    & ~(t_peakAmp >= datetime(2018,01,31) & t_peakAmp < datetime(2018,02,01));
t_peakRate = t_peakRate(sI);
t_peakAmp = t_peakAmp(sI);
peakRate = peakRate(sI);
minAfter = minAfter(sI);
peakAmp = peakAmp(sI);

TTtmp = timetable(t_peakRate,peakRate,peakAmp,minAfter,...
    'VariableNames',{'peak_rate','peak_amplitude','min. after peak rate'});
disp(TTtmp);

%%
figure('units','normalized','outerposition',[0 0 0.7 1]);
subplot(211);
loglog(peakRate,minAfter,...
    '.','markersize',20); zoom on; grid on;

subplot(212);
loglog(peakRate,peakAmp,...
    '.','markersize',20); zoom on; grid on;

[t_peakRate,sI] = sort(t_peakRate);
t_peakAmp = t_peakAmp(sI);
peakRate = peakRate(sI);
minAfter = minAfter(sI);
peakAmp = peakAmp(sI);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
axDist(1) = subplot(121);
histogram(peakAmp,...
    logspace(log10(floor(min(peakAmp))),log10(ceil(max(peakAmp))),101));
hold on;
plot([3000 3000],[0 60],'--','Color',[0.5 0.5 0.5],'linewidth',2);
ax = gca; ax.XScale = 'log'; zoom on; grid on;

axDist(2) = subplot(122);
semilogx(sort(peakAmp),length(peakAmp) - (1:length(peakAmp))','.');
grid on; zoom on; hold on;
plot([3000 3000],[1 1200],'--','Color',[0.5 0.5 0.5],'linewidth',2);
ylabel('Number $\geq$ Peak Amp.');
linkaxes(axDist,'x');

%%
% figure(1);
% hold(axWave(2),'on');
% plot(axWave(2),t_peakRate,peakRate,'p'); zoom on;
% plot(axWave(2),[min(t) max(t)],[rateThresh rateThresh],'k--','linewidth',1); zoom on;
% grid(axWave(2),'on');
% 
% figure(1);
% hold(axWave(1),'on');
% plot(axWave(1),t_peakAmp,peakAmp,'.','markersize',25); zoom on;
% plot(axWave(1),[min(t) max(t)],[ampThresh ampThresh],'k--','linewidth',1); zoom on;
% plot(axWave(1),[min(t) max(t)],[1e3 1e3],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on;
% plot(axWave(1),[min(t) max(t)],[2e3 2e3],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on;
% plot(axWave(1),[min(t) max(t)],[3e3 3e3],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on;
% plot(axWave(1),[min(t) max(t)],[4e3 4e3],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on;
% plot(axWave(1),[min(t) max(t)],[5e3 5e3],'--','linewidth',1,'Color',[0.5 0.5 0.5]); zoom on;
% 
% tdumb = getTimeVec(Sorig(1));
% vq = interp1(datenum(tdumb),Sorig(1).d,datenum(t_peakRate));
% plot(axWave(1),t_peakRate,vq,'p'); zoom on;
% 
% axWave(1).YScale = 'log';

%%
figure('units','normalized','outerposition',[0 0 1 1]);
% ax = subplot(221);
% plot(t_peakRate,cumsum(peakRate),'.'); zoom on; grid on;

ax(1) = subplot(2,2,[1 2]);
plot(t_peakAmp,cumsum(peakAmp.^exponent),'.'); zoom on; grid on;

[N,edges] = histcounts(t_peakAmp,...
    dateshift(min(t_peakAmp),'start','day'):days(1):dateshift(max(t_peakAmp),'end','day')+1);

N = N';
edges = edges(1:end-1)';
ax(2) = subplot(2,2,[3 4]);
stairs(edges,N); zoom on;
ylabel('Number');

zoom on; grid on;
% linkaxes(ax,'x');

% figure(4);
% ax = subplot(221);
% hold on;
% plot(t_peakRate(peakAmp>=3000),cumsum(peakRate(peakAmp>=3000)),'.'); zoom on; grid on;
% legend('all','$\geq$3000','location','northwest');

ax(1) = subplot(2,2,[1 2]);
hold on;
plot(t_peakAmp(peakAmp>=3000),cumsum(peakAmp(peakAmp>=3000).^exponent),'.'); zoom on; grid on;
legend('all','$\geq$3000','location','northwest');
ylabel('$\propto$ Cum. Energy');

[N,edges] = histcounts(t_peakAmp(peakAmp>=3000),...
    dateshift(min(t_peakAmp),'start','day'):days(1):dateshift(max(t_peakAmp),'end','day')+1);

N = N';
edges = edges(1:end-1)';
ax(2) = subplot(2,2,[3 4]);
hold on;
stairs(edges,N); zoom on;
ylabel('Number / Day');
legend('all','$\geq$3000','location','northwest');

zoom on; grid on;
linkaxes(ax,'x');

nDays = 90; r = t2r(t_peakAmp(peakAmp>=3000),days(nDays),[],false)/nDays;
figure(); plot(t_peakAmp(peakAmp>=3000),r,'.'); zoom on; grid on;

%%
tlaharShift = hours(t_peakAmp - dateshift(t_peakAmp,'start','day'));
figure();
plot(sort(tlaharShift),((1:length(tlaharShift))'/length(tlaharShift))); zoom on;
grid on;
xlim([0 24]);

