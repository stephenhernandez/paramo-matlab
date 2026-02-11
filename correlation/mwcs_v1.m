clear;
close all;

%%
pinoFlag = 0;
if pinoFlag
    cd ~/research/now/pichincha_nxcorr
    load sample_PINO_data_for_Hugo.mat
    tStart = 6;
    tEnd = 20;
    Fs = 100;
    plotFlag = false;
    maxDV = 5; % percent
    lfc = 2;
    hfc = 8;
else
    cd ~/research/now/dV/reventadorInfrasound/
    load azu_temps.txt
    t_temperature = datetime(azu_temps(:,1),azu_temps(:,2),azu_temps(:,3),azu_temps(:,4),azu_temps(:,5),azu_temps(:,6));
    temperature = azu_temps(:,7);
    
    %%
    if ~exist('symStack','var')
        %load AZU_microphone01_7sec11sec_singleComponent_1Hz2Hz.mat
        load AZU_microphone01_11sec15sec_singleComponent_1Hz2Hz.mat
        indiv_events = acaus;
    end
    tStart = 13;
    tEnd = 17;
    Fs = 32;
    plotFlag = false;
    maxDV = 5; % percent
end

%%
goodThreshCC = 0.7;
upsampleRate = 1; %not using
winDur = 1; %window length in seconds
overlapDur = 0.5;
detrendFlag = true;
indiv_events = detrend(indiv_events);
indiv_events = indiv_events./rssq(indiv_events);
if pinoFlag
    refTrace = pws(indiv_events(:,1:29)); %detrend(nanmean(indiv_events(:,1:29),2));
else
    refTrace = pws(indiv_events); %detrend(nanmean(indiv_events(:,1:29),2));
end
refTrace = refTrace./rssq(refTrace);
[dcutRef,startIndex,endIndex] = cutWindows(refTrace,winDur*Fs,overlapDur*Fs,detrendFlag);

%%
twinstart = tdum(startIndex);
twinend = tdum(endIndex);
goodI = twinstart >= tStart & twinend < tEnd;

%%
twinstart = twinstart(goodI);
dcutRef = resample(dcutRef(:,goodI),upsampleRate,1);
dcutRef = dcutRef./rssq(dcutRef);
[~,lags] = doCrossCorrFreqDom(dcutRef,dcutRef);
zeroIndex = find(lags == 0);

%%
nTest = size(indiv_events,2);
dt_t_slope = NaN(nTest,1);
for i = 1:nTest
    disp(i);
    indiv_events_ = indiv_events(:,i);
    dcutTest_ = cutWindows(indiv_events_,winDur*Fs,overlapDur*Fs,detrendFlag);
    dcutTest_ = resample(dcutTest_(:,goodI),upsampleRate,1);
    %dcutTest_ = normalizeWaveforms(dcutTest_);
    dcutTest_ = dcutTest_./rssq(dcutTest_);
    %[maxccp,plags,maxccn,nlags] = doccFreqCircShift(dcutTest_);
    data1 = doCrossCorrFreqDom(dcutRef,dcutTest_);
    
    [peakCC,maxI] = max(data1);
    maxI = maxI';
    %dt_t_slope_ = robustfit(twinstart(peakCC >= 0.5),maxI(peakCC >= 0.5) - 256);
    meetThresh = peakCC >= goodThreshCC;
    if sum(meetThresh)>3
        %dt_t_slope_ = polyfit(twinstart(meetThresh),maxI(meetThresh)-zeroIndex,1);
        %dt_t_slope(i) = dt_t_slope_(1);
        
        dt_t_slope_ = robustfit(twinstart(meetThresh),(maxI(meetThresh)-zeroIndex)/(Fs*upsampleRate));
        dt_t_slope(i) = 100*dt_t_slope_(2);
    end
    %figure(); ss = scatter(twinstart,maxI,[],peakCC,'filled'); zoom on; c = colorbar;
end

%%
figure();
plot(dt_t_slope,'o');
zoom on; grid on;
%zeroIndex = find(lags == 0);

%%
close all;
tic;
[prcntVector,peakCC,cc,eps,err_rms] = ...
    cwiShifts(refTrace,indiv_events,tStart,tEnd,Fs,plotFlag,maxDV,lfc,hfc);

figure();
subplot(211);
plot(prcntVector,'.'); zoom on; hold on; plot(prcntVector+err_rms,'k-'); plot(prcntVector-err_rms,'k-');

subplot(212);
plot(t,prcntVector,'.'); zoom on;
%hold on; plot(t,err_rms,'k-'); plot(t,-err_rms,'k-');
if ~pinoFlag
    hold on;
    
    [prcntVector2,peakCC,cc,eps,err_rms] = ...
        cwiShifts(refTrace,indiv_events,8,10,Fs,plotFlag,maxDV,lfc,hfc);
    
    plot(t,prcntVector2,'.'); zoom on;
    
    yyaxis right;
    pp = plot(t_temperature,temperature,'linewidth',2); pp.Color(4) = 0.8;
end
toc;
legend('AZU 13 - 17 sec.','AZU 7 - 11 sec.','temperature');