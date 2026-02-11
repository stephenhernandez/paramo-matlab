function mwcs(refTrace,testTraces,tStart,tEnd,secdur,overlapdur,goodThresh,Fs) %,lfc,hfc)
%tStart = 8;
%tEnd = 20;
%Fs = 100;

%%
%goodThresh = 0.6;
%upsampleRate = 1;
%secdur = 2;
%overlapdur = 1.99;
detrendFlag = true;
testTraces = normalizeWaveforms(detrend(testTraces));
%refTrace = nanmean(testTraces(:,1:29),2);
refTrace = normalizeWaveforms(detrend(refTrace));
[dcutRef,startIndex,endIndex] = cutWindows(refTrace,secdur*Fs,overlapdur*Fs,detrendFlag);

%%
twinstart = tdum(startIndex);
twinend = tdum(endIndex);
goodI = twinstart >= tStart & twinend < tEnd;

%%
twinstart = twinstart(goodI);
dcutRef = resample(dcutRef(:,goodI),upsampleRate,1);
dcutRef = normalizeWaveforms(dcutRef);
[~,lags] = doCrossCorrFreqDom(dcutRef,dcutRef);
zeroIndex = find(lags == 0);

%%
nTest = size(testTraces,2);
dt_t_slope = NaN(nTest,1);
for i = 1:nTest
    disp(i);
    indiv_events_ = testTraces(:,i);
    dcutTest_ = cutWindows(indiv_events_,secdur*Fs,overlapdur*Fs,detrendFlag);
    dcutTest_ = resample(dcutTest_(:,goodI),upsampleRate,1);
    dcutTest_ = normalizeWaveforms(dcutTest_);
    
    %[maxccp,plags,maxccn,nlags] = doccFreqCircShift(dcutTest_);
    data1 = doCrossCorrFreqDom(dcutRef,dcutTest_);
    
    [peakCC,maxI] = max(data1);
    maxI = maxI';
    %dt_t_slope_ = robustfit(twinstart(peakCC >= 0.5),maxI(peakCC >= 0.5) - 256);
    meetThresh = peakCC >= goodThresh;
    if sum(meetThresh)>3
        %dt_t_slope_ = polyfit(twinstart(meetThresh),maxI(meetThresh)-zeroIndex,1);
        %dt_t_slope(i) = dt_t_slope_(1);
        
        dt_t_slope_ = robustfit(twinstart(meetThresh),maxI(meetThresh)-zeroIndex);
        dt_t_slope(i) = dt_t_slope_(2);
    end
    %figure(); ss = scatter(twinstart,maxI,[],peakCC,'filled'); zoom on; c = colorbar;
end
figure(); plot(dt_t_slope,'o'); zoom on; grid on;

