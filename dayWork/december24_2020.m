% loop over dayfiles of REVN data applying Roman's harmonic tremor detector
% successively

clear; close all;

%% preallocate
maxN = 10000;
detection_times = NaT(maxN,1);
fund_freq = NaN(maxN,1);
harm_freq_1 = fund_freq;
harm_freq_2 = fund_freq;
harm_freq_3 = fund_freq;
HSI_1 = fund_freq;
HSI_2 = fund_freq;
HSI_3 = fund_freq;
fund_pow = fund_freq;

%%
waveform_dir = "~/rawdata_cotopaxi"; %fullfile("~","Desktop");
% dayStart = datetime(2025,09,12);
% dayEnd = datetime(2025,11,07);
%dayStart = datetime(2020,12,16);
dayStart = datetime(2025,10,14);
dayEnd = dayStart;

dayInc = 1;
dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);

nOvertones = 4;
lfc = 0.4;
hfc = 16;
MINMAX = 10; %10;
MINMED = 3;
MINMIN = 2;
MINFREQ = 0.3;
MINPOW = 1e-18;
MINPSD = 1e-6;
MAXINTERHARM = 3e-3;
plotFlag = true;

knetwk = "EC";
kstnm = "REVN";
khole = "";
kcmpnm = "HHZ";

n = 1;
for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,1,kstnm,kcmpnm,knetwk,khole,true,true,waveform_dir);
    if isnat(S.ref) | isempty(S)
        fprintf("skipping: %s\n",day_);
        continue;
    end

    fprintf("processing: %s\n",day_);
    delta = S.delta;
    Fs = 1/delta;
    Sf = detrendWaveforms(S);
    Sf = taperWaveforms(Sf,Fs*10);
    Sf = filterWaveforms(Sf,lfc,hfc);

    if Fs < 100
        Sf = resampleWaveforms(Sf,100);
    end
    Sf = nanGapWaveforms(Sf,0);

    %%
    [detection_times_,overtoneFreq_,harmStrength_,fundPow_,...
        t,dcut,overtonePow,interharmPow,overtoneI,powerOrig,pxx,fxx,fundPowOrig] = ...
        tremometer_control(Sf,nOvertones,MINMAX,MINMED,MINMIN,MINFREQ,MINPOW,MINPSD,plotFlag);

    %harmStrengthOrig = overtonePow(:,1:nOvertones)./interharmPow(:,1:nOvertones);

    %%
    ldetections = length(detection_times_);
    if ~ldetections
        fprintf("processed day: %s, but no detections found. moving on... \n\n",day_)
        continue;
    end

    detection_times(n:n+ldetections-1) = detection_times_;
    fund_freq(n:n+ldetections-1) = overtoneFreq_(:,1);
    harm_freq_1(n:n+ldetections-1) = overtoneFreq_(:,2);
    harm_freq_2(n:n+ldetections-1) = overtoneFreq_(:,3);
    harm_freq_3(n:n+ldetections-1) = overtoneFreq_(:,4);
    HSI_1(n:n+ldetections-1) = harmStrength_(:,1);
    HSI_2(n:n+ldetections-1) = harmStrength_(:,2);
    HSI_3(n:n+ldetections-1) = harmStrength_(:,3);
    fund_pow(n:n+ldetections-1) = fundPow_;
    n = n + ldetections;
end

%%
if n < 2
    return;
end

detection_times(n:end) = [];
fund_freq(n:end) = [];
harm_freq_1(n:end) = [];
harm_freq_2(n:end) = [];
harm_freq_3(n:end) = [];
HSI_1(n:end) = [];
HSI_2(n:end) = [];
HSI_3(n:end) = [];
fund_pow(n:end) = [];