clear; close all; clc;

Fs = 128;
secDur = 64;
lowestFreq = 0.05;
highestFreq = 50;
nfreqs = 1+1024*4;
freqVec = logspace(log10(lowestFreq),log10(highestFreq),nfreqs)';


%tStart = datetime(2018,07,01);
%tEnd = datetime(2018,07,31);
tStart = datetime(2012,01,01);
tEnd = datetime(2012,01,10);

dayVec = (tStart:tEnd)';
ldays = length(dayVec);

hvMain = NaN(nfreqs,ldays);
winlen = secDur.*Fs;
nfft = secDur*Fs;
w = kaiser(nfft/2,30);
interpMethod = "pchip";

kstnm = ["BMAS";"BULB"];
kcmpnms = ["BHZ"];
%kstnm = "CASC";
%kcmpnms = ["BHZ";"BHN";"BHE"];
%kcmpnms = ["HHZ";"HHN";"HHE"];
for i = 1:ldays
    day_ = dayVec(i);
    tic;
    S = loadWaveforms(day_,1,kstnm,kcmpnms,"EC","");
    if any(isnat(pull(S,'ref')))
        fprintf('no data on day: %s\n',datestr(day_));
        continue;
    end
    toc;

    S = [S; S(end)];
    S = differentiateWaveforms(S);
    S = resampleWaveforms(detrendWaveforms(S),Fs);
    S = syncWaveforms(detrendWaveforms(nanGapWaveforms(detrendWaveforms(S),0)));

    nOverlap = 0;
    for j = 1:3
        d_ = double(pull(S(j)));
        if j == 1
            dcut1 = doAutoCorrFreqDom(cutWindows(d_,winlen,nOverlap,true));
        elseif j == 2
            dcut2 = doAutoCorrFreqDom(cutWindows(d_,winlen,nOverlap,true));
        else
            dcut3 = doAutoCorrFreqDom(cutWindows(d_,winlen,nOverlap,true));
        end
    end
    toc;

    normers = rssq([dcut1; dcut2; dcut3])'; %sqrt(median(abs([dcut1; dcut2; dcut3])))';
    badWindowsI = normers <= 1 | peak2rms(dcut1)' >= median(peak2rms(dcut1),'omitnan');
    normers(badWindowsI) = [];
    dcut1(:,badWindowsI) = [];
    dcut2(:,badWindowsI) = [];
    dcut3(:,badWindowsI) = [];

    dcut1 = dcut1 ./ normers';
    dcut2 = dcut2 ./ normers';
    dcut3 = dcut3 ./ normers';

    %     [pxx1,~] = pwelch(dcut1,w,0,freqVec,Fs);
    %     [pxx2,~] = pwelch(dcut2,w,0,freqVec,Fs);
    %     [pxx3,~] = pwelch(dcut3,w,0,freqVec,Fs);

    [pxx1,fxx] = pwelch(dcut1,w,[],nfft,Fs);
    [pxx2,~] = pwelch(dcut2,w,[],nfft,Fs);
    [pxx3,~] = pwelch(dcut3,w,[],nfft,Fs);

    %
    pxx1 = sqrt(pxx1);
    pxx2 = sqrt(pxx2);
    pxx3 = sqrt(pxx3);

    %
    pxx1 = interp1(fxx,pxx1,freqVec,interpMethod,"extrap");
    pxx2 = interp1(fxx,pxx2,freqVec,interpMethod,"extrap");
    pxx3 = interp1(fxx,pxx3,freqVec,interpMethod,"extrap");

    %
    hv = sqrt((mean(pxx2,2,'omitnan')+mean(pxx3,2,'omitnan'))./mean(pxx1,2,'omitnan')); %ratio of mean spectra
    hvMain(:,i) = hv;
    fprintf('done with: <strong>%s</strong>\n',datestr(day_));
    toc;
end

%%
close all;
[~,maxI] = max(hvMain);
figure(); plot(dayVec,freqVec(maxI),'.'); zoom on;
figure(); loglog(freqVec,hvMain,'.');
zoom on; grid on; hold on;
loglog(freqVec,mean(hvMain,2,'omitnan'),'k-','linewidth',5); zoom on; grid on;
title(kstnm);
xlabel("frequency");
ylabel("amplification");
xlim([lowestFreq highestFreq]);
