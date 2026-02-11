%january11_2020
clear; close all; clc;
cd ~/research/now/pichincha_nxcorr/
load('GGP20082020');

%%
snippetLength = 1*5.12;
Fs = 100;
upSampleFactor = 1;
noiseStopIndex = 550*upSampleFactor;
Fs = Fs*upSampleFactor;
nfft = nfft*upSampleFactor;
indiv_events = resample(indiv_events,upSampleFactor,1);
errorInSeconds = 0;
circshiftN = errorInSeconds*Fs;

%%
mI = maxAmp <= 4;
stack = normalizeWaveforms(nanmean(normalizeWaveforms(indiv_events(:,mI)),2));
stack = circshift(stack,-circshiftN);
stack = 1*normalizeWaveforms(detrend(stack(noiseStopIndex:noiseStopIndex+snippetLength*Fs-1)));

%%
indiv_events = detrend(indiv_events(noiseStopIndex:noiseStopIndex+snippetLength*Fs-1,:));

%%
tw = 0.06;
[pxxN,fxxN] = pmtm(stack,{4,'trace'},8*nfft,Fs);
[pxxS,fxxS] = pmtm(indiv_events(:,1057),{4,'trace'},8*nfft,Fs); %12164

figure(); loglog(fxxN,pxxN); hold on; loglog(fxxS,pxxS); zoom on;
fI = fxxN > 0 & fxxN <= 40;
figure(); loglog(fxxS(fI),pxxS(fI)./pxxN(fI)); zoom on;

%%
stf = conv(stack,flipud(indiv_events(:,12164)));
[pxxS,fxxS] = pmtm(stf,{4,'trace'},4*nfft,Fs); %12164


fI = fxxS >= 3 & fxxS <= 40;
figure(); 
loglog(fxxS(fI),pxxS(fI)); zoom on;

load GGPSNR.mat
SNR = SNR.^0.5; nfft = 2^14;

U = fft(indiv_events(:,12164),nfft);
G = fft(stack,nfft);
stfWiener = fftshift(ifft((SNR.*conj(G).*U)./(SNR.*conj(G).*G + 1),[],'symmetric'));
hold on;
[pxxS,fxxS] = pmtm(stfWiener,{4,'trace'},nfft,Fs); %12164
fI = fxxS >= 3 & fxxS <= 40;
loglog(fxxS(fI),pxxS(fI)); zoom on;


stf = conv(stack,flipud(indiv_events(:,1057)));
[pxxS,fxxS] = pmtm(stf,{4,'trace'},4*nfft,Fs); %12164

fI = fxxS >= 3 & fxxS <= 40;
loglog(fxxS(fI),pxxS(fI)); zoom on;

load GGPSNR.mat
SNR = SNR.^0.5; nfft = 2^14;

U = fft(indiv_events(:,1057),nfft);
G = fft(stack,nfft);
stfWiener = fftshift(ifft((SNR.*conj(G).*U)./(SNR.*conj(G).*G + 1),[],'symmetric'));
hold on;
[pxxS,fxxS] = pmtm(stfWiener,{4,'trace'},nfft,Fs); %12164
fI = fxxS >= 3 & fxxS <= 40;
loglog(fxxS(fI),pxxS(fI)); zoom on;

%% wiener decon
clearvars -except stack indiv_events tw Fs
close all;
stack = stack / norm(stack);
nfft = 2^14;
G = fft(stack,nfft);
BMHL = 256;

load GGPSNR.mat
SNR = SNR.^0.5;
pxxS = NaN(1 + nfft/2,size(indiv_events,2));
figure('units','normalized','outerposition',[0 0 1 1]);
parfor i = 1:size(indiv_events,2)
    disp(i);
    U = fft(indiv_events(:,i),nfft);
    numerator = conj(G).*U;
    denominator = conj(G).*G;
    wienerDivision = fftshift(ifft((SNR.*conj(G).*U)./(SNR.*conj(G).*G + 1),[],'symmetric'));
    %pxxS_ = pwelch(wienerDivision,blackmanharris(BMHL),[],nfft,Fs); %12164
    pxxS_ = pmtm(wienerDivision,{4,'trace'},nfft,Fs);
    pxxS(:,i) = pxxS_;
end

%%
[~,fxxS] = pwelch(pxxS(:,1),blackmanharris(BMHL),[],nfft,Fs); %12164
fI = fxxS >= 3 & fxxS <= 50;
loglog(fxxS(fI),pxxS(fI,:)); zoom on;
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
