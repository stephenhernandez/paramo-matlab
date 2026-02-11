%january11_2020
clear; close all; clc;
cd ~/research/now/pichincha_nxcorr/
load('ggp_2008_2020_PLD');

%%
upSampleFactor = 4;
Fs = upSampleFactor*100;
nfft = nfft*upSampleFactor;
indiv_events = resample(indiv_events,upSampleFactor,1);
stack = normalizeWaveforms(nanmean(normalizeWaveforms(indiv_events(:,mI)),2));
stack = 10*normalizeWaveforms(detrend(10*stack(1:nfft)));

%%
pI = cc_all_vs_stack_2 >= 0.4;
tabs = tabs(pI);
p2p = p2p(pI);
p2rms = p2rms(pI);
indiv_events = indiv_events(:,pI);
maxAmpRMS = maxAmpRMS(pI);

%%
G = fft(stack,nfft);
tau = 1/max(abs(G));
noiseStopIndex = 575*upSampleFactor;
% snippet = normalizeWaveforms(indiv_events(noiseStopIndex:noiseStopIndex+Fs*5-1,:));
% stackSnippet = normalizeWaveforms(stack(noiseStopIndex:noiseStopIndex+Fs*5-1,:));
% [~,~,~,~,raw_shifts] = apply_vdcc(snippet);
% indiv_events  = apply_shifts(indiv_events,raw_shifts);

%%
noiseOrig = detrend(indiv_events(1:noiseStopIndex,:));
noiseSpectrum = fft(noiseOrig,nfft);
newPhaseAngle = -pi + 2*pi*rand(nfft,size(indiv_events,2));
newNoiseSpectrum = abs(noiseSpectrum).*exp(1j*newPhaseAngle);
newNoiseSpectrum(1,:) = 0;
newNoiseSpectrum((nfft/2)+1,:) = 0;
newNoise = ifft(newNoiseSpectrum,[],'symmetric');

%%
[pxxN,fxx] = pwelch(noiseOrig,[],[],nfft,Fs);
pxxS = pwelch(indiv_events(1:nfft,:),[],[],nfft,Fs);
fI = fxx > 0;
fxx = fxx(fI);
pxxN = pxxN(fI,:);
pxxS = pxxS(fI,:);
snr = pxxS./pxxN; %this gives values between 2.5 - 10 hz that i like
figure(); loglog(fxx,nanmedian(snr,2)); zoom on;

%%
tdum = (0:nfft-1)'/Fs;
figure();
a(1) = subplot(211);
plot(tdum,stack)
a(2) = subplot(212);
plot(tdum,ifft(10000*G+newNoiseSpectrum(:,1752),[],'symmetric'));
linkaxes(a,'x');
zoom on;

%%
close all;
th = 1/100;
stf1 = 1e4*getSTF(th,Fs);
tdum = (0:nfft-1)'/Fs;
figure();
a(1) = subplot(211);
plot(tdum,stack)
a(2) = subplot(212);
synth1 = conv(stf1,stack);
plot(tdum,normalizeWaveforms(synth1(1:nfft) + ifft(newNoiseSpectrum(:,1752),[],'symmetric')));
linkaxes(a,'x');
zoom on;

hold on;
th = 1/50;
stf2 = 2e4*getSTF(th,Fs);
synth2 = conv(stf2,stack);
plot(tdum,normalizeWaveforms(synth2(1:nfft) + ifft(newNoiseSpectrum(:,1752),[],'symmetric')));

hold on;
th = 1/25;
stf3 = 4e4*getSTF(th,Fs);
synth3 = conv(stf3,stack);%+ ifft(newNoiseSpectrum(:,1752),[],'symmetric');
plot(tdum,normalizeWaveforms(synth3(1:nfft) + ifft(newNoiseSpectrum(:,1752),[],'symmetric')));


%
clc
t0 = noiseStopIndex/Fs;
T = 1/10;
tdum = (0:nfft-1)'/Fs;
Gstar = conj(G);
n = 1;
maxN = 200;
stfs3 = zeros(nfft,maxN); %stfs3(1:length(stf3),1) = stfs3(1:length(stf3),1)+stf3;
u = synth3;
%u = u(th*Fs:end); 
u = detrend(u(1:nfft)); %1:nfft); 
%u = indiv_events(1:nfft,1752);
U = fft(u,nfft);
tdm = tdum - nfft/2/Fs;
while n <= maxN
    disp(n);
    f = stfs3(:,n);
    rs = conv(stack,f);
    urs = u - rs(1:nfft);
    g = fftshift(conv(tau*flipud(stack),urs));
    
    if n == 1
        figure(); plot(g); zoom on;
        g1 = g;
    end
    
    g = g(1:nfft);
    tI = tdum >= 0 & tdum <= T & g >= 0;
    g(~tI) = 0;
    n = n + 1;
    f = f+g;
    newMax = max(f);
    f = f/sum(f);
    stfs3(:,n) = newMax*f; %fftshift(g);
end


%%

