%january11_2020
clear; close all; clc;
cd ~/research/now/pichincha_nxcorr/
load('ggp_2008_2020_PLD');

%%
stack = normalizeWaveforms(detrend(stack(1:nfft)));
G = fft(stack,nfft);
tau = 1/max(abs(G));

%%
Fs = 100;
u = indiv_events(1:nfft,12164);
tdum = (0:nfft-1)'/Fs;

%%
noiseSpectrum = fft(u(1:575),nfft);
newPhaseAngle = -pi + 2*pi*rand(nfft,1);
newNoiseSpectrum = abs(noiseSpectrum).*exp(1j*newPhaseAngle);
newNoiseSpectrum(1) = 0;
newNoiseSpectrum((nfft/2)+1) = 0;

figure();
subplot(311);
semilogy(abs(noiseSpectrum),'-');
hold on;
semilogy(abs(newNoiseSpectrum),'-','linewidth',1);

subplot(312);
plot(ifft(noiseSpectrum,[],'symmetric'),'-');
hold on;
plot(ifft(newNoiseSpectrum,[],'symmetric'),'-','linewidth',1);

subplot(313);
plot(angle(noiseSpectrum),'-');
hold on;
plot(angle(newNoiseSpectrum),'-','linewidth',1);
zoom on;

%%
figure();
a(1) = subplot(211);
plot(tdum,stack)
a(2) = subplot(212);
plot(tdum,ifft(100*G+newNoiseSpectrum/10,[],'symmetric'));
linkaxes(a,'x');
zoom on;
