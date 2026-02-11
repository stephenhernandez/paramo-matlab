function [snr,fxx] = freqDomSNR(S,refTime,winlen,cmp,method,diffFlag,plotFlag)
if nargin < 4; cmp = 1; end
if nargin < 5; method = 1; end
if nargin < 6; diffFlag = false; end
if nargin < 7; plotFlag = false; end

if diffFlag
    S(cmp) = differentiateWaveforms(S(cmp));
end

Fs = 1/S(cmp).delta;
ScutNoise = cutWaveforms(S(cmp),refTime,-seconds(winlen),0);
ScutSignal = cutWaveforms(S(cmp),refTime,0,seconds(winlen));
noise = detrend(ScutNoise.d);
signal = detrend(ScutSignal.d);

nfft = nextpow2(winlen*Fs);
nfft = 2^(nfft+2);
if ~method
    [pxxN,fxx] = pmtm(noise,3,nfft,Fs);
    pxxS = pmtm(signal,3,nfft,Fs);
else
    [pxxN,fxx] = pwelch(noise,[],[],nfft,Fs);
    pxxS = pwelch(signal,[],[],nfft,Fs);
end

fI = fxx > 0;
pxxS = abs(pxxS(fI));
pxxN = abs(pxxN(fI));
snr = pxxS./pxxN;
snr = sqrt(snr);
fxx = fxx(fI);

if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    tnoise = getTimeVec(ScutNoise);
    tsignal = getTimeVec(ScutSignal);
    subplot(3,2,[1 2]);
    plot(tnoise,noise); hold on;
    h2 = plot(tsignal,signal);
    axis tight;
    ha(1) = subplot(3,2,3);
    plot(fxx,pxxN);
    ha(2) = subplot(3,2,4);
    plot(fxx,pxxN); hold on;
    plot(fxx,pxxS); %,'Color',h2.Color);
    ha(3) = subplot(3,2,[5 6]);
    semilogy(fxx,snr,'k'); hold on;
    plot([min(fxx) max(fxx)],[1 1],'--','linewidth',2,'Color',[0.5 0.5 0.5]);
    linkaxes(ha,'x');
    ha(3).XScale = 'log';
    zoom on; grid on;
end