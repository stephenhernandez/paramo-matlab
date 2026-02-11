clear; close all; clc;

tmaster = [];
dmaster = [];
maxp_master = tmaster;

S = loadWaveforms(datetime(2018,06,26),1,"SN13","HHZ","9D","",1,1,'~/data/iguana/BROADBAND/');
S = resampleWaveforms(S,50);

Fs = round(1/S.delta);
t = getTimeVec(S,0);
S = filterWaveforms(S,0.5);

nfft = 2048;
nOverlap = 0;
[dcut,startIndex,endIndex,badFlag] = cutWindows(S.d,nfft,nOverlap,true);

dcut = fft(dcut);
dcut = dcut(1:nfft/2+1,:);
dcut = abs(dcut).^2;
maxp = sum(dcut);
dcut = (1/(Fs*nfft)) * dcut;
dcut(2:end-1,:) = 2*dcut(2:end-1,:);

ttmp = t(startIndex);
maxp = maxp/(Fs*nfft);
maxp = sqrt(maxp);


tmaster = [tmaster ttmp];
dmaster = [dmaster dcut];
maxp_master = [maxp_master maxp];

decFact = 1;
nfftOrig = 2048;
nfft = nfftOrig;
nfft = nfft/decFact;
freq = 0:Fs/nfft:Fs/2;
newlen = length(freq);
if decFact > 1
    for i = 1:size(dmaster,2)
        disp(i)
        dtmp = dmaster(:,i);
        dtmp = decimate(dtmp,decFact);
        dmaster(1:newlen,i) = dtmp;
    end
end

dmaster = dmaster(1:newlen,:);
[maxdmaster,mI] = max(dmaster);

dnorm = bsxfun(@rdivide,dmaster,maxdmaster);
normp = bsxfun(@rdivide,cumsum(dmaster),sum(dmaster));
for i = 1:size(normp,2)
    mI(i) = find(normp(:,i) >= 0.5,1);
end
maxFreq = freq(mI);
boxlen = 45;
box = ones(boxlen,1)/boxlen;
smoothMaxFreq = fftfilt(box,maxFreq);

%%
figure;
imagesc(tmaster,freq,10*log10(dmaster)); axis xy; colorbar; caxis([0 70]);
hold on;
plot(tmaster,smoothMaxFreq,'k.');
