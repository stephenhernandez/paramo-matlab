function [overtoneFreq,harmStrength,fundPow,...
    overtonePow,interharmPow,overtoneI,powerOrig,pxx,fxx] = ...
    tremometer(dcut,Fs,nOvertones)

%
% [overtonePow, overtoneFreq, interharmPow, harmStrength] = tremometer(dcut,Fs,nOvertones)
%

%%
MINFRAC = 5e-2; %determined empirically
[winlen,nTimes] = size(dcut);
nfft = 2^(nextpow2(winlen)+1);
nw = 2;
nTapers = 2*nw-1;
overtoneWidth = 2*nTapers + 1;
[pxx,fxx] = pmtm(dcut,nw,nfft,Fs);

%% pre-allocate harmonic tables, nOvertones is a user-supplied argument
overtonePow = NaN(nTimes,nOvertones+1);
interharmPow = NaN(nTimes,nOvertones);
overtoneI = zeros(size(overtonePow));

%% Get an estimate of the fundamental frequency / amplitude through Harmonic Product Spectrum calculation
pxx = zpkFilter(pxx,-inf,1/overtoneWidth,1,1,1);
pxx2 = resample(pxx,1,2);
pxx3 = resample(pxx,1,3);
pxx4 = resample(pxx,1,4);
pxx5 = resample(pxx,1,5);
pxx6 = resample(pxx,1,6);

%%
ll = size(pxx6,1);
fxx = fxx(1:ll);
FREQMICROSEISM = 0.2;
microI = find(fxx >= FREQMICROSEISM,1);
pxx(1:microI,:) = 0;
pxx2(1:microI,:) = 0;
pxx3(1:microI,:) = 0;
pxx4(1:microI,:) = 0;
pxx5(1:microI,:) = 0;
pxx6(1:microI,:) = 0;

%%
psd = cat(3,pxx(1:ll,:),pxx2(1:ll,:),pxx3(1:ll,:),pxx4(1:ll,:),pxx5(1:ll,:),pxx6(1:ll,:));
psd(psd<0) = 0;
normer = sum(psd,1);
normer(normer <= 0) = 1;
psd = psd./normer;

%%
powerOrig = prod(psd,3);
normer = sum(powerOrig,1);
normer(normer <= 0) = 1;
power = powerOrig./normer;

%% zero crossings
pxx = squeeze(psd(:,:,1)); %pxx is normalized here
dprime = diff(pxx);
peakI = dprime(2:end,:);
zc = peakI.*dprime(1:end-1,:);

%%
indexConverter = ll*(0:nTimes-1)';
fundI = overtoneI(:,1);
for i = 1:nTimes
    pow_ = power(:,i);
    [PKS,LOCS] = findpeaks(pow_);

    if isempty(PKS)
        [~,i_] = max(pow_,[],1);
        fundI(i) = i_;
        continue;
    end
    i_ = find(PKS >= MINFRAC,1);
    i_ = max([i_,1]);
    fundI(i) = LOCS(i_);
end
fundPow = powerOrig(fundI+indexConverter); %testing...

%% get indices of the harmonic frequencies / amplitudes, then use to query frequency and amplitude for each harmonic
lIndex = overtoneI;
rIndex = lIndex;

%%
[rows,cols] = find(zc < 0 & peakI < 0);
prevI = zeros(size(fundI));
for i = 1:nOvertones+1
    freqI = fundI*i;
    for j = 1:nTimes
        thisI = freqI(j);                           %particular test index for this time period (column)
        rows_ = rows(cols == j);                    %return peak indices for this column (j)
        if isempty(rows_)
            continue;
        end
        [~,minI] = min(abs(rows_ - thisI),[],1);    %find index to nearest peak
        freqI(j) = rows_(minI)+1;
    end
    interWidth_ = freqI - prevI;
    fudgeFactor = 0.3*interWidth_;
    lIndex(:,i) = round(prevI + fudgeFactor);
    rIndex(:,i) = round(freqI - fudgeFactor);
    prevI = freqI;
    overtoneI(:,i) = freqI;
end

%%
powerWidth = rIndex - lIndex + 1;
overtoneFreq = fxx(overtoneI);

%%
pxxCum = cumsum(pxx);
for i = 1:nOvertones+1
    si = lIndex(:,i) + indexConverter;
    ei = rIndex(:,i) + indexConverter;
    interharmPow(:,i) = pxxCum(ei) - pxxCum(si);

    freqI = overtoneI(:,i);
    si = freqI + indexConverter - nTapers;
    ei = freqI + indexConverter + nTapers;
    overtonePow(:,i) = pxxCum(ei) - pxxCum(si);
end
overtonePow = overtonePow/overtoneWidth;
interharmPow = interharmPow./powerWidth;
harmStrength = overtonePow(:,1:nOvertones)./interharmPow(:,2:nOvertones+1);