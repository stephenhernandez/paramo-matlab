clear; %close all; clc;

Fs = 64;
secDur = 128;
lowestFreq = 0.01;
highestFreq = 25;
nfreqs = 1+1024*4;
freqVec = logspace(log10(lowestFreq),log10(highestFreq),nfreqs)';

%tStart = datetime(2009,01,01);
%tEnd = datetime(2011,12,31);
tStart = datetime(2006,01,180);
tEnd = datetime(2018,12,31);
%tStart = datetime(2009,04,28);
%tEnd = datetime(2009,04,28);
%tStart = datetime(2010,05,01);
%tEnd = datetime(2010,05,31);

dayVec = (tStart:tEnd)';
ldays = length(dayVec);

winlen = secDur.*Fs;
nfft = secDur*Fs;
w = kaiser(nfft/2,30);
interpMethod = "pchip";

kstnm = ["BRUN";"BULB";"BPAT";"BMAS";"BBIL"];
kcmpnms = "BHE";
lkstnms = length(kstnm);
hvMain = NaN(nfreqs,ldays,lkstnms);

nOverlap = 1/4;
nfft3 = nfft*2;
nfft4 = nfft+1;
tbw = 3;
[e,v] = dpss(nfft,tbw,5);
fxx = (0:Fs/nfft3:Fs/2)';

%%
lfc = 1/50;
hfc = 24;
units = 'vel';
detrendFlag = true;
R = singleSNCLFreqResponse(["BULB","BHZ","EC",""],...
    datetime(2010,01,01),...
    datetime(2010,01,02),0,Fs,units);
zeroes = R.zeros;
poles = R.poles;
constant = R.constant;
npoles = 4;
direction = true;

pwelchFlag = true;
exFlag = true;
winPerDay = NaN(ldays,1);
parfor i = 1:ldays
    day_ = dayVec(i);
    S = loadWaveforms(day_,1,kstnm,kcmpnms,"EC","");
    lS = length(S);
    if lS ~= lkstnms
        fprintf('not enough traces\n');
        continue;
    end
    if any(isnat(pull(S,'ref')))
        fprintf('no data on day: %s\n',datestr(day_));
        continue;
    end

    S = differentiateWaveforms(S);
    S = resampleWaveforms(S,Fs);
    S = interpolateWaveforms(S);
    S = transfer(S,zeroes,poles,constant,lfc,hfc,npoles,direction);
    S = scaleWaveforms(S,1e9);  % convert to nanometers
    S = nanGapWaveforms(S,0);   % in theory removes non-finites
    S = padWaveforms(S);
    S = syncWaveforms(S);

    minB = min(pull(S,'ref'));
    maxE = max(pull(S,'ref')+pull(S,'e'));

    tStart_ = dateshift(minB,'end','minute');
    tEnd_ = dateshift(maxE,'start','minute');

    dur_ = seconds(tEnd_-tStart_);
    if dur_ < secDur
        fprintf('snippet too short, skipping: %s\n',datestr(day_));
        continue;
    end

    S = cutWaveforms(S,tStart_,0,dur_);
    refs = pull(S,'ref');
    badI = isnat(refs);
    S(badI) = [];
    lS = length(S);

    if lS ~= lkstnms
        fprintf('not enough traces\n');
        continue;
    end

    try
        dcut = [];
        for j = 1:lS
            d_ = double(pull(S(j)));
            if exFlag
                dcut_ = cutWindows(d_,winlen,nOverlap,true);
            else
                dcut_ = doAutoCorrFreqDom(cutWindows(d_,winlen,nOverlap,true));
            end
            [nfft2,nWindows] = size(dcut_);
            dcut_ = dcut_(:);
            dcut = [dcut dcut_];
        end
        maxIndices = nfft+nfft2*(0:nWindows-1)';

        if exFlag
            rssqDcut = fftfilt(ones(winlen,1),dcut.^2);
            normers = median(rssqDcut(maxIndices,:),2,"omitnan");
            minNorm = min(rssqDcut(maxIndices,:),[],2);
        else
            normers = median(dcut(maxIndices,:),2,"omitnan");
            minNorm = min(dcut(maxIndices,:),[],2);
        end
        %figure(); semilogy(minNorm,'.'); zoom on; hold on; semilogy(normers,'.');
        badWindowsI = minNorm <= 2e10;

        dcut2 = [];
        for j = 1:lS
            dcut_ = dcut(:,j);
            dcut_ = reshape(dcut_,[nfft2,nWindows]);
            dcut2 = [dcut2; dcut_];
        end
        dcut = dcut2;
        dcut(:,badWindowsI) = [];
        normers(badWindowsI) = [];
        dcut = dcut ./ normers';

        normersOrig = normers;
        p2rms = peak2rms(dcut)';
        medAmp = median(p2rms,'omitnan');

        badWindowsI = normersOrig >= median(normersOrig,"omitnan") | p2rms >= medAmp;
        normers(badWindowsI) = [];
        dcut(:,badWindowsI) = [];

        nWindows = size(dcut,2);
        if ~nWindows
            fprintf('no windows to process\n');
            continue;
        end

        if lS > 1
            if exFlag
                ei = (winlen)*(1:lS)';
                si = 1 + (winlen)*(0:lS-1)';
                dcutNew = zeros(winlen,nWindows,lS);
            else
                ei = (2*winlen-1)*(1:lS)';
                si = 1 + (2*winlen-1)*(0:lS-1)';
                dcutNew = zeros(2*winlen-1,nWindows,lS);
            end

            for j = 1:lS
                dcut_ = dcut(si(j):ei(j),:);
                dcutNew(:,:,j) = dcut_;
            end
            dcut = dcutNew;
        end

        hv_ = NaN(nfreqs,1,lkstnms);

        
        for j = 1:lS
            dcut_ = squeeze(dcut(:,:,j));
            if pwelchFlag
                %[pxx1,fxx] = pwelch(dcut_,w,[],nfft3,Fs);
                [pxx1,fxx] = pmtm(dcut_,e,v,nfft2,Fs); %,'DropLastTaper',false);

                pxx1 = median(sqrt(Fs.*nfft3.*pxx1),2,"omitnan");
                hv = interp1(fxx,pxx1,freqVec,interpMethod,"extrap");
                hv_(:,1,j) = hv;
            else
                pxx1 = fft(dcut_,nfft3);
                pxx1 = (1/(Fs*nfft3))*abs(pxx1(1:nfft4,:)).^2;
                pxx1(2:end-1) = 2*pxx1(2:end-1);
                pxx1 = median(sqrt(Fs.*nfft3.*pxx1),2,"omitnan");
                hv = interp1(fxx,pxx1,freqVec,interpMethod,"extrap");
                hv_(:,1,j) = hv;
            end
        end
        hvMain(:,i,:) = hv_;
        winPerDay(i) = nWindows;
        fprintf('done with: <strong>%s</strong>\n',datestr(day_));
    catch ME
        fprintf(2,'Potentially fatal error on %s\n',datestr(day_));
        warning(ME.message);
        continue;
    end
end

%%
cd ~/research/now/hvsr/
save('tunguVerticalAmpSpectrumAcc_BHE_PMTM128sec25PrctOverlap','hvMain','freqVec','winPerDay','dayVec','kstnm');
