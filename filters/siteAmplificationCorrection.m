function S = siteAmplificationCorrection(S,amplificationCurve)

lS = length(S);
freqVec = amplificationCurve(:,1);
ampFactors = amplificationCurve(:,2);

for i = 1:lS
    %tic;
    S_ = nanGapWaveforms(detrendWaveforms(S(i)),0);
    d = S_.d;
    ld = length(d);
    nfft = 2^nextpow2(ld);

    D = fft(d,nfft);
    H = ones(nfft,1);

    %% # positive frequencies
    oddFlag = mod(nfft,2);
    if oddFlag
        midFreq = (nfft+1)/2;
    else
        midFreq = 1 + (nfft/2);
    end

    %%
    extrapValue = 1;
    Fs = 1/S_.delta;
    Fny = Fs/2;

    %
    f = (0:midFreq-1)'*Fny/midFreq;
    newAmpCurve = interp1(freqVec,ampFactors,f,"linear",extrapValue);
    H(1:midFreq) = newAmpCurve;

    % exploit symmetry
    if oddFlag
        H(midFreq+1:end) = flipud(newAmpCurve(2:midFreq));
    else
        H(midFreq+1:end) = flipud(newAmpCurve(2:midFreq-1));
    end

    if ~oddFlag
        % the Nyquist should be real
        H(midFreq) = 1; %abs(H(nfreq));
    end

    D = D./H;
    df = ifft(D,[],1,'symmetric');
    df = df(1:ld);
    S_ = dealHeader(S_,df);
    S(i) = S_;
    %toc;
end