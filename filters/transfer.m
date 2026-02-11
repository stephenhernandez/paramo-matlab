function S = transfer(S,zero,pole,constant,lfc,hfc,npoles,DECONFLAG,...
    verboseFlag,CANUSEGPU,zeroPhaseFlag)
%
% code to apply _same_ poles/zeros to _all_ traces
% assumes the sampling rate is the same across all traces
%

%%
if nargin < 5; lfc = -inf; end
if nargin < 6; hfc = -inf; end
if nargin < 7; npoles = 4; end
if nargin < 8; DECONFLAG = true; end
if nargin < 9; verboseFlag = false; end
if nargin < 10; CANUSEGPU = false; end
if nargin < 11; zeroPhaseFlag = false; end

%%
npts = pull(S,"npts");
uniqnpts = unique(npts);
luniqnpts = length(uniqnpts);

for j = 1:luniqnpts
    uniqnpts_ = uniqnpts(j);
    nI = uniqnpts_ == npts;
    S_ = S(nI);
    nn = sum(nI);

    %%
    Fs = 1/S_(1).delta;
    d = pull(S_);
    norig = length(d);

    if DECONFLAG
        if verboseFlag
            fprintf("getting response information...\n");
            fprintf("deconvolving instrument\n");
        end
        constant = 1/constant;
        pole_old = pole;
        pole = zero;
        zero = pole_old;
    else
        if verboseFlag
            fprintf("getting response information...\n");
            fprintf("applying new instrument (convolving)\n");
        end
    end

    %%
    Hdecon = cmplxResp(norig,zero,pole,constant,Fs);
    if CANUSEGPU
        Hdecon = gpuArray(Hdecon);
        d = gpuArray(d);
    end

    %%
    if any(isfinite([lfc hfc]))
        if verboseFlag
            fprintf("applying filter: %f %f\n",lfc,hfc);
        end
        Hbu = freqOperator(norig,lfc,hfc,Fs,npoles);
        if CANUSEGPU
            Hbu = gpuArray(Hbu);
        end
        if zeroPhaseFlag
            Hbu = Hbu.*conj(Hbu);
        end
        Hdecon = Hdecon.*Hbu;
    end

    D = fft(d);
    D = D.*Hdecon; % implicit expansion
    d = ifft(D,[],1,"symmetric");

    %%
    for i = 1:nn
        d_ = d(:,i);
        S_(i).d = d_;
        S_(i).npts = norig;
        S_(i).delta = 1/Fs;
        S_(i).e = seconds((norig-1)/Fs);

        %%
        [minVals,maxVals,meanVals] = minmaxmean(d_);
        S_(i).depmin = minVals;
        S_(i).depmax = maxVals;
        S_(i).depmen = meanVals;
    end
    S(nI) = S_;
end