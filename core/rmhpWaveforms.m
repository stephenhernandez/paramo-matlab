function [S,dlp] = rmhpWaveforms(S,secDur,tw)
%
% rmhpWaveforms return structure with filtered waveforms
% perform low-pass filter with indicated secDur, then substract that trace
% from the raw. Essentially returns a high-pass version.
%
% S = rmhpWaveforms(S,secDur)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Sunday, Oct 02, 2022

%%
if nargin < 2
    fprintf(2,'not enough inputs, returning empty structure\n')
    S = populateWaveforms();
    return;
end

if nargin < 3
    tw = false;
end

%%
sizeS = size(S);
S = S(:);
deltas = pull(S,'delta');
Fs = 1./deltas;
fI = Fs >= 1;
Fs(fI) = round(Fs(fI));

%% luniqnpts
uniqFs = unique(Fs);
luniqFs = sum(isfinite(uniqFs));

%%
for i = 1:luniqFs
    uniqFs_ = uniqFs(i);
    winlen = secDur*uniqFs_;
    Hd = ones(winlen,1)/winlen;
    fI = Fs == uniqFs_;
    S_ = S(fI);                             %waveforms with a uniq Fs

    %%
    npts = pull(S_,'npts');
    uniqnpts = unique(npts);
    luniqnpts = length(uniqnpts);

    for j = 1:luniqnpts
        uniqnpts_ = uniqnpts(j);
        if uniqnpts_ < winlen
            fprintf("number of pts in current trace(s): %d, number of points in filter: %d\n",uniqnpts_,winlen);
            fprintf("trace is too short, skipping!\n");
            continue;
        end
        Hbu = fft(Hd,uniqnpts_);
        Hbu = Hbu.*conj(Hbu);
        nI = find(uniqnpts_ == npts);
        nn = length(nI);

        %%
        S__ = S_(nI);

        %% taper if desired
        if tw
            S__ = taperWaveforms(S__,tw);
        end
        d0 = double(pull(S__));

        %% re-populate with filtered versions
        for k = 1:nn
            d_ = d0(:,k);
            goodI = isfinite(d_);
            d_(~goodI) = 0;

            %% get low-pass here
            dlp = Hbu.*fft(d_);
            dlp = ifft(dlp,[],1,"symmetric"); 

            %% remove low-pass here
            dOrig = d_;
            d_ = dOrig - dlp;
            d0(goodI,k) = d_(goodI);

            S__(k).d = d0(:,k);
            [minVals,maxVals,meanVals] = minmaxmean(d0(:,k));
            S__(k).depmin = minVals;
            S__(k).depmax = maxVals;
            S__(k).depmen = meanVals;
        end
        S_(nI) = S__;
    end
    S(fI) = S_;
end

%%
S = reshape(S,sizeS);