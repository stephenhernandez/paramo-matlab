function S = detrendWaveforms(S)
%
% detrendWaveforms returns structure with detrended waveforms
%
% [S,p] = detrendWaveforms(S)
%
% S: structure with detrended waveforms
% p: polynomial coefficients used to model linear trend
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019
% updated 02 JUL 2022
warning on verbose

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    gapFlag = S_.gapFlag;
    if isnat(S_.ref)
        continue;
    end

    %% case with no gaps...
    if ~gapFlag
        d = S_.d;
        d = detrend(d);
        S_ = dealHeader(S_,d);
        S(i) = S_;
        continue;
    end

    %% case where gaps exist
    npts = S_.npts;
    d = S_.d;
    gapInfo = S_.gapInfo;

    %%
    nGaps = size(gapInfo,1);
    gapStarts = gapInfo(:,1);

    newFs = 1/S_.delta;
    newRef = S_.ref;
    firstStart = gapStarts(1);
    if firstStart ~= 1
        %there are at least two distinct segments
        contStartIndex = [1; sum(gapInfo,2)];
        contEndIndex = [gapInfo(:,1)-1; npts];
    else
        contStartIndex = sum(gapInfo,2);
        if nGaps == 1
            %only one contiguous segment in entire trace
            contEndIndex = npts;
            dcontiguous = d(contStartIndex:contEndIndex);

            %t = (1:length(d))';
            %p_ = polyfit(t,d,1);
            %dSynth = p_(1)*t + p_(2);
            %d = d - dSynth;
            %S_.gapFlag = false;
            %S_.gapInfo = [];
            %newRef = i2t(contStartIndex,newRef,1./newFs);
            %S_ = dealHeader(S_,d,newFs,newRef);

            dcontiguous = detrend(dcontiguous);
            d(contStartIndex:contEndIndex) = dcontiguous;
            S_ = dealHeader(S_,d);
            S(i) = S_;
            continue;
        else
            contEndIndex = [gapInfo(2:end,1)-1; npts];
        end
    end
    nCont = length(contStartIndex);

    %%
    contCount = contEndIndex - contStartIndex + 1;
    for kk = 1:nCont
        count_ = contCount(kk);
        if count_ < 1
            continue;
        end

        si = contStartIndex(kk);
        ei = contEndIndex(kk);
        d_2 = d(si:ei);

        if isempty(d_2) || si >= ei
            continue;
        end

        if count_ >= 2^26
            warning("could not detrend, this snippet is too long");
            continue;
        end

        t = (1:length(d_2))';
        p_ = polyfit(t,d_2,1);
        dSynth = p_(1)*t + p_(2);

        d_2 = d_2 - dSynth;
        d(si:ei) = d_2;
    end

    %%
    S_ = dealHeader(S_,d,newFs,newRef);
    S(i) = S_;
end
S = reshape(S,sizeS);