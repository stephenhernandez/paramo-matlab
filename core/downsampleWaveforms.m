function S = downsampleWaveforms(S,N)
%
% downsampleWaveforms return structure with downsampled waveforms
% this is pure downsampling, with no anti-alias filter applied beforehand
%
% S = downsampleWaveforms(S,N)
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Dec 26, 2023
%

%%
if nargin < 2
    N = 2;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    ref_ = S_.ref;
    if isnat(ref_)
        S(i) = populateWaveforms(); %why do i clobber?
        continue;
    end

    gapFlag = S_.gapFlag;
    if gapFlag
        S_ = nanGapWaveforms(S_,NaN);
    end

    %%
    d = double(S_.d);

    %%
    delta = S_.delta;
    phaseOffset = 0;
    newd = downsample(d,N,phaseOffset);
    newdelta = delta*N;
    newFs = 1./newdelta;
    S(i) = dealHeader(S_,newd,newFs);

    %%
    lnew = length(newd);
    dfin = isfinite(newd);
    sumnan = sum(~dfin);
    if ~sumnan
        % no more gaps
        S(i).gapFlag = false;
        S(i).gapInfo = [];
        continue;
    end

    %%
    gapStart = find(diff(dfin)<0) + 1;
    gapEnd = find(diff(dfin)>0);

    lgs = size(gapStart,1);
    lge = size(gapEnd,1);
    if lgs > lge || ~dfin(end)
        gapEnd = [gapEnd; lnew];
        lge = lge + 1;
    end

    if lge > lgs || ~dfin(1)
        gapStart = [1; gapStart];
        lgs = lgs + 1;
    end

    %%
    lE = min([lgs lge]);
    gapDur = gapEnd(1:lE) - gapStart(1:lE) + 1;
    S(i).gapInfo = [gapStart(1:lE) gapDur];
    S(i).gapFlag = true;
end
S = reshape(S,sizeS);
