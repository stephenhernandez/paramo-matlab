function S = interpolateWaveforms(S,value,verboseFlag)
%
% linear interpolation of waveforms
%
% S = interpolateWaveforms(S)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Sep 16, 2020

%%
if nargin < 2
    value = NaN;
end

if nargin < 3
    verboseFlag = false;
end

if isempty(value) || ~isnumeric(value)
    value = NaN;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);

%%
for i = 1:lS
    S_ = S(i);
    if ~S_.gapFlag
        if verboseFlag
            fprintf("no gaps, doing nothing\n'");
        end
        continue;
    end

    %%
    d = double(S_.d);
    ld = length(d);

    %%
    gapInfo = S_.gapInfo;
    gapStart = gapInfo(:,1);
    gapEnd = sum(gapInfo,2)-1;
    gapDur = gapInfo(:,2);

    %%
    lastGap = size(gapInfo,1);
    for j = 1:lastGap
        si = gapStart(j)-1;
        ei = gapEnd(j)+1;
        if si >= 1 && ei <= ld
            skippedLength = gapDur(j);
            X = [1; 2+skippedLength];
            V1 = d(si);
            V2 = d(ei);
            if ~isfinite(V2)
                fprintf("last point in this gap is not finite\n");
            end
            V = [V1; V2];
            Xq = (2:2+skippedLength-1)';
            gap = interp1(X,V,Xq);
            d(si+1:ei-1) = gap;
        else
            si = si+1;
            ei = ei - 1;
            d(si:ei) = value;
        end
    end

    %%
    S_.d = d;
    S_.npts = ld;
    S_.e = seconds((S_.npts - 1)*S_.delta);

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S_.depmin = minVals;
    S_.depmax = maxVals;
    S_.depmen = meanVals;
    S(i) = S_;
end
S = reshape(S,sizeS);