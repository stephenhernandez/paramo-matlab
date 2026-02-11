function S = nanGapWaveforms(S,value,verboseFlag)
%
% nanGapWaveforms fill gaps with NaNs
%
% S = nanGapWaveforms(S)
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Sep 16, 2020
%

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
    gapFlag = S_.gapFlag;
    if ~gapFlag
        if verboseFlag
            fprintf("no gaps\n");
        end
        S_.gapInfo = [];
        S(i) = S_;
        continue;
    end

    %%
    d = double(S_.d);
    gapInfo = S_.gapInfo;
    gapStart = gapInfo(:,1);
    gapEnd = sum(gapInfo,2)-1;
    nGaps = size(gapInfo,1);

    %%
    if verboseFlag
        delta = S_.delta;
        tref = S_.ref;
        gapDurs = gapInfo(:,2);
        sumGaps = sum(gapDurs);
        gapDurSeconds = sumGaps*delta;
        eTot = seconds(S_.e);
        averageGapLength = delta*mean(gapDurs,"omitnan");
        fprintf(1,"%d gap(s) spanning %f seconds out of %f total seconds in trace, average gap: %f\n",...
            nGaps,gapDurSeconds,eTot,averageGapLength);
    end

    %%
    for j = 1:nGaps
        si = gapStart(j);
        ei = gapEnd(j);
        d(si:ei) = value;
        if verboseFlag && nGaps < 10
            gapStartTime = i2t(si,tref,delta);
            gapEndTime = i2t(ei,tref,delta);
            gapDurs_ = gapDurs(j);
            fprintf(1,"Gap %d: start: %s, end: %s, gap length in samples: %d\n",...
                j,gapStartTime,gapEndTime,gapDurs_);
        end
    end

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S_.depmin = minVals;
    S_.depmax = maxVals;
    S_.depmen = meanVals;

    %%
    S_.d = d;
    S(i) = S_;
end
S = reshape(S,sizeS);