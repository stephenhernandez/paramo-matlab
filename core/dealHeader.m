function [S,npts] = dealHeader(S,d,newFs,newRef)
S = S(1);
if nargin < 3
    newFs = 1/S.delta;
end

if nargin < 4
    newRef = S.ref;
end

%% redefine important parameters
S.ref = newRef;
S.gapFlag = false;
S.gapInfo = [];

%%
npts = length(d);
dfin = isfinite(d);
sumnan = sum(~dfin);
if sumnan
    gapStart = find(diff(dfin)<0) + 1;
    gapEnd = find(diff(dfin)>0);
    lgs = size(gapStart,1);
    lge = size(gapEnd,1);
    if lgs > lge || ~dfin(end)
        gapEnd = [gapEnd; npts];
        lge = lge + 1;
    end

    if lge > lgs || ~dfin(1)
        gapStart = [1; gapStart];
        lgs = lgs + 1;
    end

    %%
    lE = min([lgs lge]);
    gapStart = gapStart(1:lE);
    gapDur = gapEnd(1:lE) - gapStart + 1;
    S.gapInfo = [gapStart gapDur];
    S.gapFlag = true;
end

%%
S.npts = npts;
S.e = seconds((npts-1)/newFs);
S.delta = 1/newFs;
S.d = d;

%%
[minVals,maxVals,meanVals] = minmaxmean(d);
S.depmin = minVals;
S.depmax = maxVals;
S.depmen = meanVals;