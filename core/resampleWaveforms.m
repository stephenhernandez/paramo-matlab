function Sout = resampleWaveforms(S,newFs,verboseFlag)
%
% resampleWaveforms return structure with resampled waveforms
%
% Sout = resampleWaveforms(S,newFs,verboseFlag)
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019
%

%%
if nargin < 1
    disp();
    Sout = populateWaveforms();
    return;
end

if nargin < 2
    newFs = 100;
end

if nargin < 3
    verboseFlag = false;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
Sout = S;
for i = 1:lS
    S_ = S(i);
    if isnat(S_.ref)
        Sout(i) = populateWaveforms();
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
    newd = resampleSH(d,delta,newFs,verboseFlag);
    Sout(i) = dealHeader(S_,newd,newFs);

    %%
    lnew = length(newd);
    dfin = isfinite(newd);
    sumnan = sum(~dfin);
    if ~sumnan
        % no more gaps
        Sout(i).gapFlag = false;
        Sout(i).gapInfo = [];
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
    Sout(i).gapInfo = [gapStart(1:lE) gapDur];
    Sout(i).gapFlag = true;
end
Sout = reshape(Sout,sizeS);
