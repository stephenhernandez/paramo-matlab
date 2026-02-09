function S = scaleWaveforms(S,scalar,rmsFlag)
%
% scaleWaveforms return scaled waveforms
%
% S = scaleWaveforms(S,scalar,rmsFlag)
% S = scaleWaveforms(S,1,true) <-- scales each trace by 1*rms
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019

%%
if nargin < 2
    scalar = 1e-3;
end

if nargin < 3
    rmsFlag = false;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
if isscalar(scalar)
    scalar = scalar*ones(lS,1);
end

lscalar = length(scalar);
if lS ~= lscalar
    fprintf(2,"input dimension mismatch: length of scalar (%d) not same as length of structure (%d).\n",...
        lscalar,lS);
    return;
end

%%
refs = pull(S,"ref");
goodI = ~isnat(refs);
sumGood = sum(goodI);

%%
if ~sumGood
    fprintf(1,"no valid traces to scale. doing nothing\n");
    return;
end

%%
scalar = scalar(goodI);
S_ = S(goodI);
if ~rmsFlag
    for i = 1:sumGood
        d = S_(i).d;
        d = d*scalar(i);
        S_(i).d = d;

        %%
        [minVals,maxVals,meanVals] = minmaxmean(d);
        S_(i).depmin = minVals(1);
        S_(i).depmax = maxVals(1);
        S_(i).depmen = meanVals(1);
    end
    S(goodI) = S_;
    S = reshape(S,sizeS);
    return;
end

for i = 1:sumGood
    d = S_(i).d;
    dI = isfinite(d);
    rmsValue = rssq(d(dI));
    if rmsValue == 0
        rmsValue = 1;
    end
    d(dI) = d(dI)./(scalar(i).*rmsValue);
    S_(i).d = d;

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S_(i).depmin = minVals(1);
    S_(i).depmax = maxVals(1);
    S_(i).depmen = meanVals(1);
end