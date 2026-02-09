function S = envelopeWaveforms(S,logFlag)
%
% envelopeWaveforms return struct with original waveforms converted
% to envelopes of those waveforms
%
% S = envelopeWaveforms(S,logFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    logFlag = false;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    ref = S_.ref;
    if isnat(ref)
        continue;
    end
    
    %%
    d = S_.d;
    try
        d = abs(hilbert(d(:)));
    catch
        S_ = interpolateWaveforms(S_);
        d = S_.d;
        d = abs(hilbert(d(:)));
    end
    
    %%
    if logFlag
        d = log10(d);
    end
    S_.d = d;
    
    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S_.depmin = minVals;
    S_.depmax = maxVals;
    S_.depmen = meanVals;
    
    %%
    S(i) = S_;
end
S = reshape(S,sizeS);