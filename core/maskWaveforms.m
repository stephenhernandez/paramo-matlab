function S = maskWaveforms(S,tStart,tEnd,fillValue)
%
% maskWaveforms return structure with waveforms clobbered out within window
% specified and with specified fill value (0 default)
%
% S = maskWaveforms(S,tStart,tEnd,fillValue)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 4
    fillValue = 0;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    d = S(i).d;
    ref = S(i).ref;
    delta = S(i).delta;
    si = t2i(tStart,ref,delta);
    ei = t2i(tEnd,ref,delta);
    tI = (si:ei);
    d(tI) = fillValue;
    S(i).d = d;

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S(i).depmin = minVals;
    S(i).depmax = maxVals;
    S(i).depmen = meanVals;
end
S = reshape(S,sizeS);