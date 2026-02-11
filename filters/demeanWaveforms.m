function S = demeanWaveforms(S)
%
% demeanWaveforms returns structure with zero-mean waveforms
%
% S = demeanWaveforms(S)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    d = S(i).d;
    d = demean(d);

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S(i).depmin = minVals;
    S(i).depmax = maxVals;
    S(i).depmen = meanVals;
    S(i).d = d;
end
S = reshape(S,sizeS);