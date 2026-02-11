function S = oneBitNorm(S)
%
% S = oneBitNorm(S)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Sep 25, 2019

%% assumed the data are already filtered

%% scale each trace individually
lS = length(S);
for i = 1:lS
    d = S(i).d;
    d = sign(d);
    S(i).d = d;
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S(i).depmin = minVals;
    S(i).depmax = maxVals;
    S(i).depmen = meanVals;
end