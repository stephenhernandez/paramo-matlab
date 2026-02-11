function S = powWaveforms(S,exponent)
%
% powWaveforms return structure with waveforms raised to a power
%
% S = powWaveforms(S,exponent)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Friday, Oct 22, 2021

%%
if nargin < 2
    exponent = 2;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    S_ = S(i);
    d = pull(S_);

    %%
    nanI = isnan(d);
    d_ = d(~nanI);
    if exponent <= 1
        d_ = abs(d_);
    end

    if exponent > 0 %if not greater than zero, essentially return abs
        fprintf('exponent is: %f\n',exponent)
        d_ = d_.^exponent;
    end
    d(~nanI) = d_;

    %%
    S(i).d = d;

    %%
    [minVals,maxVals,meanVals] = minmaxmean(d);
    S(i).depmin = minVals;
    S(i).depmax = maxVals;
    S(i).depmen = meanVals;
end
S = reshape(S,sizeS);
