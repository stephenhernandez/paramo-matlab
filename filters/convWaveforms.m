function S = convWaveforms(S,n,tw,diffFlag)
%
% convWaveforms return structure with waveforms convolved with either
% n-point boxcar, or user-defined filter coefficients
%
% S = convWaveforms(S,n,tw,diffFlag)
% n: n a scalar is n-point boxcar [box = ones(n,1)/n]
% n: n a vector is user-defined filter coefficients
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    n = 3;
end

if nargin < 3
    tw = false;
end

if nargin < 4
    diffFlag = false;
end

%%
if diffFlag
    S = differentiateWaveforms(S);
end

if tw
    S = taperWaveforms(S,tw);
end

%%
if isscalar(n)
    disp(['Applying ',num2str(n),'-point boxcar moving filter (moving average)']);
    box = ones(n,1)/n;
elseif isvector(n)
    disp('convolving with user-given vector');
    box = n;
else
    disp('no filter applied');
    return;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);

for i = 1:lS
    d = double(pull(S(i)));
    
    %%
    nanI = isnan(d);
    d_ = d(~nanI);
    d_ = fftfilt(box,d_);
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
