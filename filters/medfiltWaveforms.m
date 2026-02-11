function S = medfiltWaveforms(S,N,zeroPhaseFlag,tw,diffFlag)
%
% medfiltWaveforms apply median filter to waveform structure
%
% S = medfiltWaveforms(S,N,tw,zeroPhaseFlag,diffFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    N = 3;
end

if nargin < 3
    zeroPhaseFlag = false;
end

if nargin < 4
    tw = false;
end

if nargin < 5
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
sizeS = size(S);
S = S(:);
lS = length(S);
for i = 1:lS
    d = double(pull(S(i)));
    
    %%
    nanI = isnan(d);
    d_ = d(~nanI);
    d_ = medfiltSH(d_,N,zeroPhaseFlag);
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
