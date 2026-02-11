function S = differentiateWaveforms(S,ndiff,freqFlag,verboseFlag)
%
% differentiateWaveforms returns structure with differentiated waveforms
%
% S = differentiateWaveforms(S,ndiff,freqFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    ndiff = 1;
end

if nargin < 3
    freqFlag = false;
end

if nargin < 4
    verboseFlag = false;
end

%%
sizeS = size(S);
S = S(:);
lS = length(S);
fprintf(1,"PRO TIP: detrend _before_ differentiating...\n"); %added 21 JAN 2026, delete when this advice is internalized...
if verboseFlag
    if ndiff == 1
        fprintf('Differentiating %d trace(s) in time domain once\n',lS);
    else
        fprintf('Differentiating %d trace(s) in time domain %d times\n',lS,ndiff);
    end
end

%%
if ~freqFlag
    S = nanGapWaveforms(S);
end

%%
for i = 1:lS
    S_ = S(i);
    if isnat(S_.ref)
        if verboseFlag
            fprintf('no data\n');
        end
        continue;
    end
    %%

    if freqFlag
        gapFlag = S_.gapFlag;
        if gapFlag
            S_ = interpolateWaveforms(S_);
        end

        delta = S_.delta;
        dy = S_.d;

        %%
        n = 1;
        while n <= ndiff
            dy = freqDiff(dy,1/delta);
            n = n+1;
        end

        %%
        S_ = dealHeader(S_,dy);
        S_ = nanGapWaveforms(S_);
        S(i) = S_;
        continue;
    end

    %%
    delta = S_.delta;
    dy = S_.d;

    %%
    n = 1;
    while n <= ndiff
        d1 = dy(1);
        dy = diff(dy);
        dy = [d1; dy]./delta;
        n = n + 1;
    end

    %%
    S_ = dealHeader(S_,dy);
    S(i) = S_;
end
S = reshape(S,sizeS);