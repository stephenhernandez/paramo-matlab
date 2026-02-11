function [S,normer] = normalizeWaveforms(S,detrendFlag,norm1Flag,constantRatio)
%
% normalizeWaveforms thinly veiled version of scaleWaveforms
% [S,maxAmp] = normalizeWaveforms(S,detrendFlag,norm1Flag,constantRatio)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019
% CODE TESTED: 09FEB2026

%%
if nargin < 2
    detrendFlag = false;    % detrend flag (default false)
end

if nargin < 3
    norm1Flag = false;      % normalize to 1? (default false)
end

if nargin < 4
    constantRatio = false;  %hold ratos constant between different stations
end

%%
SisStruct = isstruct(S);
if SisStruct
    sizeS = size(S);
    Sorig = S(:);
    if detrendFlag
        Sorig = detrendWaveforms(Sorig);
    end
    S = pull(Sorig);
else
    fprintf("input data are not waveform structs. treat as matrix where columns are individual traces.\n");
    if detrendFlag
        S = detrend(S);
    end
end

%%
badI = ~isfinite(S);
S(badI) = 0;
if norm1Flag
    normer = max(abs(S));
else
    normer = rssq(S);
end

%%
normer(normer==0) = 1;
if ~SisStruct
    if constantRatio
        fprintf("OPTION: constantRatio valid only for structs");
    end
    S = S./normer;
    S(badI) = NaN;
    return;
end

%%

if ~constantRatio
    lS = size(S,2); %<== number of columns of actual data matrix... has nothing to do with structs...
    S = S./normer;
    S(badI) = NaN;
    for j = 1:lS
        Sorig(j) = dealHeader(Sorig(j),S(:,j));
    end
    S = reshape(Sorig,sizeS);
    return;
end

if sizeS(2) < 2
    fprintf("OPTION: constantRatio: there is only one trace... try again!\n");
    S = populateWaveforms(1);
    return;
end

% NOTE: Constant ratio is useful for normalizing ZNE components of single
% station by a single value
normer = NaN(sizeS(1),1);
S = reshape(Sorig,sizeS);
for i = 1:sizeS(1)
    %loop here
    S_ = S(i,:)';
    d_ = pull(S_);
    sizeD = size(d_);
    d_ = d_(:); %in tontext of 3 traces (ZNE), concatenate into 1 single trace...
    badI = ~isfinite(d_);
    d_(badI) = 0;
    if norm1Flag
        normer_ = max(abs(d_));
    else
        normer_ = rssq(d_);
    end
    normer(i,1) = normer_
    d_ = d_/normer_;
    d_(badI) = NaN;
    d_ = reshape(d_,sizeD);
    for j = 1:sizeD(2)
        S__ = S_(j);
        S__ = dealHeader(S__,d_(:,j));
        S_(j) = S__;
    end
    S(i,:) = S_';
end