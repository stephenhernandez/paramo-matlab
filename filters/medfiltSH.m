function y = medfiltSH(x,n,zeroPhaseFlag)
%
% medfiltSH 1-d median filtering with or without phase lags
%
% y = medfiltSH(x,n,zeroPhaseFlag)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    n = 3;
end

if nargin < 3
    zeroPhaseFlag = false;
end

%% function for a causal median filter
y = x;
n = floor(n);
if n < 3
    disp('n must be at least 3 points, doing nothing')
    return;
end

[rx,cx] = size(x);
if rx < 3
    disp('input vector/matrix must have at least 3 elements')
    return;
end

%%
if zeroPhaseFlag
    % no phase lags present
    y = movmedian(y,n,'omitnan');
else
    % induce a phase lag
    nans = NaN(n-1,cx);
    y = [nans; y];
    y = movmedian(y,n,'omitnan','Endpoints','discard');
end
y = y(1:rx,:);
