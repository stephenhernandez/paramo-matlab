function d = demean(d,medFlag)
%
% demean return waveform with zero mean
%
% d = demean(d,medFlag)
% medflag: optional flag to remove median instead of mean (default=false)
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
if nargin < 2
    medFlag = false;
end

%%
centroid = mean(d,1,"omitnan");
if medFlag
    centroid = median(d,1,"omitnan");
end

%% for versions of matlab after 2016b, implicit expansion is done
d = d - centroid;

%% uncomment for versions before 2016b
% r = size(X,1);
% d = d - repmat(centroid,r,1);