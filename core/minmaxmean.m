function [minVals,maxVals,meanVals] = minmaxmean(d,verboseFlag)
%
% [minVals,maxVals,meanVals] = minmaxmean(d,verboseFlag)
%
% minmaxmean get min, max, and mean of columns of matrix
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019
% updated 22 sep 2021
% updated 17 may 2022
%

%%
if nargin < 2
    verboseFlag = false;
end

%%
[nr,nc] = size(d);
if nr < 2
    if verboseFlag
        fprintf('d is a row vector, reshaping to column vector\n');
    end
    d = d';
    nc = 1;
end

%%
if verboseFlag
    fprintf('processing %d trace(s)\n',nc);
end

%%
minVals = min(d,[],1,"omitnan")';
maxVals = max(d,[],1,"omitnan")';
meanVals = mean(d,1,"omitnan")';
