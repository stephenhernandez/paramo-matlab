function [dcut,startIndex,endIndex,badFlag,nwindows] = cutWindows(d,...
    winlen,nOverlap,detrendFlag)
%
% cutWindows cut windows of matrix
%
% [dcut,startIndex,endIndex,badFlag,nwindows] = cutWindows(d,winlen,nOverlap,detrendFlag)
%

%
% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019
% Major rewrites: 02FEB2026
%

%%
badFlag = true;
dcut = d;
nwindows = 0;
startIndex = 1;
endIndex = startIndex;
if nargin < 2
    fprintf(2,"not enough input arguments\n");
    return;
end

if nargin < 3
    nOverlap = 0.5;
end

if nargin < 4
    detrendFlag = true;
end

if nOverlap < 1
    % if nOverlap < 1, assume we've specified a percentage
    nOverlap = round(winlen*nOverlap);
end

%%
if nOverlap > winlen
    fprintf(2,'overlap is too long.\n');
    return;
end

dlen = size(d,1);
if winlen > dlen
    fprintf(2,'requested cut window is too long.\n');
    return;
end

%%
badFlag = false;
stride = winlen-nOverlap;
nwindows = 1+floor((dlen-winlen)/stride);
startIndex = (1+(0:nwindows-1)*stride)';
idx = repmat(startIndex', winlen, 1) + repmat((0:winlen-1)', 1, nwindows);
endIndex = idx(end,:)';
dcut = d(idx);

%%
if ~detrendFlag
    return;
end
dcut = detrend(dcut);