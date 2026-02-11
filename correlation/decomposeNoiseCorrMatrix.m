function [caus,acaus,symStack] = decomposeNoiseCorrMatrix(dayStack,detrendFlag,verboseFlag,normFlag)
%
% not the most optimized code ive ever written
%

if nargin < 2
    detrendFlag = false; %detrend only performed when normflag triggered
end

if nargin < 3
    verboseFlag = false;
end

if nargin < 4
    normFlag = false; %detrend only performed when normflag triggered
end

%%
lDayStack = size(dayStack,1);
m = (lDayStack+1)/2;
lags = (-m+1:m-1)';

% %%
% if mod(lDayStack,2) == 0 %if even
%     midpoint = lDayStack/2 + 1;
%     startInd1 = 2;
% else
%     midpoint = (lDayStack+1)/2;
%     startInd1 = 1;
% end
% 
% %%
% if verboseFlag
%     fprintf('midpoint = %d\n',midpoint);
% end

%%
% endInd1 = midpoint;
% startInd2 = midpoint;
% endInd2 = lDayStack;

nLagsI = lags <= 0;
pLagsI = lags >= 0;

%acaus = flipud(dayStack(startInd1:endInd1,:));
%caus = dayStack(startInd2:endInd2,:);

acaus = flipud(dayStack(nLagsI,:));
caus = dayStack(pLagsI,:);

if normFlag
    %detrend only performed when normflag triggered
    acaus = normalizeWaveforms(acaus,detrendFlag);
    caus = normalizeWaveforms(caus,detrendFlag);
end
symStack = (caus + acaus)/2;