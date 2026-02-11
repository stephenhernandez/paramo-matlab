function [caus,acaus,symStack,t,lags] = getDV(varargin)
%function [dVcaus,dVacaus,dVsymmetric,ccCaus,ccAcaus,ccSymmetric,dayData,t,lags] ...
%    = getDV(varargin)
%
% [dVcaus,dVacaus,dVsymmetric,ccCaus,ccAcaus,ccSymmetric,dayData,t,lags] = ...
%   getDV(tStart,tEnd,SNCL1,SNCL2,dW,lfc,hfc,newFs,secDur,maxLag,dailyFlag,pathToWaveformServerOrDirectory);
%

%
% Written by Stephen Hernandez, Instituto Geofisico, Quito, Ecuador
% Modified by Stephen Hernandez: 05 May 2020
%

%%
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2015,08,01),...    % tStart
    datetime(2015,08,31),...    % tEnd
    ["SN02","HHZ","9D",""],... 	% SNCL1
    ["SN06","HHZ","9D",""],... 	% SNCL2
    0,...                       % dW
    1,...                       % lfc
    4,...                       % hfc
    32,...                      % newFs
    2^12,...                    % secDur
    32,...                      % maxLag
    false,...                   % dailyFlag
    '~/rawdata/'};              % pathToWaveformServerOrDirectory

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,SNCL1,SNCL2,dW,lfc,hfc,newFs,secDur,maxLag,dailyFlag,...
    pathToWaveformServerOrDirectory] = deal(optsToUse{:});

%%
fprintf('%s.%s.%s.%s\n',SNCL1(1),SNCL1(3),SNCL1(4),SNCL1(2));
fprintf('%s.%s.%s.%s\n',SNCL2(1),SNCL2(3),SNCL2(4),SNCL2(2));
fprintf('\n');

%%
if isfloat(tEnd)
    days = tStart:tStart+tEnd-1;
else
    days = tStart:tEnd;
end
lDays = length(days);

%%
tw = 0.01;
totN = secDur*newFs;
mxl = totN - 1;
stackN = 2*mxl+1;
nOverlap = 0;
normFlag = true;
lags = (-mxl:mxl)';
lags = lags/newFs;
mI = abs(lags) <= maxLag;
sum_mI = sum(mI);
maxWindows = floor(86400/secDur);

%%
if dailyFlag
    getStack = true;
    dayData = NaN(stackN,lDays);
    parfor i = 1:lDays
        dayStack_ = getCF(days(i),SNCL1,SNCL2,newFs,dW,lfc,hfc,tw,totN,...
            nOverlap,normFlag,maxWindows,stackN,getStack,pathToWaveformServerOrDirectory);
        dayData(:,i) = dayStack_;
    end
    t = days;
    dayData = dayData(mI,:);
    goodI = find(sum(~isnan(dayData)));
    dayData = dayData(:,goodI);
    t = t(goodI);
    
    Nsmooth = 1; % for daily flag do not apply smoothing, leave each day as is
    if Nsmooth > 1
        box = ones(Nsmooth,1)/Nsmooth;
        dayData = flipud(convn(flipud(dayData'),box));
        dayData = dayData(Nsmooth:end,:);
        dayData = normalizeWaveforms(dayData');
    else
        dayData = normalizeWaveforms(dayData);
    end
else
    getStack = false;
    fprintf('Total number of days to process: %d\n',lDays);
    fprintf('Total number of points in truncated lag time vector (+/-): %d\n',sum_mI);
    fprintf('Number of windows per day: %d\n',maxWindows);
    
    numberElements = sum_mI*lDays*maxWindows;
    fprintf('Total number of elements: %d\n',numberElements);
    
    fprintf('Attempting to pre-allocate %d points...\n',numberElements);
    totSubWindows = sum_mI * maxWindows;
    dayData = NaN(sum_mI,lDays*maxWindows);
    
    fprintf('Attempting to pre-allocate same number (%d) of points again, just a different shape...\n',numberElements);
    dayData2 = NaN(totSubWindows,lDays);
    totalIndividualTimeSlices = maxWindows*lDays;
    si = (1:maxWindows:totalIndividualTimeSlices)';
    ei = flipud((totalIndividualTimeSlices:-maxWindows:maxWindows)');
    t = NaT(maxWindows,lDays);
    
    %%
    parfor i = 1:lDays %parfor
        disp(datestr(days(i)));
        [dayData_,t_] = getCF(days(i),SNCL1,SNCL2,newFs,dW,lfc,hfc,tw,totN,...
            nOverlap,normFlag,maxWindows,stackN,getStack,pathToWaveformServerOrDirectory);
        dayData_ = dayData_(mI,:);
        dayData_ = dayData_(:);
        dayData2(:,i) = dayData_;
        t(:,i) = t_;
    end
    
    %%
    fprintf('here, i am populating the dayData matrix (via a reshape command)\n');
    for i = 1:lDays
        dayData_ = dayData2(:,i);
        dayData(:,si(i):ei(i)) = reshape(dayData_,sum_mI,maxWindows);
    end
    t = t(:);
    
    %%
    clear dayData2 %free up some memory!!
    Nsmooth = max([4 2^(nextpow2(maxWindows)-2)]); %-3 --> approximately 4.5 hour-long windows
    fprintf('maxWindows: %d; Nsmooth: %d\n',maxWindows,Nsmooth);
    
    %%
    nanI = ~isfinite(dayData);
    dayData(nanI) = 0;
    if Nsmooth > 1
        box = ones(Nsmooth,1)/Nsmooth;
        dayData = flipud(convn(dayData',box));
        dayData = flipud(dayData(Nsmooth:end,:));
        dayData = downsample(dayData,Nsmooth/4)';
        t = downsample(t,Nsmooth/4);
        dayData = normalizeWaveforms(dayData);
    end
    
    %%
    
end
dayData = normalizeWaveforms(dayData);

%%
lags = lags(mI);
% if sum(SNCL1 == SNCL2) < 4
%     goodColumns = (sum(nanI) == 0)';
%     cutI = abs(lags) < 2;
%     dayDataCut = dayData(cutI,goodColumns);
%     raw_shifts = zeros(size(dayData,2),1);
%     [~,~,~,~,raw_shifts_2] = apply_vdcc(dayDataCut);
%     raw_shifts(goodColumns) = raw_shifts_2;
%     dayData = apply_shifts(dayData,raw_shifts);
% end

%%
[caus,acaus,symStack] = decomposeNoiseCorrMatrix(dayData);


