function varargout = cutWaveforms(S,refTime,startTime,endTime,verboseFlag,nanFlag)
%
% cutWaveforms returns structure with cut waveforms
%
% Sout = cutWaveforms(Sin,refTime,startTime,endTime)
% refTime: datetime object used as reference time
% startTime: startTime (in seconds) relative to refTime
% endTime: endTime (in seconds) relative to refTime
%
% {start,end} Time can be floats or duration objects
% if Sin has multiple traces, each trace will be looped through
% individually
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

% Modified Tuesday, Jul 23, 2019
% Fixed bug so that if no data exist is requested window, an actual empty
% file gets returned

%%
if nargin < 5
    verboseFlag = false;
end

if nargin < 6
    nanFlag = false;
end

%%
varargout{1} = populateWaveforms();
if nargin < 1
    fprintf('not enough input arguments\n');
    return;
end

%% check shapes
if ~isvector(refTime)
    fprintf('input reference time(s) must be shaped like a vector\n');
    return;
end

if ~isvector(startTime)
    fprintf('input start time(s) must be shaped like a vector\n');
    return;
end

if ~isvector(endTime)
    fprintf('input end time(s) must be shaped like a vector\n');
    return;
end

%% check types
if ~isdatetime(refTime)
    fprintf('input reference time must be datetime object\n');
    return;
end

if isfloat(startTime)
    startTime = seconds(startTime);
elseif ~isduration(startTime)
    fprintf('startTime must be float or duration object\n');
    return;
end

if isfloat(endTime)
    endTime = seconds(endTime);
elseif ~isduration(endTime)
    fprintf('endTime must be float or duration object\n');
    return;
end

%% equalize shapes
S = S(:);
lS = length(S);
ltr = length(refTime);
if nargout > lS
    fprintf('requesting more outputs than available\n');
    return;
end

lts = length(startTime);
if ltr ~= lts
    startTime = repmat(startTime,ltr,1);
end

lte = length(endTime);
if ltr ~= lte
    endTime = repmat(endTime,ltr,1);
end

%%
requestedStart = refTime + startTime;
requestedEnd = refTime + endTime;
Sout = populateWaveforms([ltr,lS]);
varargout = cell(nargout,1);

%%
for i = 1:lS
    if verboseFlag
        disp(i)
    end
    S_ = S(i);
    delta = S_.delta;
    npts = S_.npts; % apparently i never use this variable...
    trueStart = S_.ref;
    data = double(pull(S_));

    %%
    gapFlagOrig = S_.gapFlag;
    if gapFlagOrig
        gapInfoOrig = S_.gapInfo;
        gapInfoOrig = sortrows(gapInfoOrig);
        gapStartOrig = gapInfoOrig(:,1);
        gapDurationOrig = gapInfoOrig(:,2);
        gapEndOrig = gapStartOrig + gapDurationOrig - 1;
        lGapOrig = length(gapStartOrig);
    else
        gapInfoOrig = [];
    end

    %%
    for j = 1:ltr
        gapFlag = gapFlagOrig;
        Sout(j,i) = S_;
        requestedStart_ = requestedStart(j);
        requestedEnd_ = requestedEnd(j);

        %%
        if isnat(trueStart)
            fprintf(2,'input data has invalid actual start time: %s.%s.%s.%s %s\n',...
                S_.kstnm,S_.knetwk,S_.khole,S_.kcmpnm,datestr(requestedStart_));
            Sout(j,i) = populateWaveforms(1);
            continue;
        end

        %%
        if requestedEnd_ < trueStart
            if verboseFlag
                fprintf('True Start: %s\n',datestr(trueStart,'yyyy-mm-dd HH:MM:SS.FFF'));
                fprintf('Requested Start: %s\n',datestr(requestedStart_,'yyyy-mm-dd HH:MM:SS.FFF'));
                fprintf('Requested End: %s\n',datestr(requestedEnd_,'yyyy-mm-dd HH:MM:SS.FFF'));
            end
            fprintf(2,'data do not exist in requested window: %s.%s.%s.%s %s\n',S_.kstnm,S_.knetwk,S_.khole,S_.kcmpnm,datestr(requestedStart_));
            Sout(j,i) = populateWaveforms(1);
            continue;
        end

        %%
        trueEnd = trueStart + S_.e;
        if verboseFlag
            fprintf('True Start: %s\n',datestr(trueStart,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf('Requested Start: %s\n',datestr(requestedStart_,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf('Requested End: %s\n',datestr(requestedEnd_,'yyyy-mm-dd HH:MM:SS.FFF'));
            fprintf('True End: %s\n',datestr(trueEnd,'yyyy-mm-dd HH:MM:SS.FFF'));
        end

        %%
        if requestedStart_ > trueEnd || requestedEnd_ < trueStart
            fprintf(2,'data do not exist in requested window: %s.%s.%s.%s %s\n',...
                S_.kstnm,S_.knetwk,S_.khole,S_.kcmpnm,datestr(requestedStart_));
            Sout(j,i) = populateWaveforms(1);
            continue; %exit current iteration and jump to next iteration
        end

        %%
        newRef = requestedStart_;
        cutMaxN = t2i(requestedEnd_,newRef,delta);
        tei = t2i(trueEnd,newRef,delta);

        if nanFlag
            dataCut = NaN(cutMaxN,1);
        else
            dataCut = zeros(cutMaxN,1);
        end

        %%
        if requestedStart_ < trueStart %prepend gap
            %rsi = 1; % i dont use but i like knowing it
            tsi = t2i(trueStart,requestedStart_,delta);
            prependedDur = tsi-1;

            gapFlag = true;
            gapInfo = [1 prependedDur];
            if gapFlagOrig
                gapStart = gapStartOrig + prependedDur;
                gapDuration = gapDurationOrig;
                gapInfo = [gapInfo; gapStart gapDuration];
            end

            if trueEnd > requestedEnd_
                % dispose of a section
                disposeNumber = tei - cutMaxN;
                partialWindow = data(1:end-disposeNumber);
                lpw = length(partialWindow);
                dataCut(tsi:tsi+lpw-1) = partialWindow;

                % remove orphaned gaps
                if gapFlagOrig
                    gapStart = gapInfo(:,1);
                    badGS = gapStart > cutMaxN;
                    if sum(badGS)
                        gapInfo(badGS,:) = [];
                    end
                end
            else
                % append gap (requestedEnd is greater than trueEnd, need to pad!)
                appendNumber = cutMaxN - tei;
                dataCut(tsi:tsi+npts-1) = data;
                gapInfo = [gapInfo; tei+1 appendNumber];
            end
        else %requestedStart_ >= trueStart
            %tsi = 1; % i dont use but i like knowing it
            rsi = t2i(requestedStart_,trueStart,delta); % is >= 1
            cutDur = rsi-1;

            if gapFlagOrig
                gapStart = gapStartOrig - cutDur;
                gapDuration = gapDurationOrig;
                gapInfo = [gapStart gapDuration];
                newGapEndIndex = gapStart + gapDuration - 1;
                deleteI = newGapEndIndex < 1;
                if sum(deleteI)
                    gapInfo(deleteI,:) = [];
                end

                if size(gapInfo,1) > 0
                    gapFlag = true;
                    gapStart = gapInfo(:,1);
                    gapDuration = gapInfo(:,2);
                    badStarts = gapStart < 1;
                    if sum(badStarts)
                        gapEndIndex = gapStart + gapDuration - 1;   %1 (order matters)
                        gapStart(badStarts) = 1;                    %2 (order matters)
                        gapDuration = gapEndIndex - gapStart + 1;
                    end
                    gapInfo = [gapStart gapDuration];
                end
            else
                gapFlag = false;
                gapInfo = [];
            end

            if trueEnd > requestedEnd_
                % dispose of a section
                dataCut = data(rsi:rsi+cutMaxN-1);

                % adjust gap info here
                % remove any orphaned gaps/gapInfo
                if gapFlag
                    gapStart = gapInfo(:,1);
                    deleteI = gapStart > cutMaxN;
                    if sum(deleteI)
                        gapInfo(deleteI,:) = [];
                    end

                    if size(gapInfo,1) < 1
                        gapFlag = false;
                        gapInfo = [];
                    end
                end
            else
                % append gap
                appendNumber = cutMaxN - tei;
                partialWindow = data(rsi:end);
                lpw = length(partialWindow);
                dataCut(1:lpw) = partialWindow;
                gapFlag = true;
                gapInfo = [gapInfo; ...
                    lpw+1 appendNumber];
            end
        end

        %%
        dataCut = dataCut(1:cutMaxN);
        Sout(j,i).d = dataCut;
        Sout(j,i).npts = cutMaxN;
        Sout(j,i).ref = newRef;
        Sout(j,i).e = seconds((cutMaxN - 1)*delta);

        %%
        [minVals,maxVals,meanVals] = minmaxmean(dataCut);
        Sout(j,i).depmin = minVals;
        Sout(j,i).depmax = maxVals;
        Sout(j,i).depmen = meanVals;

        %%
        if gapFlag
            gapStart = gapInfo(:,1);
            gapEnd = sum(gapInfo,2);

            %%
            ge = gapEnd > cutMaxN;
            if sum(ge)
                gapEnd(ge) = cutMaxN;
            end

            %%
            gs = gapStart < 1;
            if sum(gs)
                gapStart(gs) = 1;
            end

            %%
            gapDur = gapEnd - gapStart + 1;
            sumGaps = sum(gapDur);
            if any(sumGaps >= cutMaxN)
                Sout(j,i) = populateWaveforms(1);
            else
                Sout(j,i).gapFlag = gapFlag;
                Sout(j,i).gapInfo = [gapStart gapDur];
            end
        else
            Sout(j,i).gapFlag = gapFlag;
            Sout(j,i).gapInfo = [];
        end
    end

    if ltr > 1
        varargout{i} = Sout(:,i);
    end
end

%%
if ltr == 1
    varargout{1} = Sout(:);
end