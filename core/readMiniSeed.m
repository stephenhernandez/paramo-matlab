function [S,successFlag,ENCODING,gapFlag,multiplexedFlag] = readMiniSeed(fileName,derivedRate,verboseFlag,tStart,tEnd,interpFlag)
%
% [S,successFlag,ENCODING,gapFlag] = readMiniSeed(fileName,derivedRate,verboseFlag)
%
% Code to read in miniseed file, usually called from `loadWaveforms'
%
% This code modified from version written by Martin Mityska (ReadMSEEDFast,
% Charles University), which I believe was itself originally based on code
% written by Francois Beauducel (RDMSEED,IPGP). Ive taken it upon myself to
% modify the code to suit my personal needs.
%
% S has a structure similar to SAC, see `populateWaveforms' for field names
%
% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague
%

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019
% Major revamps in June/July 2025

%% set default inputs
if nargin < 2
    derivedRate = false;
end

if nargin < 3
    verboseFlag = false;
end

if nargin < 4
    tStart = NaT;
end

if nargin < 5
    tEnd = NaT;
end

if nargin < 6
    interpFlag = false;
end

if derivedRate && verboseFlag
    fprintf("deprecated option derivedRateFlag, not using...\n");
end

%%
tStart = datenum(tStart);
tEnd = datenum(tEnd);

%% preallocate outputs
S = populateWaveforms();
successFlag = false;
ENCODING = [];
gapFlag = false;
multiplexedFlag = false;

%% Opening file and loading raw data
s = dir(fileName);
if isempty(s)
    if verboseFlag
        fprintf("File not found.\n");
    end
    return;
end

if ~s.bytes || isempty(s)
    if verboseFlag
        fprintf("File found, but empty; returning empty struct.\n");
    end
    return;
end

m = memmapfile(fileName,"Format","uint8");
raw = m.Data;

%%
MinRecordLength = 8;
offsetPositionTooBig = 48;
S = populateWaveforms();
nIter = 0;
while ~isempty(raw)
    nIter = nIter + 1;
    sizeData = size(raw,1);
    if sizeData < 21
        if verboseFlag
            fprintf("cannot do this iteration, size is %d (too small)\n",sizeData);
        end
        return;
    end

    testYear = typecast(raw(20:21),"uint16");
    isBigEndian = testYear < 2056;
    if sizeData < offsetPositionTooBig
        if verboseFlag
            fprintf("cannot do this iteration, size is %d (too small)\n",sizeData);
        end
        return;
    end

    OffsetFirstBlockette = typecastArray(raw(47:48),"uint16");
    if isBigEndian
        OffsetFirstBlockette = swapbytes(OffsetFirstBlockette);
    end

    if OffsetFirstBlockette > offsetPositionTooBig
        if verboseFlag
            fprintf("cannot do this iteration (%d), OffsetFirstBlockette is > %d. you get what you get...\n",...
                nIter,offsetPositionTooBig);
        end
        return;
    end

    dataRecordLength = typecastArray(raw(OffsetFirstBlockette+7,:)',"uint8");
    if isBigEndian
        dataRecordLength = swapbytes(dataRecordLength);
    end

    %fprintf("%d %d %d %d\n",nIter,dataRecordLength,isBigEndian,numel(raw));
    if dataRecordLength < MinRecordLength
        if verboseFlag
            fprintf("cannot do this iteration (%d), dataRecordLength is %d\n",...
                nIter,2^dataRecordLength);
        end
        return;
    end
    dataRecordLength = 2.^double(dataRecordLength);

    %%
    ncols = ceil(sizeData/dataRecordLength);
    newSize = dataRecordLength*ncols;
    if sizeData < newSize
        fprintf("expanding size of raw data..\n");
        raw_ = zeros(newSize,1,"uint8");
        raw_(1:sizeData) = raw;
        raw = raw_;
    end

    %%
    raw = reshape(raw,dataRecordLength,ncols);
    b_ = raw(9:20,:);
    goodI = all(b_ <= 127 & b_ >= 21)';
    sumBadI = sum(~goodI);

    %%
    rawb = [];
    if sumBadI
        %fprintf("some records are bad: %d %d %d\n",nIter,sumBadI,ncols);
        raw = [raw(:,goodI) raw(:,~goodI)]; %move good records to the beginning
        if sumBadI == ncols
            if nIter == 1
                %if this condition met on first iteration, then just return
                %default empty struct...
                if verboseFlag
                    fprintf(2,"this file (%s) is bad through and through. nothing returned.\n",fileName);
                end
            elseif verboseFlag
                fprintf(1,"this file (%s) has bad records, but did return some data\n",fileName);
                fprintf("\n");
            end
            return;
        end
        goodI = [true(sum(goodI),1); false(sum(~goodI),1)];
        raw = raw(:);
        firstHeader = readMiniSeedHeaders(raw,isBigEndian);
        dataRecordLength = firstHeader.dataRecordLength;
        sizeData = length(raw);
        ncols = ceil(sizeData/dataRecordLength);
        newSize = dataRecordLength*ncols;
        if sizeData < newSize
            fprintf("expanding size of raw data..\n");
            raw_ = zeros(newSize,1,"uint8");
            raw_(1:sizeData) = raw;
            raw = raw_;
        end

        raw = reshape(raw,dataRecordLength,ncols);
        OffsetFirstBlockette = typecastArray(raw(47:48),"uint16");
        if isBigEndian
            OffsetFirstBlockette = swapbytes(OffsetFirstBlockette);
        end

        allRecordLengths = typecastArray(raw(OffsetFirstBlockette+7,:)',"uint8");
        if isBigEndian
            allRecordLengths = swapbytes(allRecordLengths);
        end
        allRecordLengths = 2.^double(allRecordLengths);
        goodI = goodI & allRecordLengths == dataRecordLength;

        if ~sum(goodI) %none are good...
            if verboseFlag
                fprintf(2,"iteration: %d produced no good records.\n",nIter);
            end
            raw = [];
            continue;
        end
        rawb = raw(:,~goodI);
        raw = raw(:,goodI);
    end

    %
    if sizeData >= 1e9
        fprintf(2,"sizeData > 1e9. No code to handle this case. Add code? Leaving...\n");
        return;
    end

    [S_,successFlag,ENCODING,gapFlag,multiplexedFlag] = cutMiniSeed(tStart,tEnd,...
        dataRecordLength,raw,isBigEndian,verboseFlag,interpFlag); %typical file: 0.6 seconds

    %
    refs = pull(S_,"ref");
    if any(isnat(refs))
        if verboseFlag
            fprintf("something went wrong 1\n");
        end
        return;
    end

    %
    raw = rawb(:);
    if nIter == 1
        S = S_;
        continue;
    end

    [S,mergeSuccess] = mergeWaveforms([S; S_]);
    if ~mergeSuccess
        S = [S; S_]; %#ok<AGROW>
        multiplexedFlag = true;
        if verboseFlag
            fprintf("merge unsuccessful, aborting...");
        end
        return;
    end
end