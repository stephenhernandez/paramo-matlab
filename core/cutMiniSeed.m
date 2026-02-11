function [S,successFlag,ENCODING,gapFlag,multiplexedFlag] = cutMiniSeed(tStart,tEnd,dataRecordLength,...
    raw,isBigEndian,verboseFlag,interpFlag)
%
% [S,successFlag,ENCODING,gapFlag] = cutMiniSeed(tStart,tEnd,dataRecordLength,...
%    raw,isBigEndian,sizeData,verboseFlag,interpFlag)
%

%%
% Reading all headers of the file at once of a possibly multiplexed file...
headerInfoAllVolumes = readMiniSeedHeaders(raw,isBigEndian);
sampleRateOrig = double(headerInfoAllVolumes.sampleRate);
nSamples = double(headerInfoAllVolumes.nSamples);
stationCode = headerInfoAllVolumes.stationCode;
networkCode = headerInfoAllVolumes.networkCode;
channelCode = headerInfoAllVolumes.channelId;
locationCode = headerInfoAllVolumes.locationCode;
encodingFormat = headerInfoAllVolumes.encoding;
dataBeginOffset = headerInfoAllVolumes.dataBeginOffset;
startTimesOrig = headerInfoAllVolumes.startTime;
[startTimesOrig,sortI] = sort(startTimesOrig);

if ~issorted(sortI)
    % if not sorted, sort them
    encodingFormat = encodingFormat(sortI);
    dataBeginOffset = dataBeginOffset(sortI);
    sampleRateOrig = sampleRateOrig(sortI);
    nSamples = nSamples(sortI);
    stationCode = stationCode(sortI);
    networkCode = networkCode(sortI);
    channelCode = channelCode(sortI);
    locationCode = locationCode(sortI);
    raw = raw(:,sortI);
end

%%
dataStart = startTimesOrig(1);
dataEnd = startTimesOrig(end) + nSamples(end)/sampleRateOrig(end)/86400;
if tEnd < dataStart || tStart > dataEnd
    fprintf(2,'Error: Requested data does not span available data\n');
    successFlag = false;
    return;
end

%datetime([tStart tEnd],"ConvertFrom","datenum")
badI = false(size(sampleRateOrig)); %start with none being bad...
fints = ~isnan([tStart tEnd]);
filterCondition = sum(fints);
if filterCondition == 1
    if fints(1)
        tEndChunks = startTimesOrig + (nSamples-1)./sampleRateOrig./86400;
        lowIndex = find(tEndChunks < tStart,1,'last');
        badI(1:lowIndex) = true;
    else
        highIndex = find(startTimesOrig > tEnd,1,'first');
        badI(highIndex:end) = true;
    end
elseif filterCondition == 2
    tEndChunks = startTimesOrig + (nSamples-1)./sampleRateOrig./86400;
    lowIndex = find(tEndChunks < tStart,1,"last");
    highIndex = find(startTimesOrig > tEnd,1,"first");
    badI(1:lowIndex) = true;
    badI(highIndex:end) = true;
end

%% cut! (file is still possibly muliplexed)
sampleRateOrig(badI) = [];
nSamples(badI) = [];
stationCode(badI) = [];
networkCode(badI) = [];
channelCode(badI) = [];
locationCode(badI) = [];
startTimesOrig(badI) = [];
raw(:,badI) = [];

%%
nRecords = size(raw,2);
allSNCLs = strcat(networkCode,stationCode,locationCode,channelCode);
strLength = strlength(allSNCLs);
badI = strLength == 0;
sumbad = sum(badI);
if sumbad
    sampleRateOrig(badI) = [];
    nSamples(badI) = [];
    stationCode(badI) = [];
    networkCode(badI) = [];
    channelCode(badI) = [];
    locationCode(badI) = [];
    startTimesOrig(badI) = [];
    raw(:,badI) = [];
    allSNCLs(badI) = [];
end

if sumbad == nRecords
    fprintf(2,'Error: Data doesnt have valid metadata!\n');
    successFlag = false;
    return;
end

uniqSNCLs = unique(allSNCLs);
lUniqSNCLs = length(uniqSNCLs);
multiplexedFlag = lUniqSNCLs > 1;
if ~multiplexedFlag
    if verboseFlag
        fprintf("file is not multiplexed\n");
    end

    networkCode = networkCode(1);
    stationCode = stationCode(1);
    locationCode = locationCode(1);
    channelCode = channelCode(1);
    encodingFormat = encodingFormat(1);         % assume theyre all the same for each record in each sncl...
    dataBeginOffset = dataBeginOffset(1);       % assume theyre all the same for each record in each sncl...
    [S,successFlag,gapFlag] = readmseed_single_sncl(networkCode,...
        stationCode,locationCode,channelCode,startTimesOrig,encodingFormat,...
        raw,isBigEndian,verboseFlag,dataRecordLength,...
        dataBeginOffset,nSamples,sampleRateOrig,interpFlag);
    ENCODING = encodingFormat;
    return;
end

%%
if verboseFlag
    fprintf("this file is multiplexed, %d channels identified...\n",lUniqSNCLs);
end
S = populateWaveforms(lUniqSNCLs);
n = 0;
successFlag = false;
gapFlag = false(lUniqSNCLs,1);
ENCODING = gapFlag;
for i = 1:lUniqSNCLs
    thisSNCL = uniqSNCLs(i);
    snclI = allSNCLs == thisSNCL;
    snclI = find(snclI);

    networkCode_ = networkCode(snclI(1));
    stationCode_ = stationCode(snclI(1));
    locationCode_ = locationCode(snclI(1));
    channelCode_ = channelCode(snclI(1));
    raw_ = raw(:,snclI);
    startTimesOrig_ = startTimesOrig(snclI);
    encodingFormat_ = encodingFormat(snclI);
    dataBeginOffset_ = dataBeginOffset(snclI);
    encodingFormat_ = encodingFormat_(1);
    dataBeginOffset_ = dataBeginOffset_(1);
    nSamples_ = nSamples(snclI);
    sampleRateOrig_ = sampleRateOrig(snclI);

    [S_,successFlag_,gapFlag_] = readmseed_single_sncl(networkCode_,... %typical file: 0.5 seconds
        stationCode_,locationCode_,channelCode_,startTimesOrig_,...
        encodingFormat_,raw_,isBigEndian,verboseFlag,dataRecordLength,...
        dataBeginOffset_,nSamples_,sampleRateOrig_,interpFlag);

    if verboseFlag
        fprintf("success: %d, encoding: %d, gapFlag: %d\n",...
            successFlag_,encodingFormat_,gapFlag_);
    end

    if successFlag_
        %success with this particular sncl...
        n = n+1;
        S(n) = S_(1);
        ENCODING(n) = encodingFormat_;
        gapFlag(n) = gapFlag_;
    end
end

%%
if n < 1
    return;
end
S = S(1:n);
ENCODING = ENCODING(1:n);
gapFlag = gapFlag(1:n);
successFlag = true;