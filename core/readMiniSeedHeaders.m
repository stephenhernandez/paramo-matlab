function [headerInfo,OffsetFirstBlockette] = readMiniSeedHeaders(raw,isBigEndian)
%
% headerInfo = readMiniSeedHeaders(raw,isBigEndian)
%
% reads values from all headers at once
% swaps to little endian if necessary
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
sequenceNumber = strip(string(char(raw(1:6,:))'),"0");
dataQualityCode = strip(string(char(raw(7,:))'));
stationCode = deblank(strip(string(char(raw(9:13,:))')));
locationCode = deblank(strip(string(char(raw(14:15,:))')));
channelCode = deblank(strip(string(char(raw(16:18,:))')));
networkCode = deblank(strip(string(char(raw(19:20,:))')));
startTime = readBeginTime(raw(21:30,:),isBigEndian);

%%
nSamples = typecastArray(raw(31:32,:),"uint16");
sampleRateFactor = typecastArray(raw(33:34,:),"int16");
sampleRateMultiplier = typecastArray(raw(35:36,:),"int16");
nFollowingBlockettes = typecastArray(raw(40,:),"uint8");
timeCorrection = typecastArray(raw(41:44,:),"uint32");
dataBeginOffset = typecastArray(raw(45:46,:),"uint16");
OffsetFirstBlockette = typecastArray(raw(47:48,:),"uint16");

%%
if isBigEndian % if big endian, then swap bytes
    nSamples = swapbytes(nSamples);
    sampleRateFactor = swapbytes(sampleRateFactor);
    sampleRateMultiplier = swapbytes(sampleRateMultiplier);
    nFollowingBlockettes = swapbytes(nFollowingBlockettes);
    timeCorrection = swapbytes(timeCorrection);
    dataBeginOffset = swapbytes(dataBeginOffset);
    OffsetFirstBlockette = swapbytes(OffsetFirstBlockette);
end

%% vector of sample rate values for every block
if sampleRateFactor > 0
    if sampleRateMultiplier >= 0
        sampleRate = sampleRateFactor.*sampleRateMultiplier;
    else
        sampleRate = -sampleRateFactor./sampleRateMultiplier;
    end
else
    if sampleRateMultiplier >= 0
        sampleRate = -sampleRateMultiplier./sampleRateFactor;
    else
        sampleRate = 1./(sampleRateFactor.*sampleRateMultiplier);
    end
end

%%
[encoding,wordOrder,dataRecordLength,blocketteType] = ...
    readMiniSeedBlockettes(raw,OffsetFirstBlockette,isBigEndian);

%%
headerInfo = struct("stationCode",stationCode,...
    "locationCode",locationCode,...
    "channelId",channelCode,...
    "networkCode",networkCode,...
    "nSamples",nSamples,...
    "sampleRate",sampleRate,...
    "nFollowingBlockettes",nFollowingBlockettes,...
    "timeCorrection",timeCorrection,...
    "dataBeginOffset",dataBeginOffset,...
    "blockettBeginOffset",OffsetFirstBlockette,...
    "startTime",startTime,...
    "encoding",encoding,...
    "wordOrder",wordOrder,...
    "dataRecordLength",dataRecordLength,...
    "sequenceNumber",sequenceNumber,...
    "dataQualityCode",dataQualityCode,...
    "blocketteType",blocketteType);