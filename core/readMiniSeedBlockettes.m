function [encoding,wordOrder,dataRecordLength,blocketteType] ...
    = readMiniSeedBlockettes(raw,OffsetFirstBlockette,isBigEndian)
%
% [encoding,wordOrder,dataRecordLength] ...
% = readMiniSeedBlockettes(raw,blocketteBeginOffset,isBigEndian)
%
% reads first blockette in stream to get encoding format
% assumes only one encoding format used in entire file
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019

%%
OffsetFirstBlockette_ = OffsetFirstBlockette(1);
blocketteType = typecastArray(raw(OffsetFirstBlockette_+1:OffsetFirstBlockette_+2,:),"int16");
if isBigEndian
    blocketteType = swapbytes(blocketteType);
end

%%
blocketteType1 = blocketteType(1);
if blocketteType1 ~= 1000
    fprintf("Cannot read blockette 1000. Actual blockette returned: %d. Trying STEIM2 encoding...\n",blocketteType1);
    sizeB = size(blocketteType);
    encoding = 11*ones(sizeB);                  % 11 is Steim2 code
    dataRecordLength = 2^09*ones(sizeB);        % default record length
    wordOrder = ones(sizeB);                    % big endian as default
    return;
end

%% Data only blockette
encoding = typecastArray(raw(OffsetFirstBlockette_+5,:)',"uint8");
wordOrder = typecastArray(raw(OffsetFirstBlockette_+6,:)',"uint8");
dataRecordLength = typecastArray(raw(OffsetFirstBlockette_+7,:)',"uint8");
if ~isBigEndian
    dataRecordLength = 2.^double(dataRecordLength);
    return;
end
encoding = swapbytes(encoding);
wordOrder = swapbytes(wordOrder);
dataRecordLength = swapbytes(dataRecordLength);
dataRecordLength = 2.^double(dataRecordLength);