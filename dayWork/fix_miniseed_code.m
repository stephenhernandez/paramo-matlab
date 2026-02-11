function [S,successFlag,ENCODING,gapFlag] = readMiniSeed(fileName,derivedRate,verboseFlag,tStart,tEnd,interpFlag)
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

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019
% major revamps in June 2025

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

%% preallocate outputs
S = populateWaveforms();
successFlag = false;
ENCODING = [];
gapFlag = false;

%% Opening file and loading raw data
s = dir(fileName);
if ~s.bytes
    disp('file found, but empty; returning empty struct.');
    return;
end

m = memmapfile(fileName,'Format','uint8');
raw = m.Data(1:56);
testYear = typecast(raw(21:22),'uint16');
testDay = typecast(raw(23:24),'uint16');
isBigEndian = testDay >= 512 || (ismember(testDay,[1,256]) && (testYear > 2311 || testYear < 1799));
raw = m.Data;

%%
S = [];
while ~isempty(raw)
    firstHeader = readMiniSeedHeaders(raw,isBigEndian);
    dataRecordLength = firstHeader.dataRecordLength;
    sizeData = length(raw);
    n_records = floor(sizeData/dataRecordLength);
    raw = reshape(raw,dataRecordLength,n_records);

    b_ = raw(9:20,:);
    badI = any(b_ > 127 | b_ < 21);

    headerInfoAllVolumes = readMiniSeedHeaders(raw,isBigEndian);
    allRecordLengths = headerInfoAllVolumes.dataRecordLength;
    goodRecordLength = allRecordLengths(~badI);
    uniqRecordLengths = unique(goodRecordLength);
    dataRecordLength_ = uniqRecordLengths(1);
    goodI = allRecordLengths == dataRecordLength_;
    raw_ = raw(:,goodI);

    if sizeData < 1e9
        [S_,successFlag,gapFlag] = cutMiniSeed(tStart,tEnd,dataRecordLength_,...
            raw_,isBigEndian,verboseFlag,interpFlag);
    end
    refs = pull(S_,"ref");
    if any(isnat(refs))
        fprintf("something went wrong 1\n");
        break;
    end
    S = [S; S_];

    raw(:,goodI) = [];
    raw = raw(:);
end