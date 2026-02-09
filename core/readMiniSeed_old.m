function [S,successFlag,ENCODING,gapFlag] = readMiniSeed_old(fileName,derivedRate,verboseFlag,tStart,tEnd,interpFlag)
% clear; close all; clc;
% fileName = 'EC.SAGA..HHZ.D.2020.170';
% derivedRate = false;
% verboseFlag = false;

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

%% preallocate outputs
S = populateWaveforms();
successFlag = false;
ENCODING = [];
gapFlag = false;

%% Opening file and loading raw data
try
    s = dir(fileName);
    m = memmapfile(fileName,'Format','uint8');
catch
    if verboseFlag
        fprintf(1,'ERROR: %s\n',fileName);
    end
    return;
end

if ~s.bytes
    disp('file found, but empty; returning empty struct.');
    return;
end

%%
raw = m.Data(1:55);
testYear = typecast(raw(20:21),'uint16');
isBigEndian = testYear < 2056;

%%
try
    firstHeader = readMiniSeedHeaders(raw,isBigEndian);
catch ME
    warning(ME.message);
    return;
end

%%
dataRecordLength = firstHeader.dataRecordLength;
sizeData = size(m.Data,1);

%%
% experimental, short circuit
% if ~isnat(tStart)
%     ncols = sizeData/dataRecordLength;
%     offset = dataRecordLength*(ncols - 1);
%     lastHeader = readMiniSeedHeaders(m.Data(offset + (1:55))',isBigEndian);
%     lastStamp = double(lastHeader.startTime);
%     lastStamp = datetime(lastStamp(1),01,...
%         lastStamp(2),...
%         lastStamp(3),...
%         lastStamp(4),...
%         lastStamp(5)+lastStamp(6)*1e-4);
%     lastStamp = lastStamp + seconds(double(lastHeader.nSamples - 1)/double(lastHeader.sampleRate));
%
%     if lastStamp <= tStart
%         if verboseFlag
%             fprintf('%s %s\n',lastStamp,tStart);
%         end
%         return;
%     end
% end

%%
tStart = datenum(tStart);
tEnd = datenum(tEnd);

%%
if derivedRate && verboseFlag
    fprintf("deprecated option derivedRateFlag, not using...\n");
end

if sizeData < 1e9
    raw = m.Data;
    [S,successFlag,ENCODING,gapFlag] = cutMiniSeed(tStart,tEnd,dataRecordLength,...
        raw,isBigEndian,verboseFlag,interpFlag);
    return;
end

sizeDataOld = sizeData;
sizeData = sizeDataOld/dataRecordLength;
facts = factor(sizeData);
maxFact = max(facts);
numIt = sizeData/maxFact;

Stry = populateWaveforms(numIt);
trialIndices = maxFact*(0:numIt)';
si = trialIndices(1:end-1);
ei = trialIndices(2:end);

if verboseFlag
    fprintf('<strong>\t%s\t%s\t%s\t%s\t%s</strong>\n','i','numIt','sizeData','maxFact','sizeDataOld');
end

for i = 1:numIt
    if verboseFlag
        fprintf("%d %d %d %d %d\n",i,numIt,sizeData,maxFact,sizeDataOld);
    end
    si_ = si(i);
    ei_ = ei(i);
    raw = m.Data((si_*dataRecordLength)+1:(ei_*dataRecordLength));
    [S_,successFlag,ENCODING,gapFlag] = cutMiniSeed(tStart,tEnd,dataRecordLength,...
        raw,isBigEndian,sizeData,verboseFlag,interpFlag);
    Stry(i) = S_;
end
S = mergeWaveforms(Stry);