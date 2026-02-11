function [trace,gapInfo,gapFlag,nOverlaps,delta,dr] = ...
    decodeMiniSeedSignal(record_start_times,encodedSignalMatrix,...
    encodingFormat,isBigEndian,dataRecordLength,origSampleRate,...
    dataBeginOffset,interpFlag,velFlag,averageFlag)
%
% [trace,gapInfo,gapFlag,delta,dr,acc,nGoodPerColumn] = decodeMiniSeedSignal(tStarts,...
%    encodedSignalMatrix,encodingFormat,isBigEndian,BLOCK_RECORD_SIZE,...
%    origSampleRate,nSamples)
%
% Only encoding formats 4, 5, 10, and 11 have been tested
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019
% Major mods: June 2025

%% set defaults
if nargin < 9
    velFlag = true;
end

if nargin < 9
    averageFlag = true;
end

trace = [];
gapInfo = [];
gapFlag = false;
delta = 1/origSampleRate;
dr = [];
nOverlaps = [];

%% Decoding routine for STEIM2 only
switch encodingFormat
    case 0      % --- decoding format: ASCII text
        trace = typecastArray(encodedSignalMatrix,"int8");
        return;
    case 1      % --- decoding format: 16-bit integers
        trace = typecastArray(encodedSignalMatrix,"int16");
        if isBigEndian
            trace = swapbytes(trace);
        end
        return
    case 2      % --- decoding format: 24-bit integers
        % unsupported
        fprintf("This option not supportedn");
        return;
    case 3      % --- decoding format: 32-bit integers
        trace = typecastArray(encodedSignalMatrix,"int32");
        if isBigEndian
            trace = swapbytes(trace);
        end
        return;
    case 4      % --- decoding format: IEEE floating point
        trace = typecastArray(encodedSignalMatrix,"single"); %this is a column vector
        if isBigEndian
            trace = swapbytes(trace);
        end
        return;
    case 5      % --- decoding format: IEEE double precision floating point
        trace = typecastArray(encodedSignalMatrix,"double");
        if isBigEndian
            trace = swapbytes(trace);
        end
        return;
    case {10,11,19}
        steim = find(encodingFormat == [10,11,19]);
        n_bytes = 4;
        cut_record_length = size(encodedSignalMatrix,1);
        number_chunks_of_n_bytes_per_record = cut_record_length/n_bytes;
        signal_matrix = typecastArray(encodedSignalMatrix,"uint32");

        if isBigEndian
            signal_matrix = swapbytes(signal_matrix);
        end

        % each original record of length 448 bytes gets reduced to 112 chunks of 4 bytes (32 bits)...
        signal_matrix = reshape(signal_matrix,number_chunks_of_n_bytes_per_record,[]); % each record now has 112 "floats"
        numelSignalMatrix = numel(signal_matrix);

        % read first int from every 64 int chunk (contains encoded nibbles)
        Q = signal_matrix(1:16:number_chunks_of_n_bytes_per_record,:); % keep only 7 "floats"? why discard > 90%?
        Qcol = numel(Q);

        % reshape to one long row vector
        Q = reshape(Q,1,[]);

        % prepare matrix Q for bitshift
        Q = repmat(Q,16,1); %dup data 16 times? why?

        % prepare bit mask for bitshift
        bshiftMask = repmat(-30:2:0,Qcol,1)';
        bQ = bitshift(Q,bshiftMask); %the original 7 "special floats" have each been shifted 16 times to give 16 new variations

        % decode nibbles
        nibbles = bitand(bQ,oldbitcmp32(0,2));

        % forward integration constant
        x0 = bitsign(signal_matrix(2,:),32); % returns signed double value from unsigned n-bit number x

        nibbles = nibbles(:);
        signal_matrix = signal_matrix(:);
        switch steim
            case 1
                % STEIM-1: 3 cases following the nibbles
                maxValuesInChunk = 4;
                acc = NaN(maxValuesInChunk,numelSignalMatrix);	% initiates array with NaN

                k = nibbles == 1;			% nibble = 1 : four 8-bit differences
                acc(1:4,k) = bitsplit(signal_matrix(k),32,8);
                if ~isBigEndian
                    acc(1:4,k)=flipud(acc(1:4,k));
                end

                k = nibbles == 2;			% nibble = 2 : two 16-bit differences
                acc(1:2,k) = bitsplit(signal_matrix(k),32,16);
                if ~isBigEndian
                    acc(1:4,k)=flipud(acc(1:4,k));
                end

                k = nibbles == 3;			% nibble = 3 : one 32-bit difference
                acc(1,k) = bitsign(signal_matrix(k),32);
            case 2
                % STEIM-2: 7 cases following the nibbles and dnib
                maxValuesInChunk = 7;
                acc = NaN(maxValuesInChunk,numelSignalMatrix);	% initiates array with NaN
                dnib = bitshift(signal_matrix,-30);

                kk = nibbles == 1;
                acc(1:4,kk) = bitsplit(signal_matrix(kk),32,8);
                if ~isBigEndian %isLittleEndian
                    acc(1:4,kk) = flipud(acc(1:4,kk));
                end

                k2 = nibbles == 2;
                k3 = nibbles == 3;
                dnib0 = dnib == 0;
                dnib1 = dnib == 1;
                dnib2 = dnib == 2;
                dnib3 = dnib == 3;

                kk = k2 & dnib1;   % dnib = 1 : one 30-bit difference
                acc(1,kk) = bitsign(signal_matrix(kk),30);

                kk = k2 & dnib2;   % dnib = 2 : two 15-bit differences
                acc(1:2,kk) = bitsplit(signal_matrix(kk),30,15);

                kk = k2 & dnib3;   % dnib = 3 : three 10-bit differences
                acc(1:3,kk) = bitsplit(signal_matrix(kk),30,10);

                kk = k3 & dnib0;   % dnib = 0 : five 6-bit difference
                acc(1:5,kk) = bitsplit(signal_matrix(kk),30,6);

                kk = k3 & dnib1;   % dnib = 1 : six 5-bit differences
                acc(1:6,kk) = bitsplit(signal_matrix(kk),30,5);

                kk = k3 & dnib2;   % dnib = 2 : seven 4-bit differences (28 bits!)
                acc(1:7,kk) = bitsplit(signal_matrix(kk),28,4);
            case 3
                % unsupported
                fprintf("This option not supported\n");
                return;
        end
    otherwise
        errorMsg = sprintf("Error: unknown data encoding.\n Supported formats are 4, 5, 10, or 11.\n");
        fprintf(2,errorMsg);
        return;
end

%% Building matrix of column vectors of signal from data blocks
acc = acc(:);
lacc = size(acc,1);
expanded_record_length = maxValuesInChunk*((dataRecordLength-dataBeginOffset)/n_bytes); % the parentheses matter, do not delete!!!!
mod_lddd_rowsize = mod(lacc,expanded_record_length);
if mod_lddd_rowsize %if not eqal to zero...
    acc = [acc; NaN(expanded_record_length-mod_lddd_rowsize,1)];
end
acc = reshape(acc,expanded_record_length,[]);
%disp(size(acc));

%%
t_i_s = 86400*(record_start_times - record_start_times(1)); %start times in seconds
%t_i_s(1:20)
bmI = isfinite(acc);
nGoodPerColumn = sum(bmI)';
t_i_e = t_i_s + delta*(nGoodPerColumn-1);
seams = round((t_i_s(2:end) - t_i_e(1:end-1))/delta)-1;
diffts = diff(t_i_s);
sps = nGoodPerColumn(1:end-1)./diffts;
imperfections = seams ~= 0;
nImperfections = sum(imperfections);
%figure(); plot(t_i_s,t_i_e,'.'); zoom on; grid on;
%figure(); plot(t_i_s,nGoodPerColumn,'.'); zoom on; grid on;

if ~nImperfections %if no gaps/overlaps
    %best case scenario: no imperfections, no overlaps, no gaps
    if velFlag
        firstI = find(bmI,1);
        acc(firstI(1),:) = x0; %potentially dangerous, assume all columns start at same index...
        acc = cumsum(acc,"omitnan"); %integrate
    end
    trace = acc(bmI); %<-- this step takes long but i cant speed up
    delta = 1./median(sps);
    return;
end

%%
if velFlag
    firstI = find(bmI,1);
    acc(firstI(1),:) = x0; %potentially dangerous, assume all columns start at same index...
    acc = cumsum(acc,1,"omitnan"); %integrate
end

%% delete completely overlapped windows (special case, if they exist)
overlapI = seams < 0;
nOverlaps = sum(overlapI);
% figure(); plot(t_i_s(1:end-1),seams,'.'); zoom on; grid on;
% table(t_i_s(1:8),seams(1:8),nGoodPerColumn(1:8))
% fprintf("%d %d\n",length(t_i_s),length(unique(t_i_s)));

if nOverlaps
    uniq_t_i_s = unique(t_i_s);
    gc_t_i_s = groupcounts(t_i_s);
    multi_t = uniq_t_i_s(gc_t_i_s > 1);
    lia = ismember(t_i_s,multi_t);
    delI = lia & [seams; 0] < 0 & -[seams; 0] == nGoodPerColumn;
    if sum(delI)
        %fprintf("there are repeated records...\n");
        record_start_times(delI) = [];  % delete row (100% necessary)
        acc(:,delI) = [];               % delete row (100% necessary)
        bmI(:,delI) = [];               % delete row (100% necessary)
        t_i_s = 86400*(record_start_times - record_start_times(1)); %start times in seconds
        nGoodPerColumn = sum(bmI)';
        t_i_e = t_i_s + delta*(nGoodPerColumn-1);
        seams = round((t_i_s(2:end) - t_i_e(1:end-1))/delta)-1;
        diffts = diff(t_i_s);
        sps = nGoodPerColumn(1:end-1)./diffts;
    end
end

overlapI = seams < 0;
nOverlaps = sum(overlapI);
if nOverlaps
    delI = [-seams; 0] >= nGoodPerColumn | [0; -seams] >= nGoodPerColumn;
    sumdel = sum(delI);
    while sumdel
        %fprintf("there are completely overlapped windows...\n");
        % %find(delI)
        %disp([sumdel length(delI)])

        delI = find(delI);
        record_start_times(delI) = [];  % delete row (100% necessary)
        acc(:,delI) = [];               % delete row (100% necessary)
        bmI(:,delI) = [];               % delete row (100% necessary)
        t_i_s = 86400*(record_start_times - record_start_times(1)); %start times in seconds
        nGoodPerColumn = sum(bmI)';
        t_i_e = t_i_s + delta*(nGoodPerColumn-1);
        seams = round((t_i_s(2:end) - t_i_e(1:end-1))/delta)-1;
        diffts = diff(t_i_s);
        sps = nGoodPerColumn(1:end-1)./diffts;
        delI = [-seams;0] >= nGoodPerColumn | [0;-seams] >= nGoodPerColumn;
        sumdel = sum(delI);
    end
end

% fix partially overlapping windows (if they exist)
imperfections = seams ~= 0;
nImperfections = sum(imperfections);
overlapI = seams < 0;
nOverlaps = sum(overlapI);

if nOverlaps
    %fprintf("there are partially overlapped windows...\n");
    overlapI = find(overlapI);
    for i = 1:nOverlaps
        currentI = overlapI(i);
        overlapLength = -seams(currentI);
        nGoodCurrent_ = nGoodPerColumn(currentI);

        if overlapLength >= nGoodCurrent_
            fprintf("something seriously wrong, this case should have been handled already.\n");
            return;
        end

        bmI_orig = bmI(:,currentI);
        bmI_current = find(bmI_orig);
        dCurrent = acc(bmI_current,currentI);
        if averageFlag
            dOverlapped = dCurrent(end-overlapLength+1:end);
            nextI = currentI + 1;
            bmInext = bmI(:,nextI);
            dNext = acc(bmInext,nextI);
            dNext_ = dNext(1:overlapLength);
            overlappedData = mean([dOverlapped dNext_],2,"omitnan");
            dNext(1:overlapLength) = overlappedData;
            acc(bmInext,nextI) = dNext;
        end

        dCurrent(end-overlapLength+1:end) = [];
        bmI_orig = false(size(bmI_orig));
        bmI_current(end-overlapLength+1:end) = [];
        bmI_orig(bmI_current) = true;
        acc(bmI_orig,currentI) = dCurrent;
        bmI(:,currentI) = bmI_orig;
    end
    %fprintf("processed %d overlaps\n",nOverlaps);
end

trace = acc(bmI); %<-- this step takes long but i cant speed up
%delta = 1./median(sps(~imperfections))
if nOverlaps == nImperfections
    return;
end

%% inject gaps (if they exist)
nGoodPerColumn = sum(bmI)';
t_i_e = t_i_s + delta*(nGoodPerColumn-1);
seams = round((t_i_s(2:end) - t_i_e(1:end-1))/delta)-1;
seams = [0; seams];
ei = cumsum(nGoodPerColumn)+cumsum(seams);
si = [1; ei(1:end-1)+1];

gapI = seams > 0;
nGaps = sum(gapI);
gapFlag = true;
gapInfo = NaN(nGaps,2);
gapI = find(gapI);
traceTmp = trace; %<-- has no gaps, havent injected them yet
gei = seams;
gei(gapI) = si(gapI) + seams(gapI) - 1;
tei = cumsum(nGoodPerColumn);
trace = zeros(sum(nGoodPerColumn)+sum(seams),1);
mainStartI = 1;
tmpStartI = 1;
for i = 1:nGaps
    iNow = gapI(i);
    iPrev = iNow - 1;

    gapStartI = si(iNow);
    tmpEndI = tei(iPrev);
    gapEndI = gei(iNow);
    mainEndI = ei(iPrev);
    skippedLength = seams(iNow);
    gapInfo(i,:) = [gapStartI skippedLength];

    %fprintf("skipped length: %d, column index: %d\n",skippedLength,iNow);
    gap = NaN(skippedLength,1);
    if interpFlag
        V = [traceTmp(tmpEndI); traceTmp(gapStartI)];
        X = [0; 1+skippedLength];
        Xq = (1:skippedLength)';
        gap = interp1(X,V,Xq);
    end
    trace(mainStartI:mainEndI) = traceTmp(tmpStartI:tmpEndI);
    trace(gapStartI:gapEndI) = gap;
    mainStartI = gapEndI+1;
    tmpStartI = tmpEndI+1;
end
trace(mainStartI:end) = traceTmp(tmpStartI:end);
%fprintf("processed %d gaps\n",nGaps);