function [trace,gapInfo,tStarts,gapFlag,delta,dr,blockMatrix,numelColumn] = ...
    decodeMiniSeedSignal(tStarts,encodedSignalMatrix,encodingFormat,...
    isBigEndian,dataRecordLength,origSampleRate,nSamples,interpFlag)

%
% [trace,gapInfo,tStarts,gapFlag,delta] = decodeMiniSeedSignal(tStarts,...
%    encodedSignalMatrix,encodingFormat,isBigEndian,BLOCK_RECORD_SIZE,...
%    origSampleRate,nSamples)
%
% Only encoding formats 4,10, and 11 have been tested
%

% Original code by Martin Mityska (2014)
% Faculty of Science
% Charles University in Prague

% Modified by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Monday, Jul 22, 2019
% More mods: Friday, Mar 21, 2025

%% set defaults
trace = [];
gapInfo = [];
gapFlag = false;
delta = 1/origSampleRate;

%% Decoding routine for STEIM2 only
switch encodingFormat
    %         Todo add support for next encoding formats
    case 0      % --- decoding format: ASCII text
        signalMatrix=typecastArray((reshape(encodedSignalMatrix',1,[])')','int8');
        trace = signalMatrix;

    case 1      % --- decoding format: 16-bit integers
        signalMatrix=typecastArray((reshape(encodedSignalMatrix',2,[])')','int16');
        if isBigEndian
            signalMatrix=swapbytes(signalMatrix(:));
        end
        trace = signalMatrix;

    case 2      % --- decoding format: 24-bit integers
        % unsupported
        disp('This option not supported');

    case 3      % --- decoding format: 32-bit integers
        signalMatrix = typecastArray((reshape(encodedSignalMatrix',4,[])')','int32');
        if isBigEndian
            signalMatrix = swapbytes(signalMatrix(:));
        end
        trace = signalMatrix;

    case 4      % --- decoding format: IEEE floating point
        %D.EncodingFormatName = 'FLOAT32';
        %retype file to float32
        signalMatrix = typecastArray((reshape(encodedSignalMatrix',4,[])')','single'); %this is a column vector
        if isBigEndian
            signalMatrix = swapbytes(signalMatrix(:));
        end
        totSamples = sum(nSamples);
        trace = double(signalMatrix(1:totSamples));

    case 5      % --- decoding format: IEEE double precision floating point
        signalMatrix = typecastArray((reshape(encodedSignalMatrix',8,[])')','double');
        if isBigEndian
            signalMatrix = swapbytes(signalMatrix(:));
        end
        trace = signalMatrix;

    case {10,11,19}
        steim = find(encodingFormat == [10,11,19]);
        %fprintf('steim encoding format: %d\n',steim);
        signalMatrix = typecastArray((reshape(encodedSignalMatrix',4,[])),'uint32');
        if isBigEndian
            signalMatrix = swapbytes(signalMatrix);
        end
        signalMatrix = reshape(signalMatrix,numel(encodedSignalMatrix(1,:))/4,[]);
        [nr,nc] = size(signalMatrix);
        numelSignalMatrix = nr*nc; %numel(signalMatrix);

        % read first int from every 64 int chunk (contains encoded nibbles)
        Q = signalMatrix(1:16:nr,:);

        % reshape to one long row vector
        Q = reshape(Q,1,size(Q,1)*size(Q,2));

        % prepare matrix Q for bitshift
        Q = repmat(Q,16,1);

        % prepare bit mask for bitshift
        bshiftMask = repmat(-30:2:0,size(Q,2),1)';
        bQ = bitshift(Q,bshiftMask);

        % decode nibbles
        nibbles = bitand(bQ,oldbitcmp32(0,2));

        % forward integration constant
        x0 = bitsign(signalMatrix(2,:),32);

        switch steim
            case 1
                % STEIM-1: 3 cases following the nibbles
                maxValuesInChunk = 4;
                ddd = NaN(maxValuesInChunk,numelSignalMatrix);	% initiates array with NaN

                k = find(nibbles == 1);			% nibble = 1 : four 8-bit differences
                if ~isempty(k)
                    ddd(1:4,k) = bitsplit(signalMatrix(k),32,8);
                    if ~isBigEndian
                        ddd(1:4,k)=flipud(ddd(1:4,k));
                    end
                end

                k = find(nibbles == 2);			% nibble = 2 : two 16-bit differences
                if ~isempty(k)
                    ddd(1:2,k) = bitsplit(signalMatrix(k),32,16);

                    if ~isBigEndian
                        ddd(1:4,k)=flipud(ddd(1:4,k));
                    end

                end

                k = find(nibbles == 3);			% nibble = 3 : one 32-bit difference
                if ~isempty(k)
                    ddd(1,k) = bitsign(signalMatrix(k),32);
                end
            case 2
                % STEIM-2: 7 cases following the nibbles and dnib
                maxValuesInChunk = 7;
                ddd = NaN(maxValuesInChunk,numelSignalMatrix);	% initiates array with NaN
                dnib = bitshift(signalMatrix(:),-30);
                nibblesVec = nibbles(:);

                k = nibblesVec == 1;
                if sum(k)
                    try
                        ddd(1:4,k) = bitsplit(signalMatrix(k),32,8);
                    catch
                        dr = [];
                        return;
                    end

                    if ~isBigEndian %isLittleEndian
                        ddd(1:4,k) = flipud(ddd(1:4,k));
                    end
                end

                k = nibblesVec == 2;  % nibble = 2 : must look in dnib
                if sum(k)
                    kk = k & (dnib == 1);   % dnib = 1 : one 30-bit difference
                    if sum(kk)
                        ddd(1,kk) = bitsign(signalMatrix(kk),30);
                    end
                    kk = k & (dnib == 2);   % dnib = 2 : two 15-bit differences
                    if sum(kk)
                        ddd(1:2,kk) = bitsplit(signalMatrix(kk),30,15);
                    end
                    kk = k & (dnib == 3);   % dnib = 3 : three 10-bit differences
                    if sum(kk)
                        ddd(1:3,kk) = bitsplit(signalMatrix(kk),30,10);
                    end
                end

                k = nibblesVec == 3;  % nibble = 3 : must look in dnib
                if sum(k)
                    kk = k & (dnib == 0);   % dnib = 0 : five 6-bit difference
                    if sum(kk) %~isempty(kk)
                        ddd(1:5,kk) = bitsplit(signalMatrix(kk),30,6);
                    end
                    kk = k & (dnib == 1);   % dnib = 1 : six 5-bit differences
                    if sum(kk)
                        ddd(1:6,kk) = bitsplit(signalMatrix(kk),30,5);
                    end
                    kk = k & (dnib == 2);   % dnib = 2 : seven 4-bit differences (28 bits!)
                    if sum(kk)
                        ddd(1:7,kk) = bitsplit(signalMatrix(kk),28,4);
                    end
                end
        end

        %% Building matrix of column vectors of signal from data blocks
        ddd = reshape(ddd,[],1);
        if mod(size(ddd,1),((dataRecordLength-64)/4)*maxValuesInChunk) ~= 0
            ddd = [ddd; ...
                NaN(((dataRecordLength-64)/4)*maxValuesInChunk-mod(size(ddd,1),...
                ((dataRecordLength-64)/4)*maxValuesInChunk),1)];
        end
        blockMatrix = reshape(ddd,((dataRecordLength-64)/4)*maxValuesInChunk,[]);
        cBlockMatrix = size(blockMatrix,2);

        %% sort tstarts and block matrix according to time
        [tStarts,tI] = sort(tStarts);
        blockMatrix = blockMatrix(:,tI);
        x0 = x0(tI);

        %%
        isfinddd = isfinite(ddd);
        sumisfinddd = sum(isfinddd);
        trace = zeros(sumisfinddd,1);

        %%
        si = 1;
        gapInfo = [];

        %%
        % can i modify code so that i dont have to run loop? question posed
        % on 21 march 2025

        bmI = isfinite(blockMatrix);
        numelColumn = sum(bmI)';
        %tStarts2 = datetime(tStarts,'ConvertFrom','datenum');
        currentEnd2 = tStarts + (numelColumn-1)./(origSampleRate*86400);
        timeDifference2 = 86400*(tStarts(2:end) - currentEnd2(1:end-1));
        skippedLength2 = round(origSampleRate*timeDifference2) - 1;

        % gapIndices = skippedLength2 > 0;
        % overlapIndices = skippedLength2 < 0;
        seamlessIndices = skippedLength2 == 0;
        column = blockMatrix(:,1);
        bmI_ = isfinite(column);
        column = column(bmI_);
        dNext = cumsum([x0(1); column(2:numelColumn(1))]);
        lSeamless = sum(seamlessIndices);
        %last_index = t2i(datetime(currentEnd2(end),'ConvertFrom','datenum'),tStarts2(1),delta); %use 'round' version from 01APR2025

        numelColumn_ = numelColumn(1);
        trace(si:si+numelColumn_-1,1) = dNext;
        si = si + numelColumn_;
        if lSeamless == cBlockMatrix-1
            %no gaps, no overlaps
            %fprintf("no gaps, no overlaps\n");
            dr = 86400*diff(tStarts)./(numelColumn(1:end-1));
            trace = miniseed_seamless(trace,bmI,cBlockMatrix,blockMatrix,x0,...
                numelColumn,si);
        else
            %overlaps or gaps exist
        end


        %% identify whether or not gaps have been found
        if ~isempty(gapInfo)
            gapFlag = true;
        end
        delta = median(dr);

    otherwise
        errorMsg = sprintf('Error: unknown data encoding.\n Supported formats are 4, 10, or 11.\n');
        fprintf(2,errorMsg);
        return;
end
