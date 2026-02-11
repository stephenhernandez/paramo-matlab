function [S,successFlag,gapFlag] = cut_miniseed_single_sncl(startTimesOrig,... %headerInfoAllVolumes,...
    encodingFormat,raw,isBigEndian,verboseFlag,dataRecordLength,...
    dataBeginOffset,nSamples,sampleRateOrig,interpFlag)
%
%
%

%%
allSNCLs = strcat(networkCode,stationCode,locationCode,channelCode);
uniqSNCLs = unique(allSNCLs);
lUniqSNCLs = length(uniqSNCLs);
S = populateWaveforms(lUniqSNCLs);
startCol = dataBeginOffset+1;
n = 0;
for aa = 1:lUniqSNCLs
    thisSNCL = uniqSNCLs(aa);
    uniqSnclI = allSNCLs == thisSNCL;

    if ~sum(uniqSnclI)
        if verboseFlag
            fprintf(2,"empty sncl\n");
        end
        continue;
    end

    if verboseFlag
        fprintf('Found: %s\n',thisSNCL);
    end

    %%
    stnm = unique(stationCode(uniqSnclI));
    ntwk = unique(networkCode(uniqSnclI));
    chan = unique(channelCode(uniqSnclI));
    hole = unique(locationCode(uniqSnclI));

    %%
    encodedSignalMatrix = raw(uniqSnclI,startCol:end); %sliceColumns);
    nSamplesBeforeDecompression = nSamples(uniqSnclI);
    startTimes = startTimesOrig(uniqSnclI);
    sampleRateChunks = abs(sampleRateOrig(uniqSnclI)); % HACK: take abs() because sometimes a negative is returned

    if isempty(sampleRateChunks)
        if verboseFlag
            fprintf(2,'sampleRateChunks is empty\n');
        end
        return;
    end

    %%
    if any(diff(sampleRateChunks))
        if verboseFlag
            fprintf(2,'not all samplerate chunks are the same\n');
        end
        %figure(); plot(dn2dt(startTimes),sampleRateChunks,'.'); zoom on; grid on;
        return;
    end

    %%
    origSampleRate = double(sampleRateChunks(1));
    if ~origSampleRate
        if verboseFlag
            fprintf('using derived rate (temporarily)\n');
        end
        %
        origSampleRate = median(1./((diff(startTimes)*86400)./nSamplesBeforeDecompression(1:end-1)));
    end

    % get the signal
    %try
    successFlag = true;
    firstTime = startTimes(1);
    [trace,gapInfo,gapFlag] = decodeMiniSeedSignal(startTimes,...
        encodedSignalMatrix,encodingFormat,isBigEndian,dataRecordLength,...
        origSampleRate,nSamplesBeforeDecompression,interpFlag);

    if isempty(trace)
        if verboseFlag
            fprintf(2,'empty trace\n');
        end
        successFlag = false;
        return;
    end

    %% compare both rates and output some verbosity (if desired)
    % sigDigs = 5;
    % dr = round((1/delta)*(10^sigDigs))/(10^sigDigs);
    % clashFlag = ~(dr == origSampleRate);
    % if clashFlag
    %     if verboseFlag
    %         disp('beware, derived rate from time stamps and nominal rate are not the same');
    %     end
    % end

    % if derivedRate
    %     if verboseFlag
    %         if clashFlag
    %             disp(['rates clash. using derived sample rate (',num2str(dr,sigDigs),') versus nominal rate (',num2str(origSampleRate,sigDigs),').']);
    %         else
    %             disp(['rates do not clash. using derived sample rate (',num2str(dr,sigDigs),') versus nominal rate (',num2str(origSampleRate,sigDigs),').']);
    %         end
    %     end
    %     delta = 1/dr;
    % else
    %     if verboseFlag
    %         if clashFlag
    %             disp(['rates clash. using nominal sample rate (',num2str(origSampleRate,sigDigs),') versus derived rate (',num2str(dr,sigDigs),').']);
    %         else
    %             disp(['rates do not clash. using nominal sample rate (',num2str(origSampleRate,sigDigs),') versus derived rate (',num2str(dr,sigDigs),').']);
    %         end
    %     end
    delta = 1/origSampleRate;
    %end

    % populate structures
    n = n + 1;
    % if n > 1
    %     % grow S
    %     %S(n,1) = populateWaveforms();
    % end

    npts = size(trace,1);
    S(n).('kstnm') = stnm;
    S(n).('ref') = dn2dt(firstTime);
    S(n).('kcmpnm') = chan;
    S(n).('knetwk') = ntwk;
    S(n).('khole') = hole;
    S(n).('delta') = delta;
    S(n).('d') = double(trace);
    S(n).('b') = seconds(0);
    S(n).('npts') = npts;
    S(n).('e') = seconds((npts - 1)*S(n).delta); % assumes b = 0
    S(n).('gapFlag') = gapFlag;
    S(n).('gapInfo') = gapInfo;

    %%
    [minVals,maxVals,meanVals] = minmaxmean(trace);
    S(n).('depmin') = minVals;
    S(n).('depmax') = maxVals;
    S(n).('depmen') = meanVals;

    %catch ME
    %warning(ME.message)
    %fprintf('some type of error in cutMiniSeed...\n');
    %successFlag = false;
    %rethrow(ME)
    %return;
    %end
end