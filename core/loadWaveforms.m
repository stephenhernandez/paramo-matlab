function [W,n] = loadWaveforms(varargin)
%
% loadWaveforms returns structure with requested waveform data
%
% S = loadWaveforms(datetime(2015,08,01),1,"BREF","BHZ","EC","",derivedRate,verboseFlag,rawDataDir);
% derivedRate: false
% verboseFlag: false
% rawDataDir: defaults are ["~/rawdata/";"~/rawdata_cotopaxi/";"~/rawdata_old/"]
%
% returns at least one empty struct should no data be found
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2015,08,13),...
    2,...
    "BREF",...
    "BHZ",...
    "EC",...
    "",...
    false,...
    false,...
    ["~/rawdata/";"~/rawdata_prod/";"~/rawdata_old/"]};

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,stnm,chan,net,locID,derivedRate,verboseFlag,rawDataDirs] = deal(optsToUse{:});

%%
if isfloat(tStart)
    tStart = dn2dt(tStart);
end

if isfloat(tEnd)
    dayList = tStart:tStart+tEnd-1;
else
    dayList = tStart:tEnd;
end

%%
lstnm = length(unique(stnm));
lchan = length(unique(chan));
lnet = length(unique(net));
lloc = length(unique(locID));
totalCombinations = lstnm*lchan*lnet*lloc;

%%
W = populateWaveforms(totalCombinations);
PWD = pwd; %always tack on current directory to list of folders to look through...
rawDataDirs = [rawDataDirs; string(PWD)];
lRawDirs = length(rawDataDirs);

%% get the job done
n = 0;
for i = 1:lstnm
    stnm_ = char(stnm(i));
    for j = 1:lchan
        chan_ = char(chan(j));
        for k = 1:lnet
            net_ = char(net(k));
            for l = 1:lloc
                locID_ = char(locID(l));

                successFlag = false;
                m = 0;
                while ~successFlag & m < lRawDirs
                    m = m+1;
                    raw_dir = rawDataDirs(m);
                    [S_,nSuccess] = loadSingleSNCL(dayList,stnm_,chan_,net_,locID_,derivedRate,verboseFlag,raw_dir);

                    %%
                    if ~nSuccess
                        if verboseFlag
                            fprintf(1,'Could not import: %s.%s.%s.%s for day beginning: %s within dir: %s\n',...
                                stnm_,net_,locID_,chan_,string(dayList(1)),raw_dir);
                        end
                        continue;
                    else
                        successFlag = true;
                    end

                    %%
                    S_ = S_(1:nSuccess);
                    if nSuccess > 1
                        try
                            [S_,successFlag] = mergeWaveforms(S_);
                        catch
                            successFlag = false;
                        end
                    end

                    %%
                    if ~successFlag || isnat(S_.ref)
                        % unsuccessful, say so
                        fprintf(2,"error: %s.%s.%s.%s, %s\n",net_,stnm_,locID_,chan_,dayList(1));
                        continue;
                    end
                    n = n + 1;
                    W(n) = S_;
                end
            end
        end
    end
end

%% return at least one empty struct should no data be found
if ~n
    W = W(1);
    return;
end
W = W(1:n);
