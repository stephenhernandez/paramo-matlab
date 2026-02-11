function [stack,t,S,d1,d2] = getCF(varargin)
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2018,06,01),...    % dayStart
    ["SN02","HHZ","9D",""],... 	% SNCL1
    ["SN06","HHZ","9D",""],... 	% SNCL2; true,...                    % allFlag
    32,...                      % newFs
    0,...                       % dW
    1,...                       % lfc
    4,...                       % hfc
    false,...                   % twFlag
    32768,...                   % totN
    0,...                       % nOverlap
    false,...                   % normFlag
    84,...                      % maxWindows
    65535,...                   % stackN
    true,...                    % getStack
    '~/rawdata/'};              % pathToWaveformServerOrDirectory

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[dayStart,SNCL1,SNCL2,newFs,dW,lfc,hfc,tw,totN,nOverlap,normFlag,...
    maxWindows,stackN,getStack,pathToWaveformServerOrDirectory] = deal(optsToUse{:});

%% pre-allocate
if getStack
    stack = NaN(stackN,1);
    t = NaT(1);
else
    stack = NaN(stackN,maxWindows);
    t = NaT(maxWindows,1);
end

%%
verboseFlag = false;
stnm1 = char(SNCL1(1)); chan1 = SNCL1(2); net1 = SNCL1(3); locid1 = SNCL1(4);
stnm2 = char(SNCL2(1)); chan2 = SNCL2(2); net2 = SNCL2(3); locid2 = SNCL2(4);

%%
rawDataDir = pathToWaveformServerOrDirectory;

if strcmp(stnm1(1:2),'SN')
    rawDataDir1 = '~/data/iguana/BROADBAND/';
else
    rawDataDir1 = rawDataDir;
end

if strcmp(stnm2(1:2),'SN')
    rawDataDir2 = '~/data/iguana/BROADBAND/';
else
    rawDataDir2 = rawDataDir;
end
stnm1 = string(stnm1);
stnm2 = string(stnm2);
autoCorrFlag = sum(SNCL1 == SNCL2) == 4;

%%
if verboseFlag
    fprintf('%s (%d)\n',datestr(dayStart,'yyyy-mm-dd'),day(dayStart,'doy'));
end

%%
S = loadWaveforms(dayStart,1,stnm1,chan1,net1,locid1,false,false,rawDataDir1);
startOne = S(1).ref;
endOne = startOne + S(1).e;
nptsOne = S(1).npts;

if isnat(startOne) || nptsOne <= stackN
    if verboseFlag
        fprintf('<strong>no data or not enough data</strong>\n');
        fprintf('\n');
    end
    d1 = []; d2 = [];
    S = populateWaveforms();
    return;
end

if verboseFlag
    fprintf('Trace One: %s, %s, %d\n',datestr(startOne),datestr(endOne),nptsOne);
end

%%
transferFlag = false; %temporary on 02feb2021 to test long period signal
diffFlag = false;
npoles = 4;
zeroPhaseFlag = false;
if isempty(tw)
    tw = false;
end

%%
if autoCorrFlag
    % the second is the same as the first
    % we've determined that S has data of sufficient length, preprocess...
    
    [d1,t_] = nxcorrPreprocess(S,lfc,hfc,newFs,npoles,zeroPhaseFlag,tw,diffFlag,transferFlag);
    [d1,~,endIndex,badFlag] = cutWindows(d1,totN,nOverlap,false);
    
    if badFlag
        if verboseFlag
            fprintf('<strong>no data or not enough data</strong>\n');
            fprintf('\n');
        end
        
        d1 = []; d2 = [];
        S = populateWaveforms();
        return;
    end
    d1 = doAutoCorrFreqDom(d1,false);
else
    %% the two traces are different...
    S(2,1) = loadWaveforms(dayStart,1,stnm2,chan2,net2,locid2,false,false,rawDataDir2);
    startTwo = S(2).ref;
    endTwo = startTwo + S(2).e;
    nptsTwo = S(2).npts;
    
    %
    if isnat(startTwo) || nptsTwo <= stackN
        if verboseFlag
            fprintf('<strong>no data or not enough data</strong>\n');
            fprintf('\n');
        end
        
        d1 = []; d2 = [];
        S = populateWaveforms();
        return;
    end
    
    % if verbose, print summary of both traces
    if verboseFlag
        fprintf('Trace Two: %s, %s, %d\n',datestr(startTwo),datestr(endTwo),nptsTwo);
    end
    
    % check if they overlap
    if endTwo < startOne || startTwo > endOne
        if verboseFlag
            fprintf('<strong>data do not overlap, skipping</strong>\n');
            fprintf('\n');
        end
        
        d1 = []; d2 = [];
        S = populateWaveforms();
        return;
    end
    
    %%
    [d,t_] = nxcorrPreprocess(S,lfc,hfc,newFs,npoles,zeroPhaseFlag,tw,diffFlag,transferFlag);
    [d1,~,endIndex,badFlag] = cutWindows(d(:,1),totN,nOverlap,false);
    
    %%
    if badFlag
        if verboseFlag
            fprintf('<strong>no data or not enough data</strong>\n');
            fprintf('\n');
        end
        
        d1 = []; d2 = [];
        S = populateWaveforms();
        return;
    end
    
    %%
    d2 = cutWindows(d(:,2),totN,nOverlap,false);
    d1 = fdWhiten(d1,lfc,hfc,dW,newFs);
    d2 = fdWhiten(d2,lfc,hfc,dW,newFs);
    
    %%
    if normFlag
        d1 = normalizeWaveforms(d1);
        d2 = normalizeWaveforms(d2);
    end
    d1 = doCrossCorrFreqDom(d1,d2,false);
end

%%
d1 = normalizeWaveforms(d1);

%%
if getStack
    stack = pws(d1);
    %stack = nanmean(d1,2);
else
    c = size(d1,2);
    if c > maxWindows
        c = maxWindows;
    end
    stack(:,1:c) = d1(:,1:c);
    t(1:c) = t_(endIndex(1:c));
end

%%
if ~exist('d2','var')
    d2 = d1;
end
end
