function [tabs,snr] = eventDetector(varargin)
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2015,08,01),...    % tStart
    datetime(2015,08,31),...    % tEnd
    2,...                       % mph
    5,...                       % short term average [sec.]
    10,...                      % long term average [sec.]
    1,...                       % lfc
    8,...                       % hfc
    "VCH1",...                  % stnm
    "HHZ",...                   % chan
    "EC",...                    % network
    "",...                      % locID
    false,...                   % diffFlag
    4,...                       % npoles
    '~/rawdata/'};              % raw data directory

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,mph,sta,lta,lfc,hfc,stnm,chan,net,locID,diffFlag,npoles,rawDataDir] = deal(optsToUse{:});

%%
if isfloat(tEnd)
    days = tStart:tStart+tEnd-1;
else
    days = tStart:tEnd;
end
ldays = length(days);

if isempty(mph)
    mph = 2;
end

if isempty(sta)
    sta = 5;
end

if isempty(lta)
    lta = 10;
end

if isempty(lfc)
    lfc = -inf;
end

if isempty(hfc)
    hfc = -inf;
end
if isempty(npoles)
    npoles = 4;
end

%%
tabs = [];
snr = [];
cornersfin = isfinite([lfc hfc]);

%%
for i = 1:ldays
    disp(days(i));
    S_ = loadWaveforms(days(i),1,stnm,chan,net,locID,0,0,rawDataDir);
    if isnat(S_.ref)
        disp([datestr(days(i)),' is empty. skipping.'])
    else
        if diffFlag
            S_ = differentiateWaveforms(S_);
        end
        
        if any(cornersfin)
            S_ = filterWaveforms(S_,lfc,hfc,npoles,[],false,false);
        end
        
        % [locs_,snr_] = stalta(S_,sta,lta,mph,true,0);
        hFlag = false;
        plotFlag = false;
        envFiltFlag = false;
        [locs_,snr_] = stalta(S_,sta,lta,mph,hFlag,plotFlag,envFiltFlag);
        
        ll = length(locs_);
        
        if ll
            t_ = getTimeVec(S_);
            tabs = [tabs; t_(locs_)];
            snr = [snr; snr_];
        end
    end
end
