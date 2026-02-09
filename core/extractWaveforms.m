function [S,successFlag] = extractWaveforms(varargin)
%
% extractWaveforms read SNCL(s), then cut. returns waveform structure
% associated with individual time cuts (if more than one requested).
%
% S = extractWaveforms(tStart,tEnd,"BREF","BHZ","EC","",derivedRate,verboseFlag,maxDays,filterFlag,[lfc hfc waFlag transferFlag diffFlag newFs]);
% derivedRate: false
% verboseFlag: false
% rawDataDir: default is '~/rawdata/'
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019
% Modified Saturday, Jul 27, 2019
% Modified Saturday, Aug 03, 2019
% Modified Saturday, Aug 03, 2020
% Modified Tuesday, Jan 28, 2025
% Modified Monday March 24 2025 removed option to specify directory (undoing changes implemented on jan. 28)

%
%[S,successFlag] = extractWaveforms(tStart,tEnd,stnm,chan,net,locID,derivedRate,verboseFlag,maxDays,filterFlag,[lfc hfc waFlag transferFlag diffFlag newFs]);

%%
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2015,08,13),...    % tStart
    1,...                       % tEnd
    "BREF",...                  % stnm
    "BHZ",...                   % chan
    "EC",...                    % net
    "",...                      % locID
    false,...                   % derivedRate
    false,...                   % verboseFlag
    2,...                       % maxDays
    false,...                   % filterFlag
    []};                        % filterObject [lfc hfc waFlag transferFlag diffFlag newFs]

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[tStart,tEnd,stnm,chan,net,locID,derivedRate,verboseFlag,maxDays,filterFlag,filterObject] = ...
    deal(optsToUse{:});

%%
maxDur = seconds(days(maxDays));
if ~isdatetime(tStart)
    disp('tStart must be datetime array');
    return;
end

%%
successFlag = false;
lts = length(tStart);
lte = length(tEnd);

%%
if lts == lte
    if isdatetime(tEnd)
        dur = tEnd - tStart;
    elseif isduration(tEnd)
        dur = tEnd;
    else
        dur = seconds(tEnd);
    end
elseif lte == 1
    if isfloat(tEnd)
        tEnd = seconds(tEnd);
    end
    dur = tEnd;
    dur = repmat(dur,lts,1);
else
    fprintf(2,"something went wrong\n");
    return;
end

%%
if isempty(derivedRate)
    derivedRate = false;
end

if isempty(verboseFlag)
    verboseFlag = false;
end

%%
if any(dur >= maxDur)
    fprintf(2,"some durations are too long. make smaller?\n");
    return;
end

%%
roundDaysStart = dateshift(tStart,'start','day');
roundDaysEnd = dateshift(tStart+dur,'start','day');
bigDur = days(roundDaysEnd - roundDaysStart);
loadDays = unique(roundDaysStart);
lUniqLoadDays = length(loadDays);

%%
if verboseFlag
    fprintf(1,'cutting %d window(s) from %d day(s)\n',lts,lUniqLoadDays);
end

%%
stnm = string(stnm);
chan = string(chan);
net = string(net);
locID = string(locID);

lstnm = length(unique(stnm));
lchan = length(unique(chan));
lnet = length(unique(net));
lloc = length(unique(locID));

%%
S = populateWaveforms([lts lstnm*lnet*lchan*lloc]);

%%
npoles = 4;
cumEvents = 0;
for i = 1:lUniqLoadDays
    %%
    loadDay = loadDays(i);
    dI1 = loadDay == roundDaysStart & bigDur < 1;   % window encompassed in single day
    dI2 = loadDay == roundDaysStart & bigDur >= 1;  % window straddles two days

    sumdi1 = sum(dI1);
    sumdi2 = sum(dI2);

    tic;
    if sumdi2
        [S_,nLoadedChannels] = loadWaveforms(loadDay,max(bigDur)+1,stnm,chan,net,locID,derivedRate);
        dI = dI1 | dI2;
    elseif sumdi1
        [S_,nLoadedChannels] = loadWaveforms(loadDay,1,stnm,chan,net,locID,derivedRate);
        dI = dI1;
    end

    dI = find(dI);
    if nLoadedChannels
        successFlag = true;
        if verboseFlag
            cumEvents = cumEvents + length(dI);
            fprintf("cutting: %d of %d events, day %d of %d, (%s)\n",cumEvents,lts,i,lUniqLoadDays,loadDay);
        end
    else
        if verboseFlag
            fprintf("no channels loaded for day: %s; error of some sort, fix\n",string(loadDay));
        end
        continue;
    end

    try
        if nLoadedChannels > 1
            S_ = syncWaveforms(S_,true,true,true);
        end
    catch ME
        warning(ME);
        return;
    end

    %%
    % no filtering desired
    if ~filterFlag
        for j = 1:nLoadedChannels
            Stmp = S_(j);
            Scut = cutWaveforms(Stmp,tStart(dI),0,dur(dI),verboseFlag,true);
            S(dI,j) = Scut;
        end
        continue; %goto next day
    end


    %%
    % filtering desired, execute logic
    lfo = length(filterObject);
    if ~lfo
        %filterObject is empty, set defaults
        lfc = 0.5;
        hfc = 8;
        S_ = filterWaveforms(S_,lfc,hfc,npoles);
    else
        %filterObject not empty
        if lfo < 4
            filterObject(4) = true; % set transferFlag = true
        end
        lfc = filterObject(1);
        hfc = filterObject(2);
        waFlag = filterObject(3);
        transferFlag = filterObject(4);
        lfo = length(filterObject);

        if lfo < 5
            diffFlag = false;
        else
            diffFlag = filterObject(5);
        end

        if lfo < 6
            del_ = pull(S_,'delta');
            newFs = 1./del_;
            delI = del_ < 1;
            newFs(delI) = round(1/del_(delI));
        else
            newFs = filterObject(6);
        end

        if transferFlag
            units = 'disp'; %this is internal only
            if diffFlag
                units = 'vel'; %this is internal only
            end

            if waFlag
                if diffFlag
                    units = 'acc'; %this is internal only
                else
                    units = 'vel'; %this is internal only
                end
                S_ = resampleWaveforms(S_,newFs);
            else
                tw = 2*ceil(max([1 newFs/lfc]));
                S_ = differentiateWaveforms(S_);
                S_ = resampleWaveforms(S_,newFs);
                S_ = nanGapWaveforms(S_,0);
                S_ = taperWaveforms(S_,tw);
            end

            S_ = transferWaveforms(S_,lfc,hfc,npoles,newFs,units,true,waFlag);

            if ~waFlag
                S_ = scaleWaveforms(S_,1e9);
            end
        else
            % no deconvolution!
            if diffFlag
                S_ = differentiateWaveforms(S_);
            end
            S_ = resampleWaveforms(detrendWaveforms(S_),newFs);
            S_ = filterWaveforms(S_,lfc,hfc,npoles);
        end
    end

    %%
    for j = 1:nLoadedChannels
        Stmp = S_(j);
        Scut = cutWaveforms(Stmp,tStart(dI),0,dur(dI),false,true);
        S(dI,j) = Scut;
    end
end

if ~successFlag
    if verboseFlag
        fprintf("no data found, returning empty struct\n");
    end
    S = populateWaveforms();
end
