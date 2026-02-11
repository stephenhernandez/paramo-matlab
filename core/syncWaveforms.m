function [S,varargout] = syncWaveforms(S,verboseFlag,strictFlag,legacyWarning)
%
% syncWaveforms return structure with waveforms synced in time
%
% S = syncWaveforms(S)
% Sync traces so that they all have same begin and end
% Some traces can be shorter than their original length
%
% [S,ref,e] = syncWaveforms(S); %alternate form

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Wednesday, Jul 24, 2019
% Updated 11 JUNE 2025

%%
if nargin < 2
    verboseFlag = false;
end

if nargin < 3
    strictFlag = true;
end

if nargin < 4
    legacyWarning = false;
end

if verboseFlag && legacyWarning
    fprintf("extendFlag is deprecated. Use padWaveforms() instead.\n");
end

%%
sizeS = size(S);
S = S(:);

%%
if isscalar(S)
    if verboseFlag
        fprintf("only one trace, nothing to sync.\n");
    end
    if nargout > 1
        origins = pull(S,'ref');
        varargout{1} = origins;
        if nargout > 2
            e = pull(S,'e');
            varargout{2} = e;
        end
    end
    return;
end

%%
ref = pull(S,'ref');
goodI = find(~isnat(ref));
lgood = length(goodI);
if lgood < 2
    if verboseFlag
        fprintf("none or only 1 is good, nothing to sync.\n");
    end

    if nargout > 1
        origins = pull(S,'ref');
        varargout{1} = origins;
        if nargout > 2
            e = pull(S,'e');
            varargout{2} = e;
        end
    end
    return;
end

%%
S_ = S(goodI);
ref = ref(goodI);
b = pull(S_,'b');
e = pull(S_,'e');
deltas = pull(S_,'delta');
refs = ref+b;
etimes = refs + e;

%%
GLOBALSTARTTIME = min(refs);        % will potentially lengthen traces, no zero-padding
GLOBALENDTIME = max(etimes);        % will potentially lengthen traces, no zero-padding

%%
if GLOBALSTARTTIME >= GLOBALENDTIME
    S = populateWaveforms();
    if nargout > 1
        origins = S.ref;
        varargout{1} = origins;
        if nargout > 2
            e = S.e;
            varargout{2} = e;
        end
    end
    return;
end

%% GET THE JOB DONE
totalDur = GLOBALENDTIME-GLOBALSTARTTIME;
if ~strictFlag
    S_ = cutWaveforms(S_,GLOBALSTARTTIME,0,totalDur,verboseFlag);
    S(goodI) = S_;
    S = reshape(S,sizeS);

    %%
    if nargout > 1
        origins = pull(S,'ref');
        varargout{1} = origins;
        if nargout > 2
            e = pull(S,'e');
            varargout{2} = e;
        end
    end
    return;
end

%%
minDelta = min(deltas,[],1);
totalDur = seconds(totalDur);
nSamples = floor(totalDur/minDelta);
tMain = minDelta*(0:nSamples)'; %units in seconds

%%
for i = 1:lgood
    S__ = S_(i);
    thisref = refs(i);
    thisDelta = deltas(i);
    npts_ = S__.npts;
    localStart = seconds(thisref - GLOBALSTARTTIME);
    mainI = find(tMain >= localStart,1);
    baseVec = thisDelta*(0:npts_-1)';
    tOrig = baseVec + localStart; %with local sampling rate

    newRef = tMain(mainI);
    d1 = S__.d;
    %if thisDelta == minDelta, then baseVec == tMain
    %if thisDelta > minDelta, then baseVec different than tMain
    %in either case, tQuery is correctly defined!
    tQuery = baseVec + newRef; %with local sampling rate

    % the interp1 command has potential to chop data off fom one end, while
    % exprapolating "new" data on the far end. but at least time stamps are
    % synced between all stations...
    d2 = interp1(tOrig,d1,tQuery,"linear","extrap");
    S__ = dealHeader(S__,d2,1./thisDelta,...
        GLOBALSTARTTIME+seconds(newRef));
    S_(i) = S__;
end

%%
S(goodI) = S_;
S = reshape(S,sizeS);

%%
if nargout > 1
    origins = pull(S,'ref');
    varargout{1} = origins;
    if nargout > 2
        e = pull(S,'e');
        varargout{2} = e;
    end
end