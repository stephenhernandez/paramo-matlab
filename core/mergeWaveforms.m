function [Sout,successFlag] = mergeWaveforms(S,interpFlag,verboseFlag)
%
% mergeWaveforms merge (in time) a sequence of waveform structures
% properly handles gaps and overlaps
%
% [Sout,successFlag] = mergeWaveforms(S,interpFlag,verboseFlag)
% interpFlag: linear interpolcation between gaps (defdault = false = fill with zeros)
% verboseFlag: be wordy [default = false]
% successFlag: was i successful or not?
%

% Written by Stephen Hernandez
% Instituto Geofisico, Quito, Ecuador
% Tuesday, Jul 23, 2019

%%
Sout = populateWaveforms();
successFlag = false;

%%
if nargin < 1
    fprintf("not enough input arguments\n");
    return;
end

if nargin < 2
    interpFlag = false;
end

if nargin < 3
    verboseFlag = false;
end

%%
successFlag = true;
kstnm = pull(S,'kstnm');
kcmpnm = pull(S,'kcmpnm');
khole = pull(S,'khole');
knetwk = pull(S,'knetwk');
allSNCLs = strcat(knetwk,kstnm,khole,kcmpnm);
uniqSNCLs = unique(allSNCLs);
lUniqSNCLs = length(uniqSNCLs);
Sout = populateWaveforms(lUniqSNCLs);
n = 0;
for i = 1:lUniqSNCLs
    uniqSNCLs_ = uniqSNCLs(i);
    lia = ismember(allSNCLs,uniqSNCLs_);
    Sout_ = S(lia);
    if sum(lia) > 1
        [Sout_,successFlag] = merge_(Sout_,interpFlag,verboseFlag);
    end

    if ~successFlag
        continue;
    end
    n = n+1;
    Sout(n) = Sout_;
end
Sout = Sout(1:n);
end

function [Sout,successFlag] = merge_(S,interpFlag,verboseFlag)
Sout = populateWaveforms();
successFlag = false;
if nargin < 1
    fprintf('not enough input arguments');
    return;
end

if nargin < 2
    interpFlag = false;
end

if nargin < 3
    verboseFlag = false;
end

%%
S = S(:);
ref = pull(S,'ref');
badI = isnat(ref);
ref(badI) = [];
S(badI) = [];

%%
lS = length(S);
Sout = S;
successFlag = false;
if lS == 0
    return;
elseif lS == 1
    successFlag = true;
    return;
end

%%
sigDigs = 3;
[ref,sI] = sort(ref);
S = S(sI);

%%
e = pull(S,'e');
e = ref + e;

%%
delta = pull(S,'delta');
Fs = 1./delta;
if Fs(1) < 1 && Fs(1) > 0
    delta = round((delta)*(10^sigDigs))/(10^sigDigs);
elseif Fs(1) >= 1
    delta = 1./(round((Fs)*(10^sigDigs))/(10^sigDigs));
else
    fprint(2,'something has gone wrong\n');
    return;
end

%%
kstnm = pull(S,'kstnm');
kcmpnm = pull(S,'kcmpnm');
khole = pull(S,'khole');
knetwk = pull(S,'knetwk');

%%
badI = ~ismember(kstnm(2:end),kstnm(1));
if any(badI) %check kstnm
    fprintf(2,'Error merging: not all have same sensor name\n');
    badI = find(badI)+1;
    lbad = length(badI);
    for ll = 1:lbad
        fprintf('wanted: %s, got: %s\n',kstnm(1),kstnm(badI(ll)));
    end
    return
end

%%
badI = ~ismember(kcmpnm(2:end),kcmpnm(1));
if any(badI) %check kcmpnm
    fprintf(2,'Error merging: not all have same component name\n');
    badI = find(badI)+1;
    lbad = length(badI);
    for ll = 1:lbad
        fprintf('wanted: %s.%s, got: %s.%s\n',kcmpnm(1),kstnm(1),kcmpnm(badI(ll)),kstnm(badI(ll)));
    end
    return
end

%%
badI = ~ismember(khole(2:end),khole(1));
if any(badI) %check khole
    fprintf(2,'Error merging: not all have same loc ID\n');
    badI = find(badI)+1;
    lbad = length(badI);
    for ll = 1:lbad
        fprintf('wanted: %s.%s.%s, got: %s.%s.%s\n',kstnm(1),khole(1),kcmpnm(1),...
            kstnm(badI(ll)),khole(badI(ll)),kcmpnm(badI(ll)));
    end
    return
end

%%
badI = ~ismember(knetwk(2:end),knetwk(1));
if any(badI) %check network
    fprintf(2,'Error merging: not all have same knetwk\n');
    badI = find(badI)+1;
    lbad = length(badI);
    for ll = 1:lbad
        fprintf('wanted: %s.%s.%s.%s, got: %s.%s.%s.%s\n',knetwk(1),khole(1),kstnm(1),kcmpnm(1),...
            knetwk(badI(ll)),khole(badI(ll)),kstnm(badI(ll)),kcmpnm(badI(ll)));
    end
    return
end

%%
if any(diff(delta)) %check sampling rates
    disp(delta);
    fprintf(2,'Error merging files from SNCL %s.%s.%s.%s: files must have same sampling rate\n',knetwk(1),khole(1),kstnm(1),kcmpnm(1));
    ME = MException('mergeWaveforms:SampleRateError','Error merging: files must have same sampling rate');
    throw(ME)
end

%% made it this far, all is well, set delta/Fs to a scalar
delta = delta(1);
Fs = 1/delta;
if delta <= 1
    Fs = round(Fs);
end

%% remove redundant files
lengthOfNegatives = 1;
while lengthOfNegatives
    trueI = true(lS,1);
    for i = 2:lS
        redundancyCheck = round(Fs*seconds(e(i) - e(i-1)));
        if redundancyCheck <= 0
            %we've identified a redundant file
            if verboseFlag
                disp(['redundant file found: ',num2str(i)]);
            end
            trueI(i) = false;
            continue; % short-circuit loop so that i can remove file...
        end
    end
    % remove redundant file
    lS = sum(trueI);
    e = e(trueI);
    S = S(trueI);

    % check for redundant files again
    ediff = seconds(diff(e));
    lengthOfNegatives = sum(ediff <= 0);
end

%%
data = cell(lS,1);
for i = 1:lS
    data{i} = S(i).d;
end
ref = pull(S,'ref');
npts = pull(S,'npts');
b = pull(S,'b');
stla = pull(S,'stla');
stlo = pull(S,'stlo');
evla = pull(S,'evla');
evlo = pull(S,'evlo');
evdp = pull(S,'evdp');
dist = pull(S,'dist');
az = pull(S,'az');
baz = pull(S,'baz');
gcarc = pull(S,'gcarc');
cmpaz = pull(S,'cmpaz');
cmpinc = pull(S,'cmpinc');
kstnm = pull(S,'kstnm');
khole = pull(S,'khole');
kcmpnm = pull(S,'kcmpnm');
knetwk = pull(S,'knetwk');
gapFlags = pull(S,'gapFlag');
gapFlagMain = any(gapFlags);

%% populate first trace
i = 1;
si = 1;
main = NaN(sum(npts),1); %this is just an initial estimate, its possible the final trace will be longer (gaps) or shorter (overlaps)
nCurrentTrace = npts(i);
currentEndIndex = e(i);

%disp([length(data{i}) nCurrentTrace])
main(si:si+nCurrentTrace-1) = data{i};

si = si + nCurrentTrace;    % start index of next block
ei = si - 1;                % end index of this block
gapStartIndex = si;         % potential start index of gap if skipped samples found (else ignored)

%% merge the rest of the traces
for i = 2:lS
    gapFlagNow = gapFlags(i);
    if gapFlagNow
        gapInfoNow = S(i).gapInfo;
        gapInfoNow(:,1) = ei + gapInfoNow(:,1);
        S(i).gapInfo = gapInfoNow;
    end
    dNext = double(data{i});

    %%
    prevEnd = currentEndIndex;
    nCurrentTrace = npts(i);
    currentEndIndex = e(i);
    timeGap = seconds(ref(i) - prevEnd);    % duration in seconds between successive files
    skippedLength = round(Fs*timeGap) - 1; 	% number of skipped (or overlapped) samples between successive files

    %% a gap exists, either fill with constant or interpolate (linear)
    if skippedLength > 0
        if ~gapFlagMain
            gapFlagMain = true;
        end

        %%
        gapFlagNow = gapFlags(i);
        if gapFlagNow
            gapInfoNow = S(i).gapInfo;
            gapInfoNow(:,1) = skippedLength + gapInfoNow(:,1);
            S(i).gapInfo = gapInfoNow;
        end

        %%
        if interpFlag
            if verboseFlag
                disp('gap exists, interpolating.');
            end
            X = [1; 2+skippedLength];
            V = [main(ei); dNext(1)];
            Xq = (2:2+skippedLength-1)';
            Vq = interp1(X,V,Xq);
        else
            if verboseFlag
                disp('gap exists, filling gap with NaNs.');
            end
            Vq = NaN(skippedLength,1);
        end

        %%
        gapInfo_ = [gapStartIndex skippedLength];   % index where gap starts and how long it lasts
        dNext = [Vq; dNext];                        %#ok<AGROW> % prepend this block with gap-value or interpolated values
        main(si:si+nCurrentTrace+skippedLength-1,1) = dNext;

        gapInfo = S(i-1).gapInfo;
        gapInfo = [gapInfo; gapInfo_];              %#ok<AGROW>
        S(i-1).gapInfo = gapInfo;

        si = si + nCurrentTrace + skippedLength;    % start index of next block
        ei = si - 1;                                % end index of this block
    elseif skippedLength < 0
        %%
        gapFlagNow = gapFlags(i);
        if gapFlagNow
            gapInfoNow = S(i).gapInfo;
            gapInfoNow(:,1) = skippedLength + gapInfoNow(:,1);
            S(i).gapInfo = gapInfoNow;
        end

        %% files overlap in time, averaging them
        nCurrentTrace = nCurrentTrace + skippedLength;
        overlapLength = abs(skippedLength);
        if verboseFlag
            fprintf("overlap exists, taking average.\n");
        end
        overlapOrig = main(ei-overlapLength+1:ei);
        overlapNext = dNext(1:overlapLength);
        overlappedData = mean([overlapOrig overlapNext],2,"omitnan");

        main(ei-overlapLength+1:ei) = overlappedData;
        dNext = dNext(1+overlapLength:end);         % we assume length(dNext) is greater than 1
        main(si:si+nCurrentTrace-1,1) = dNext;

        si = si + nCurrentTrace;                	% start index of next block
        ei = si - 1;                                % end index of this block
    else
        %% neither gaps nor overlaps present (seamless)
        if verboseFlag
            fprintf("neither gaps nor overlaps present (seamless)\n");
        end
        main(si:si+nCurrentTrace-1,1) = dNext;
        si = si + nCurrentTrace;                	% start index of next block
        ei = si - 1;                                % end index of this block
    end
    gapStartIndex = si;                             % potential start index of future gap if skipped samples found (else ignored)
end

%%
main = main(1:ei);
Sout = populateWaveforms();
Sout.ref = ref(1);
Sout.b = b(1);
Sout.kstnm = kstnm(1);
Sout.khole = khole(1);
Sout.kcmpnm = kcmpnm(1);
Sout.knetwk = knetwk(1);
Sout.stla = stla(1);
Sout.stlo = stlo(1);
Sout.evla = evla(1);
Sout.evlo = evlo(1);
Sout.evdp = evdp(1);
Sout.dist = dist(1);
Sout.az = az(1);
Sout.baz = baz(1);
Sout.gcarc = gcarc(1);
Sout.cmpaz = cmpaz(1);
Sout.cmpinc = cmpinc(1);

%%
lnew = ei;
dfin = isfinite(main);
sumnan = sum(~dfin);
if ~sumnan
    if verboseFlag
        fprintf('no NaNs detected...\n');
    end
else
    gapStart = find(diff(dfin)<0) + 1;
    gapEnd = find(diff(dfin)>0);
    lgs = size(gapStart,1);
    lge = size(gapEnd,1);
    if lgs > lge || ~dfin(end)
        gapEnd = [gapEnd; lnew];
        lge = lge + 1;
    end

    if lge > lgs || ~dfin(1)
        gapStart = [1; gapStart];
        lgs = lgs + 1;
    end

    %%
    lE = min([lgs lge]);
    gapDur = gapEnd(1:lE) - gapStart(1:lE) + 1;
    Sout.gapInfo = [gapStart(1:lE) gapDur];
    Sout.gapFlag = true;
end

%% assign final fields
Sout.npts = ei;
Sout.delta = delta;
Sout.e = seconds((Sout.npts - 1)*Sout.delta);
Sout.d = main;

%%
[minVals,maxVals,meanVals] = minmaxmean(main);
Sout.depmin = minVals;
Sout.depmax = maxVals;
Sout.depmen = meanVals;
successFlag = true;
end