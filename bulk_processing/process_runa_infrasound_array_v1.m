function [vI,nUsed,tMain,ampMain,bazMain,velMain,meanCCsMain,medCCsMain] = ...
    process_runa_infrasound_array_v1(S)
warning off signal:findpeaks:largeMinPeakHeight

%
thresh = 0.75;
initialThresh = 0.6;
newFs = 100;
Fs2 = 2e3;
secDur = 6;                                 % 12 seconds
winlen = round(secDur*newFs);
noiseDur = 3;
lfc = 1;
hfc = lfc*4;
sensitivity = 1/3600;
ampThresh = 1e-2;

%
refEllipse = referenceEllipsoid('wgs84');

% -78.4670855	-0.5934007 %p1
% -78.4669994	-0.5949173 %p2
% -78.4669167	-0.5947333 %pv

%pablo version, station order: 1,2,3,4,5
elementLats = [-0.5943264; -0.5944882; -0.5943; -0.5940137; -0.5942238];
elementLons = [-78.4673559; -78.4670238; -78.4670167; -78.4670507; -78.4667395];
elementElevs = [1236; 1236; 1234; 1233; 1238];

%orig version, station order: 2,3,4,1,5 (likely wrong)
% elementLats = [-0.594564; -0.5943; -0.5941; -0.5943; -0.594316];
% elementLons = [-78.4670; -78.4667; -78.4671; -78.4674; -78.46706];
% elementElevs = [1236; 1236; 1234; 1233; 1238];

maxN = 1e4;
meanCCsMain = NaN(maxN,1);
medCCsMain = meanCCsMain;
ampMain = meanCCsMain;
velMain = meanCCsMain;
bazMain = meanCCsMain;
tMain = NaT(maxN,1);

lS = length(elementLats);
dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
    dz = [dz; elementElevs(i+1:end) - elementElevs(i)];
end
Gorig = [dx dy];
[~,Gcolumns] = size(Gorig);

minNumStations = 3;
n = 0;
refs = pull(S,'ref');
badComponents = isnat(refs);
goodComponents = ~badComponents;
lS = sum(goodComponents);
if lS < minNumStations
    fprintf(2,"not enough data, returning nothing\n");
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end

shiftsMain = NaN(lS,maxN);
S(badComponents) = [];
elementLats(badComponents) = [];
elementLons(badComponents) = [];
elementElevs(badComponents) = [];

dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);
    stel = elementElevs(i);
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
    dz = [dz; -stel + elementElevs(i+1:end)];
end

tw = 0.0008;
Sf = syncWaveforms(S,0,1,0);
Sf = detrendWaveforms(...
    intWaveforms(...
    taperWaveforms(...
    filterWaveforms(...
    taperWaveforms(...
    syncWaveforms(...
    detrendWaveforms(...
    differentiateWaveforms(Sf))),tw),lfc,hfc),tw)));

%%
Sf = resampleWaveforms(Sf,newFs);
tunif = getTimeVec(Sf);
if isempty(tunif)
    fprintf(2,"Timevector is invalid, Sf is probably empty\n");
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end

allLocs = [];
for i = 1:lS
    locs1 = stalta(Sf(i),noiseDur,noiseDur,2,false,false,false,10,true);
    allLocs = [allLocs; locs1];
end

allLocs = sort(allLocs);
t1 = tunif(allLocs);
rate = t2r(t1,seconds(2));
i2 = rate >= minNumStations; %<-- experiment with removeDuplicateMatches here!!!
t1 = t1(i2);
rate = rate(i2);

t1 = removeRepeatedMatches(t1,rate,noiseDur/2);
if isempty(t1)
    fprintf('something went wrong with t1\n');
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end
clear rate i2;

Fs = round(1./median(pull(Sf,'delta')));
if Fs ~= newFs
    fprintf('something went wrong with resampling step\n');
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end

cutEndTimes = t1-seconds(noiseDur/2)+seconds(secDur);
badTimesI = cutEndTimes >= min(pull(Sf,'ref') + pull(Sf,'e'));
t1(badTimesI) = [];
if isempty(t1)
    fprintf('t1 is empty, nothing to do\n');
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end

dstack = [];
allAmps = [];
for i = 1:lS
    S1cut = cutWaveforms(Sf(i),t1-seconds(noiseDur/2),0,seconds(secDur),false);
    try
        d1cut = pull(S1cut);
        d1cut = d1cut(1:winlen,:);
    catch
        lsf = length(Sf); clear ax;
        for kk = 1:lsf
            figure(1); hold on;
            ax(kk) = subplot(lsf,1,kk);
            plot(getTimeVec(Sf(kk)),Sf(kk).d);
        end
        linkaxes(ax,'x');

        disp(winlen)
        t1-seconds(noiseDur/2)
        seconds(secDur)
        pull(Sf,'ref')
        pull(Sf,'ref') + pull(Sf,'e')
        fprintf("i: %d; lS: %d; len(t1): %d\n",i,lS,length(t1));
        t1(i)
        Sf(i).gapInfo
        return;
    end
    allAmps = [allAmps; sensitivity*0.5*peak2peak(d1cut)];
    dstack = [dstack; d1cut];
end

allAmps = median(allAmps)';
tdumcut = pull(S1cut,'ref');
winI = allAmps >= ampThresh & ~isnat(tdumcut);
startIndex = t2i(tdumcut,Sf(1).ref,1/Fs);

%%
clear tdumcut;
nWindows = sum(winI);
if ~nWindows
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    fprintf("no events detected with STALTA method\n");
    return;
end

%%
allAmps(~winI) = [];
dstack(:,~winI) = [];
startIndex(~winI) = [];

t1 = tunif(startIndex);
allAmpsOrig = allAmps;
[t1,allAmps,removeIndices] = removeRepeatedMatches(t1,allAmpsOrig,secDur/2);
startIndex = t2i(t1,Sf(1).ref,1/Fs);
clear tunif;

for i = 1:length(removeIndices)
    rI = removeIndices{i};
    dstack(:,rI) = [];
end

nWindows = length(startIndex);
shifts = NaN(lS,nWindows); % lS can be either 3, 4, or 5 (if minNumStations = 3)
meanCCs = NaN(nWindows,1);
medCCs = meanCCs;
Ndiff = 0.5*lS*(lS-1);
plags = NaN(Ndiff,nWindows);
slownii = NaN(Gcolumns,nWindows);
baz = meanCCs;
vel = meanCCs;

for i = 1:nWindows
    dslice = dstack(:,i);
    ampTmp = allAmps(i);
    if ampTmp < ampThresh
        %fprintf("amp too small, skipping...\n");
        continue;
    end

    dslice = reshape(dslice,[winlen,lS]);
    [~,maxccp_] = apply_vdcc(dslice,[],false,false,false);
    sqmaxccp = squareform(maxccp_);
    sqmaxccp(sqmaxccp==0) = NaN;
    indivAvgCC = median(sqmaxccp,"omitnan")';

    goodStationsI = indivAvgCC >= initialThresh;
    nGoodStations = sum(goodStationsI);
    if nGoodStations < minNumStations
        %fprintf(2,"Window %d: Skipping, not enough good data\n",i);
        continue;
    end

    [d_,az_] = distance(mean(elementLats(goodStationsI)),mean(elementLons(goodStationsI)),...
        elementLats(goodStationsI),elementLons(goodStationsI),refEllipse);
    G = [d_.*cosd(90-az_),d_.*sind(90-az_)];
    G = [getDD(G(:,1)) getDD(G(:,2))];

    dslice = detrend(dstack(:,i));
    dslice = reshape(dslice,[winlen,lS]);
    dsliceOrig = dslice;
    dslice = resample(dsliceOrig,Fs2,newFs);
    dsliceResampled = dslice(:,goodStationsI);
    amp_ = sensitivity*0.5*median(peak2peak(dsliceResampled));
    allAmps(i) = amp_;

    Gshift = Gvdcc(nGoodStations);

    [maxccp,deltaT1] = doccFreqCircShift(dsliceResampled,false);
    W = maxccp.*speye(nGoodStations*(nGoodStations-1)*0.5);
    medcc = median(maxccp,"all","omitnan");
    meancc = mean(maxccp,"all","omitnan");
    raw_shifts = Gshift \ [deltaT1; 0];
    raw_shifts = raw_shifts / Fs2;
    dt = getDD(raw_shifts);
    shifts(goodStationsI,i) = raw_shifts;

    meanCCs(i) = meancc;
    medCCs(i) = medcc;
    slow_ = ((G'*W*G)^(-1))*G'*W*dt;
    slownii(:,i) = slow_;

    baz(i,1) = 90-atan2d(slow_(2),slow_(1)) + 180;
    vel(i,1) = 1./rssq(slow_);
    if baz(i) > 360
        baz(i) = baz(i)-360;
    end
    plags(1:length(dt),i) = dt;
end

goodTimeWindowsI = medCCs >= thresh & allAmps >= ampThresh;
summi = sum(goodTimeWindowsI);
if ~summi
    fprintf("no windows match medCC and/or amplitude threshold\n");
    vI = [];
    nUsed = [];
    tMain = [];
    ampMain = [];
    bazMain = [];
    velMain = [];
    meanCCsMain = [];
    medCCsMain = [];
    return;
end

%
t = getTimeVec(Sf);
goodTimeWindowsI = find(goodTimeWindowsI);

n = n + 1;
allAmps = allAmps(goodTimeWindowsI);
shifts = shifts(:,goodTimeWindowsI);
meanCCs = meanCCs(goodTimeWindowsI);
medCCs = medCCs(goodTimeWindowsI);
vel = vel(goodTimeWindowsI);
baz = baz(goodTimeWindowsI);

shiftsMain(:,n:n+summi-1) = shifts;
ampMain(n:n+summi-1) = allAmps; %sensitivity*0.5*peak2peak(dstack)';
meanCCsMain(n:n+summi-1) = meanCCs;
medCCsMain(n:n+summi-1) = medCCs;
tMain(n:n+summi-1) = t(startIndex(goodTimeWindowsI));
bazMain(n:n+summi-1) = baz;
velMain(n:n+summi-1) = vel;

n = n+summi-1;
shiftsMain(:,n+1:end) = [];
ampMain(n+1:end) = [];
meanCCsMain(n+1:end) = [];
medCCsMain(n+1:end) = [];
tMain(n+1:end) = [];
bazMain(n+1:end) = [];
velMain(n+1:end) = [];

nUsed = sum(isfinite(shiftsMain))';
vI = velMain >= 240 & velMain <= 460 & ampMain >= 0.04 & ...
    medCCsMain >= thresh & nUsed >= minNumStations; % & bazMain >= 270 & bazMain <= 350;
