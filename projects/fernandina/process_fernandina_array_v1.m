function [vI,nUsed,tMain,ampMain,bazMain,velMain,meanCCsMain,medCCsMain] = ...
    process_fernandina_array_v1(S)
warning off signal:findpeaks:largeMinPeakHeight

%
staltaThresh = 4;
thresh = 0.5;
initialThresh = 0.4;
newFs = 100;
Fs2 = 5e3;
secDur = 4; %3;
winlen = round(secDur*newFs);
noiseDur = 2;
lfc = 2; %2;
hfc = 8; %lfc*4;
sensitivity = 1;
ampThresh = 1e1;

%
% kstnms = ["FARC"; "FARE"; "FARN"; "FARS"; "FARW"];
refEllipse = referenceEllipsoid('wgs84');
elementLats = [-0.300175; -0.300077; -0.298125; -0.30236; -0.300363];
elementLons = [-91.56933; -91.567347; -91.569528; -91.568687; -91.571017];
elementElevs = [50.4; 59.1; 26.5; 83.8; 58.3];

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
Sf = interpolateWaveforms(Sf);

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
    %locs1 = stalta(Sf(i),2,10,staltaThresh,false,false,false);
    %allLocs = [allLocs; locs1];

    %temporary, so that all subwindows are processed, change back when done!
    ti = getTimeVec(Sf(i));
    locs1 = (1:winlen/2:length(ti))';
    allLocs = [allLocs; locs1(2:end-1)];
    clear ti;
end

allLocs = sort(allLocs);
t1 = tunif(allLocs);
rate = t2r(t1,seconds(noiseDur));
%figure(); plot(t1,rate,'.'); zoom on;
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
    S1cut = cutWaveforms(Sf(i),t1-seconds(noiseDur),0,seconds(secDur),false);
    try
        d1cut = pull(S1cut);
        d1cut = d1cut(1:winlen,:);
        %plotWaveforms(S1cut);
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
startIndex = t2i(tdumcut,Sf(1).ref,1/Fs)
winI = allAmps >= ampThresh & isfinite(startIndex);

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
[t1,allAmps,removeIndices] = removeRepeatedMatches(t1,allAmps,secDur/4);
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
%     figure(1); subplot(nWindows,1,i);
%     plot(dslice); zoom on;
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
    if baz(i) > 180
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
vI = velMain >= 235 & velMain <= 400 & ampMain >= 0.01 & ...
    medCCsMain >= thresh & nUsed >= minNumStations & bazMain >= 300 & bazMain <= 320;
