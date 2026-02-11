function [tMain,ampMain,velMain,bazMain,medCCsMain,shiftsMain,waveformMain] = ...
    infrasoundArrayProcessing(S,varargin)

nVarargin = length(varargin);
functionDefaults = {...
    1/4,...         % lfc
    2,...           % hfc
    1e-2,...        % ampThresh
    0.9,...         % thresh
    0.6,...         % initialThresh
    1e3,...         % Fs2
    12,...          % secDur
    2/3,...         % nOverlap
    2,...           % nStrides
    1/3600,...      % sensitivity
    []};            % Fs

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[lfc,hfc,ampThresh,thresh,initialThresh,Fs2,secDur,nOverlap,...
    nStrides,sensitivity,Fs] = deal(optsToUse{:});

tic;
if isempty(Fs)
    Fs = max(round(1./pull(S,'delta')));
end
winlen = round(secDur*Fs);
stride = round(winlen*(1-nOverlap))/Fs;

%
refEllipse = referenceEllipsoid('wgs84');
elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
elementElevs = [1236; 1236; 1234; 1233; 1238];
nComps = length(elementLats);

maxN = 1e4;
badComponents = [];
shiftsMain = NaN(nComps-length(badComponents),maxN);
meanCCsMain = NaN(maxN,1);
medCCsMain = meanCCsMain;
ampMain = meanCCsMain;
ampTaperDiff = ampMain;
p2rmsMain = ampMain;
medratioMain = meanCCsMain;
stdratioMain = meanCCsMain;
velMain = meanCCsMain;
bazMain = meanCCsMain;
belMain = meanCCsMain;
tMain = NaT(maxN,1);

elementLats(badComponents) = [];
elementLons(badComponents) = [];
elementElevs(badComponents) = [];

lS = length(elementLats);
waveformMain = NaN(winlen*lS,maxN);
dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);

    %
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
    dz = [dz; elementElevs(i+1:end) - elementElevs(i)];

    %
    %disp([elementLats(i),elementLons(i)])
    %disp([d_ az])
end

Gorig = [dx dy]; % dz];
[Grows,Gcolumns] = size(Gorig);
plagsMain = NaN(Grows,maxN);
slownessMain = NaN(Gcolumns,maxN);
minNumStations = 3;

%
n = 0;
elementLats = [-2.193637; -2.193408; -2.193701; -2.193944; -2.193641];
elementLons = [-78.099766; -78.099465; -78.099223; -78.099517; -78.099506];
elementElevs = [1236; 1236; 1234; 1233; 1238];
day_ = dateshift(S(1).ref,'start','day');

refs = pull(S,'ref');
badComponents = isnat(refs);
goodComponents = ~badComponents;
lS = sum(goodComponents);
if lS < minNumStations
    fprintf(2,'not enough channels for: %s\n',datestr(day_));
    return;
end

elementLats(badComponents) = [];
elementLons(badComponents) = [];
elementElevs(badComponents) = [];

dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
end
Gorig = [dx dy];
S(badComponents) = [];

%
tw = 0.004;
if any(isfinite([lfc hfc]))
    fprintf(1,'processing: %s\n',datestr(day_));
    Sf = detrendWaveforms(...
        intWaveforms(...
        detrendWaveforms(...
        filterWaveforms(...
        taperWaveforms(...
        detrendWaveforms(...
        differentiateWaveforms(S)),tw),lfc,hfc))));

    fprintf('done filtering...\n');
else
    Sf = S;
end
Sf = syncWaveforms(detrendWaveforms(resampleWaveforms(Sf,Fs)));

%%
maxRange = 2;
sta = 5;
lta = 5;
mph = 1.5;
hFlag = false;
envFiltFlag = false;
hfc = -inf;

locsMain = [];
for i = 1:length(Sf)
    locs_ = stalta(Sf(i),sta,lta,mph,hFlag,false,envFiltFlag,hfc,false);
    locsMain = [locsMain; locs_];
    locs_ = stalta(Sf(i),sta*2,lta*2,mph*2,hFlag,false,envFiltFlag,hfc,false);
    locsMain = [locsMain; locs_];
end

%
lS = length(Sf);
t = getTimeVec(Sf);
tTrial = sort(t(locsMain));
rate = t2r(tTrial,seconds(maxRange));
fI = rate >= minNumStations;
tPot = tTrial(fI);
tPot = tPot - seconds(maxRange);

delI = false(size(tPot));
diff_tPot = seconds(diff(tPot));
fI2 = find(diff_tPot < maxRange);

delI(fI2+1) = true;
tPot(delI) = [];

nI = tPot >= Sf(1).ref & tPot+seconds(secDur) < Sf(1).ref+Sf(1).e;
nWindows = sum(nI);
if ~nWindows
    return;
end
tPot = tPot(nI);

%%
dstack = [];
% Fs = round(1./median(pull(Sf,'delta')));
% if Fs ~= newFs
%     fprintf('something went wrong with resampling step\n');
%     return;
% end

for i = 1:lS
    Scut = cutWaveforms(Sf(i),tPot,0,seconds(secDur));

    dcut = double(pull(Scut));
    %[dcut,startIndex] = cutWindows(d,winlen,nOverlap,true);
    dstack = [dstack; dcut];
end
ampNoTaper = sensitivity*0.5*peak2peak(dstack)';

tw = 0.10;
dstack = [];
%tdum = (0:winlen-1)';
taperWin = tukeywin(winlen,tw);
%taperWin(tdum > winlen/2) = 1;
for i = 1:lS
    Scut = cutWaveforms(Sf(i),tPot,0,seconds(secDur));

    dcut = double(pull(Scut));
    %[dcut,startIndex] = cutWindows(d,winlen,nOverlap,true);
    dstack = [dstack; dcut.*taperWin];
end
allAmps = sensitivity*0.5*peak2peak(dstack)';

winI = allAmps >= ampThresh;

%t_ = getTimeVec(Sf(1));
tdumcut = tPot; %t_(startIndex);
%clear t_;
nWindows = sum(winI); %size(dstack,2);
if ~nWindows
    return;
end
tPot(~winI) = [];
ampNoTaper(~winI) = [];
allAmps(~winI) = [];
dstack(:,~winI) = [];
si = 1 + winlen*(0:lS-1)';
ei = winlen*(1:lS)';

shifts = NaN(lS,nWindows);
meanCCs = NaN(nWindows,1);
medCCs = meanCCs;
medratio = meanCCs;
stdratio = meanCCs;

%
dispIndex = [];
prcWins = round(100*(1:nWindows)/nWindows);
for i = [1 10 20 30 40 50 60 70 80 90]
    dispIndex = [dispIndex; find(prcWins>i,1)];
end

%
Ndiff = 0.5*lS*(lS-1);
plags = NaN(Ndiff,nWindows);
slownii = NaN(Gcolumns,nWindows);
baz = NaN(nWindows,1);
vel = baz;
if Gcolumns > 2
    bel = baz;
end

%
fprintf("running through: %d potential windows\n",nWindows);
p2rms = NaN(size(ampNoTaper));
for i = 1:nWindows
    if ismember(i,dispIndex)
        fprintf('%g\n',floor(100*i/nWindows));
    end

    %
    G = Gorig;
    dslice = dstack(:,i);

    dslice = reshape(dslice,[winlen,lS]);
    [~,maxccp_] = apply_vdcc(dslice,[],false,false,false);
    sqmaxccp = squareform(maxccp_);
    sqmaxccp(sqmaxccp==0) = NaN;
    indivAvgCC = median(sqmaxccp,"omitnan")';
    newGood = indivAvgCC >= initialThresh;

    if sum(newGood) < minNumStations
        continue;
    end

    G1 = squareform(G(:,1));
    G2 = squareform(G(:,2));
    G1 = G1(:,newGood);
    G2 = G2(:,newGood);
    G1 = G1'; G2 = G2';
    G1 = squareform(G1(:,newGood));
    G2 = squareform(G2(:,newGood));
    G = [G1(:) G2(:)];

    dslice = detrend(dstack(:,i));
    dslice = reshape(dslice,[winlen,lS]);
    dsliceOrig = dslice;
    dslice = resample(dsliceOrig,Fs2,Fs);
    dsliceResampled = dslice(:,newGood);
    amp_ = sensitivity*0.5*median(peak2peak(dsliceResampled));
    allAmps(i) = amp_;
    p2rms(i) = median(peak2rms(dsliceResampled));

    [~,~,~,~,raw_shifts,meancc,medcc] = apply_vdcc(dsliceResampled,[],false,false,false);
    raw_shifts = raw_shifts - min(raw_shifts);
    shifts(newGood,i) = raw_shifts;
    dslice = dsliceOrig;

    meanCCs(i) = meancc;
    medCCs(i) = medcc;
    slow_ = lscov(G,(getDD(raw_shifts)/Fs2));
    slownii(:,i) = slow_;

    if Gcolumns > 2
        [azimuth,elevation,r] = cart2sph(slow_(1),slow_(2),slow_(3));
        vel(i,1) = r;
        baz(i,1) = 90 - rad2deg(azimuth);
        bel(i,1) = rad2deg(elevation);
        if bel(i) > 90
            bel(i) = bel(i) - 90;
        end
        fprintf('%g %g %g\n',vel(i),baz(i),bel(i));
    else
        baz(i,1) = 90 - atan2d(slow_(1),slow_(2)); %+360;
        vel(i,1) = 1./rssq(slow_);
        if baz(i) > 360
            baz(i) = baz(i)-360;
        end
    end

    if baz(i) > 180
        baz(i) = baz(i) - 360;
    end

    if baz(i) < -180
        baz(i) = baz(i) + 360;
    end

    % super dirty, must be a way to speed up
    medampratios = NaN(1,Ndiff);
    dslice = abs(dslice);
    refBlock = dslice;

    ri = 1;
    ln = size(refBlock,2)-1;
    for j = 1:lS-1
        d1 = dslice(:,j);
        refBlock = circshift(refBlock,-1,2);
        refBlock = refBlock(:,1:end-1);
        ln = size(refBlock,2);
        medampratios(ri:ri+ln-1) = median(log10(d1./refBlock),'omitnan');
        ri = ri + ln;
    end
    medratio(i) = median(medampratios,'omitnan');
    stdratio(i) = mad(medampratios,1,2); %'omitnan');
end

%medCCsOrig = medCCs;
goodWindowsI = medCCs >= thresh & allAmps >= ampThresh;
summi = sum(goodWindowsI);
if ~summi
    return;
end

%
%t = getTimeVec(Sf);
goodWindowsI = find(goodWindowsI);
if summi > 1
    ccGood = medCCs(goodWindowsI);
    tGood = tPot(goodWindowsI); %t(startIndex(goodWindowsI));
    tdiff = seconds(diff(tGood));
    lstride = tdiff < nStrides*stride;
    sumlstride = sum(lstride);
    while sumlstride
        lstride = find(lstride,1);
        i1 = lstride(1);
        i2 = i1 + 1;
        cc1 = ccGood(i1);
        cc2 = ccGood(i2);
        if cc1 > cc2
            goodWindowsI(i2) = [];
        else
            goodWindowsI(i1) = [];
        end
        ccGood = medCCs(goodWindowsI);
        tGood = tPot(goodWindowsI); %t(startIndex(goodWindowsI));
        tdiff = seconds(diff(tGood));
        lstride = tdiff < nStrides*stride;
        sumlstride = sum(lstride);
    end
    summi = length(goodWindowsI);
end

%
n = n + 1;
ampNoTaper = ampNoTaper(goodWindowsI);
allAmps = allAmps(goodWindowsI);
p2rms = p2rms(goodWindowsI);
dstack = dstack(:,goodWindowsI);
shifts = shifts(:,goodWindowsI);
meanCCs = meanCCs(goodWindowsI);
medCCs = medCCs(goodWindowsI);
plags = plags(:,goodWindowsI);
vel = vel(goodWindowsI);
baz = baz(goodWindowsI);
waveformMain(:,n:n+summi-1) = dstack;
if Gcolumns > 2
    bel = bel(goodWindowsI);
end

shiftsMain(:,n:n+summi-1) = shifts;
ampMain(n:n+summi-1) = allAmps;
p2rmsMain(n:n+summi-1) = p2rms;
ampTaperDiff(n:n+summi-1) = (ampNoTaper - allAmps)./allAmps;
meanCCsMain(n:n+summi-1) = meanCCs;
medCCsMain(n:n+summi-1) = medCCs;
tMain(n:n+summi-1) = tGood; %t(startIndex(goodWindowsI));
medratioMain(n:n+summi-1) = medratio(goodWindowsI);
stdratioMain(n:n+summi-1) = stdratio(goodWindowsI);
slownessMain(:,n:n+summi-1) = slownii(:,goodWindowsI);
plagsMain(:,n:n+summi-1) = plags;
bazMain(n:n+summi-1) = baz;
velMain(n:n+summi-1) = vel;
if Gcolumns > 2
    belMain(n:n+summi-1) = bel;
end

elapsedTime = toc;
fprintf("-----------------------------------------------\n");
fprintf("processed <strong>%d</strong> events on day: %s\n",summi,datestr(day_));
fprintf("Elapsed Time: <strong>%f</strong>\n",elapsedTime);
fprintf("-----------------------------------------------\n");
n = n+summi-1;

shiftsMain(:,n+1:end) = [];
ampMain(n+1:end) = [];
meanCCsMain(n+1:end) = [];
medCCsMain(n+1:end) = [];
medratioMain(n+1:end) = [];
stdratioMain(n+1:end) = [];
slownessMain(:,n+1:end) = [];
waveformMain(:,n+1:end) = [];
ampTaperDiff(n+1:end) = [];
p2rmsMain(n+1:end) = [];
tMain(n+1:end) = [];
bazMain(n+1:end) = [];
velMain(n+1:end) = [];

%%
plotFlag = false;
if plotFlag
    close all;
    nUsed = sum(isfinite(shiftsMain))';
    goodI = velMain >= 220 & velMain <= 470 & ampMain >= ampThresh & ...
        medCCsMain >= thresh & nUsed >= 3;

    figure();
    plot(tMain(goodI),velMain(goodI),'o'); zoom on; grid on; hold on;
    plot(tMain(~goodI),velMain(~goodI),'.');title("Velocity");


    figure();
    plot(tMain(goodI),medCCsMain(goodI),'o'); zoom on; grid on; hold on;
    %plot(tMain(vI),meanCCsMain(vI),'.'); zoom on; grid on;
    legend("Median","Mean");

    figure();
    bazMain(bazMain < 0) = bazMain(bazMain < 0) + 360;
    plot(tMain(goodI),bazMain(goodI),'o'); zoom on; grid on; hold on;
    plot(tMain(~goodI),bazMain(~goodI),'.'); title("BAZ");
    ax = gca;
    ax.YLim = [0 360];
    ax.YTick = (0:30:360)';

    figure();
    subplot(211);
    plot(tMain(goodI),(0:sum(goodI)-1)','o'); zoom on; grid on; title("Cumulative Number");
    subplot(212);
    nHours = 1;
    rate = t2r(tMain(goodI),hours(nHours));
    plot(tMain(goodI),rate,'.'); zoom on; grid on; title("Hourly Rate");

    figure();
    plot(tMain(goodI),ampMain(goodI),'o'); zoom on; grid on; title("Amplitudes");
    hold on;
    plot(tMain(~goodI),ampMain(~goodI),'.');
    ax = gca; ax.YScale = 'log';

    t2 = tMain(goodI);
    [~,ax] = plotWaveforms(Sf);
    ax = ax(end);
    hold(ax,'on');
    ylim_ = ax.YLim;
    for i = 1:length(t2)
        t_ = t2(i);
        plot(ax,[t_ t_],ylim_,'k','linewidth',2);
        aa = area(ax,[t_; t_+seconds(secDur)],max(ylim_)*ones(2,1),min(ylim_));
        aa.FaceAlpha = 0.5;
    end

    figure();
    ax__(1) = subplot(211);
    plot(tMain,p2rmsMain,'o'); zoom on; grid on; hold on; title("peak-to-rms");
    ax__(2) = subplot(212);
    plot(tMain,ampTaperDiff,'o'); zoom on; grid on; hold on;
    plot(tMain(goodI),ampTaperDiff(goodI),'.'); zoom on; grid on;
    linkaxes(ax__,'x');

    %figure();
    %plot(tTrial,(0:length(tTrial)-1)','.'); zoom on; grid on;
    fprintf("QC'd Events: <strong>%d</strong>, Discarded Events: <strong>%d</strong>, \n",sum(goodI),sum(~goodI));
end
