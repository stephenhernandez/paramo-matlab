clear; close all;
lfc = 0.25;
hfc = 1;
ampThresh = 2e-2;
thresh = 0.9;
initialThresh = 0.6;
Fs2 = 1e3;
secDur = 12;
% nOverlap = 2/3;
% nStrides = 2;
% sensitivity = 1/3600;
plotFlag = true;

dayStart = datetime(2021,10,14);
dayEnd = datetime(2022,08,31);
dInc = 1;
dayVec = (dayStart:dInc:dayEnd)';
lDays = length(dayVec);

nComps = 5;
minNumStations = 3;
maxN = 2e5;
badComponents = [];
shiftsMain = NaN(nComps-length(badComponents),maxN);
medCCsMain = NaN(maxN,1);
ampMain = medCCsMain;
velMain = medCCsMain;
bazMain = medCCsMain;
tMain = NaT(maxN,1);

kholelist = ["01";"02";"03";"04";"05"];
Fs = 50;
if isempty(Fs)
    Fs = max(round(1./pull(S,'delta')));
end
winlen = round(secDur*Fs);
waveformMain = NaN(winlen*nComps,maxN);

n = 0;
for i = 1:lDays
    day_ = dayVec(i);
    tic;
    S = populateWaveforms(nComps);
    for j = 1:nComps
        khole_ = kholelist(j);
        S_ = loadWaveforms(day_,dInc,"SAG1","BDF","EC",khole_); %this order cannot change
        S(j) = S_;
    end
    toc;

    refs = pull(S,'ref');
    badComponents = isnat(refs);
    goodComponents = ~badComponents;
    lS = sum(goodComponents);
    if lS < minNumStations
        fprintf(2,'no data for: %s\n',datestr(day_));
        continue;
    end

    Sf = resampleWaveforms(detrendWaveforms(intWaveforms(detrendWaveforms(...
        filterWaveforms(taperWaveforms(detrendWaveforms(...
        differentiateWaveforms(S)),0.0004),lfc,hfc)))),Fs);

    [tMain_,ampMain_,velMain_,bazMain_,medCCsMain_,shiftsMain_,waveformMain_] = ...
        infrasoundArrayProcessing(Sf,-inf,-inf,ampThresh,thresh,initialThresh,...
        Fs2,secDur); %,nOverlap,nStrides,sensitivity,Fs);

    sumGood = length(tMain_);
    if ~sumGood
        continue;
    end

    %%
    n = n + 1;
    nComps_ = size(shiftsMain_,1);
    newWinlen = size(waveformMain_,1);
    shiftsMain(1:nComps_,n:n+sumGood-1) = shiftsMain_;
    ampMain(n:n+sumGood-1) = ampMain_;
    medCCsMain(n:n+sumGood-1) = medCCsMain_;
    tMain(n:n+sumGood-1) = tMain_;
    bazMain(n:n+sumGood-1) = bazMain_;
    velMain(n:n+sumGood-1) = velMain_;
    waveformMain(1:newWinlen,n:n+sumGood-1) = waveformMain_;

    %%
    n = n + sumGood - 1;
end

%%
shiftsMain(:,n+1:end) = [];
ampMain(n+1:end) = [];
medCCsMain(n+1:end) = [];
waveformMain(:,n+1:end) = [];
tMain(n+1:end) = [];
bazMain(n+1:end) = [];
velMain(n+1:end) = [];

%%
if plotFlag
    close all;
    nUsed = sum(isfinite(shiftsMain))';
    goodI = velMain >= 205 & velMain <= 475 & ampMain >= ampThresh & ...
        medCCsMain >= thresh & nUsed >= 3;

    figure();
    plot(tMain(goodI),velMain(goodI),'o'); zoom on; grid on; hold on;
    plot(tMain(~goodI),velMain(~goodI),'.');title("Velocity");


    figure();
    plot(tMain(goodI),medCCsMain(goodI),'o'); zoom on; grid on; hold on;
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
    subplot(211);
    SS = scatter(tMain(goodI),ampMain(goodI),[],bazMain(goodI),'filled');
    c1 = colorbar;
    caxis([0 360]);
    zoom on; grid on; title("Amplitudes");
    ax = gca; ax.YScale = 'log';
    SS.MarkerFaceAlpha = 0.5;
    SS.MarkerEdgeColor = 'w';
    SS.MarkerEdgeAlpha = 0.5;
    ylabel("$Pa$");
    c1.Label.String = "$back-azimuth$";
    c1.Label.Interpreter = 'latex';

    subplot(212);
    nHours = 1;
    rate = t2r(tMain(goodI),hours(nHours));
    plot(tMain(goodI),rate,'.'); zoom on; grid on; title("Hourly Rate");
    c2 = colorbar;
    c2.Visible = 'off';

    fprintf("QC'd Events: <strong>%d</strong>, Discarded Events: <strong>%d</strong>, \n",sum(goodI),sum(~goodI));
end
