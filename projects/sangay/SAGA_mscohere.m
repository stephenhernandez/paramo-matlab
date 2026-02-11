clear; close all; clc;
%dayStart = datetime(2018,06,26);
%dayEnd = datetime(2018,06,26);
dayStart = datetime(2025,10,12);
dayEnd = datetime(2025,10,12);

dInc = 2;
dayVec = (dayStart:dayEnd)';
lDays = length(dayVec);
value = 0;

%
newFs = 40;
secDur = 5*60;
winlen = secDur*newFs;
nOverlap = 3/5;

%kstnm = ["BVC2";"BTAM"]; %"BREF";
kstnm = "VCH1";
kstnm = "REVS";
nf = 451;
lfc = 1/50;
hfc = 20;

for i = lDays:-1:1
    day_ = dayVec(i);
    S = loadWaveforms(day_,dInc,kstnm,["HHZ";"BDF"],"EC",["";"01"]);
    %S = loadWaveforms(day_,dInc,kstnm,["BHE"],"EC",["";"01"]);
    lS = length(S);
    if any(isnat(pull(S,'ref'))) || lS < 2
        fprintf('error with %s, skipping...\n',datestr(day_));
        continue;
    end

    fprintf('processing: %s\n',datestr(day_,29));
    S = syncWaveforms(S,false,true,true);
    S = resampleWaveforms(S,newFs);
    S = nanGapWaveforms(S,0);
    S = padWaveforms(S);
    %S = resampleWaveforms(S,newFs);
    S = interpolateWaveforms(S,value);

    %
    d1 = S(1).d;
    d2 = S(2).d;

    ld = min([length(d1) length(d2)]);
    if ld < winlen
        fprintf('error with %s, length not long enough, skipping...\n',datestr(day_));
        continue;
    end

    [d1,~,endIndex] = cutWindows(d1,winlen,nOverlap,true);
    disp(size(d1));
    d1 = detrend((d1));

    d2 = cutWindows(d2,winlen,nOverlap,true); 
    disp(size(d2)); 
    d2 = detrend((d2));

    t = getTimeVec(S(1));
    tdumb = t(endIndex);

    minT = min(tdumb);
    maxT = max(tdumb);

    maxTicks = 13;
    totDur = hours(maxT - minT);
    if 0.65*floor(totDur) >= 24
        units = 'day';
        unitDur = 1;
    elseif 0.65*floor(totDur) >= 1
        units = 'hour';
        unitDur = 1/24;
    elseif 0.65*floor(totDur) > 1/60
        units = 'minute';
        unitDur = 1/60;
    else
        units = 'second';
        unitDur = 1/3600;
    end

    tickStart = datenum(dateshift((minT),'start',units,'nearest'));
    tickEnd = datenum(dateshift((maxT),'start',units,'nearest'));
    nUnits = round((tickEnd - tickStart)/unitDur);
    while nUnits > maxTicks
        nUnits = ceil(nUnits/2);
    end
    nTicks = nUnits + 1;
    uniformTicks = linspace(tickStart,tickEnd,nTicks);
    %sss = seconds(tEnd-min(tStart));

    %%
    disp('start mscohere');
    tic; [cxy,f] = mscohere(d1,d2,[],[],logspace(log10(lfc),log10(hfc),nf),newFs);
    toc;
    disp('stop mscohere');

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig.Visible = 'on';
    h = pcolor(datenum(tdumb),f,(cxy)); 
    axis xy; 
    c2 = colorbar;
    h.EdgeColor = 'none'; 
    ax = gca;
    ax.YScale = 'log'; 
    set(ax,'ColorScale','log')
    clim([1/10 1]); zoom on;
    datetick('x','keeplimits');

    ax(1).XTick = uniformTicks;
    ax(1).XTickLabel = datestr(ax(1).XTick);
    ax(1).XTickLabelRotation = 16;
    title(strcat(datestr(minT)," - ",datestr(maxT)));

    fname_ =strcat(kstnm,'_mscohere_',datestr(day_,'yyyy.mm.dd'),'.jpg');
    fname = fullfile('~','research','now','sangay','mscohere',fname_);
    fprintf('figure name: %s\n',fname);
    %print('-djpeg',fname);
end
