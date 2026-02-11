%recent eruptions at galapagos volcanoes
clear; close all;
clc;

% approximate eruption start times in Galapagos local time
axisMode = 'log';

% 2015, 2017, 2018, 2018, 2020, 2022
kstnms = ["FER1"; "FER2"; "FER1"; "VCH1"; "FER1"; "FER1"];
kcmpnms = ["BHZ"; "HHZ"; "BHZ"; "HHZ"; "BHZ"; "BHZ"];
lfc = [5; 2; 2; 1/2; 2; 5];
hfc = [7; 8; 8; 16; 8; 7];
dInc = 4;

%%
approximateEruptionStartTimes = [datetime(2015,05,25,01,00,00),...
    datetime(2017,09,04,12,25,00),...
    datetime(2018,06,16,11,10,00),...
    datetime(2018,06,26,13,35,00),...
    datetime(2020,01,12,18,05,00),...
    datetime(2022,01,06,23,15,00)];

lstarts = length(approximateEruptionStartTimes);

for i = 1:lstarts
    kstnm_ = kstnms(i);
    kcmpnm_ = kcmpnms(i);
    startTime_ = approximateEruptionStartTimes(i);
    lfc_ = lfc(i);
    hfc_ = hfc(i);

    %%
    startTime2 = dateshift(startTime_,'start','day')-1;
    S = loadWaveforms(startTime2,dInc,kstnm_,kcmpnm_,"EC");
    Sf = scaleWaveforms(transferWaveforms(S,lfc_,hfc_,4,100,'vel',true),1e6);

    %%
    [ampVec,winlen] = amplitudeVector(Sf.d,100,60,true,false,1);
    ampVec = downsample(ampVec,60*100);
    [S2,npts] = dealHeader(Sf,ampVec,1/60);
    t = getTimeVec(S2);
    t = t - hours(6);

    %%
    %figure();
    tI = t >= startTime_ - hours(24) & t < startTime_ + hours(24);
    maxValue = max(ampVec(tI));
    minValue = min(ampVec(tI));
    maxValue = ceil(log10(maxValue)) + 1;
    minValue = floor(log10(minValue)) - 1;

    figure('units','normalized','outerposition',[0 0 3/4 1]);
    lineWidth = 2;
    ll = plot([startTime_ startTime_],[10.^minValue 10.^maxValue],'k--','linewidth',lineWidth);
    ll.Color(4) = 0.75;
    hold on;

    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(t,ampVec,'.');
    zoom on; hold on; grid on;
    ll = plot(t,zpkFilter(ampVec,-inf,1/30,1,1,true),'linewidth',4);
    ll.Color(4) = 0.8;
    ax = gca;
    ax.YScale = axisMode;

    ylim([10^minValue 10^maxValue]);
    xlim([startTime_ - hours(24) startTime_ + hours(24)]);
end
