%function SAG1_infrasound_array_processing_v2(dayStart,dayEnd,dayInc)
clear; close all; clc;
dayStart = datetime(2024,07,09);
dayEnd = datetime(2024,07,09);
dayInc = 1;

saveFlag = false;
plotFlag = true;

dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);

vI = [];
nUsed = [];
tMain = [];
ampMain = [];
bazMain = [];
velMain = [];
meanCCsMain = [];
medCCsMain = [];

minNumStations = 3;
kholelist = ["01";"02";"03";"04";"05"];
31for jj = 1:lDays
    day_ = dayVec(jj);

    fprintf("loading day: %s\n",day_);
    S = populateWaveforms(nComps);
    for i = 1:nComps
        khole_ = kholelist(i);
        S_ = loadWaveforms(day_,dayInc,"SAG1","BDF","EC",khole_); %this order cannot change
        if day_ >= datetime(2024,07,01) && day_ < datetime(2025,03,13)
            if i == 4
                S_ = scaleWaveforms(S_,-150);
            end
        end
        S(i) = S_;
    end

    refs = pull(S,'ref');
    badComponents = isnat(refs);
    goodComponents = ~badComponents;
    lS = sum(goodComponents);
    if lS < minNumStations
        fprintf(2,'no data for: %s\n',day_);
        continue;
    end

    [vI_,nUsed_,tMain_,ampMain_,bazMain_,velMain_,meanCCsMain_,...
        medCCsMain_,~,badShifts_] = ...
        process_sag1_infrasound_array_v1(S);
    vI = [vI; vI_];
    vI = logical(vI);
    nUsed = [nUsed; nUsed_];
    tMain = [tMain; tMain_];
    ampMain = [ampMain; ampMain_];
    bazMain = [bazMain; bazMain_];
    velMain = [velMain; velMain_];
    meanCCsMain = [meanCCsMain; meanCCsMain_];
    medCCsMain = [medCCsMain; medCCsMain_];
end

%%
if saveFlag
    cd ~/masa/old/research/now/sangay/
    % clearvars -except *Main lfc hfc secDur ampThresh thresh Fs Fs2 ...
    %     minNumStations nUsed vI Gcolumns
    % clear waveformMain

    ampMain_ = ampMain;
    bazMain_ = bazMain;
    meanCCsMain_ = meanCCsMain;
    medCCsMain_ = medCCsMain;
    nUsed_ = nUsed;
    tMain_ = tMain;
    velMain_ = velMain;
    vI_ = vI;

    load('SAG1_infrasoundArrayCatalog_400mHz_6400mHz.mat');
    ampMain = [ampMain; ampMain_];
    bazMain = [bazMain; bazMain_];
    meanCCsMain = [meanCCsMain; meanCCsMain_];
    medCCsMain = [medCCsMain; medCCsMain_];
    nUsed = [nUsed; nUsed_];
    tMain = [tMain; tMain_];
    velMain = [velMain; velMain_];
    vI = [vI; vI_];

    [~,sortI] = sort(tMain);
    ampMain = ampMain(sortI);
    bazMain = bazMain(sortI);
    meanCCsMain = meanCCsMain(sortI);
    medCCsMain = medCCsMain(sortI);
    nUsed = nUsed(sortI);
    tMain = tMain(sortI);
    velMain = velMain(sortI);
    vI = logical(vI(sortI));

    [tMain,ampMain,meanCCsMain,bazMain,medCCsMain,velMain,vI,nUsed] = ...
        filterCatalog(tMain,ampMain,3,meanCCsMain,bazMain,medCCsMain,velMain,vI,nUsed);

    % clearvars -except *Main lfc hfc secDur ampThresh thresh Fs Fs2 ...
    %     minNumStations nUsed vI Gcolumns plotFlag
    clear waveformMain
    save('SAG1_infrasoundArrayCatalog_400mHz_6400mHz.mat','tMain','ampMain',...
        'meanCCsMain','bazMain','medCCsMain','velMain','vI','nUsed',...
        'minNumStations','-v7.3');
end

%%
close all;
tic;
causalFlag = true;
ampThresh = 0.01;
thresh = 0.7;
vI = velMain >= 290 & velMain <= 390 & ampMain >= ampThresh & ...
    medCCsMain >= thresh & nUsed >= minNumStations & bazMain >= 275 & ...
    bazMain <= 345;

if plotFlag
    %vI = vI & ampMain >= ampThresh;
    figure();
    plot(tMain(vI),velMain(vI),'.'); zoom on; grid on;
    title("Velocity");

    figure();
    plot(tMain(vI),medCCsMain(vI),'.'); zoom on; grid on; hold on;
    %plot(tMain(vI),meanCCsMain(vI),'.'); zoom on; grid on;
    legend("Median","Mean");

    figure();
    bazMain(bazMain < 0) = bazMain(bazMain < 0) + 360;
    SS = scatter(tMain(vI),bazMain(vI),20*exp(log10(ampMain(vI))),velMain(vI),'filled');
    SS.MarkerFaceAlpha = 0.5;
    SS.MarkerEdgeColor = "k";
    SS.MarkerEdgeAlpha = 0.5;
    zoom on; grid on; %title("BAZ");
    ax = gca;
    ax.YLim = [0 360];
    colorbar;
    %ax.YTick = (0:30:360)';

    figure('units','normalized','outerposition',[0 0 0.6 1]);
    nSubplots = 2;
    tiledlayout(nSubplots,1, 'Padding', 'compact', 'TileSpacing', 'none');
    ax(1) = nexttile(); %subplot(211);
    semilogy(tMain(vI),ampMain(vI),'.'); zoom on; grid on; ylabel('Amplitud [$Pa$]');
    ax(2) = nexttile(); %subplot(212);
    nHours = 1;
    rate = t2r(tMain(vI),hours(nHours),[],causalFlag);
    plot(tMain(vI),rate/nHours,'.'); zoom on; grid on;
    ylabel("Rate [Hourly]");
    ax(1).XTickLabel = [];
    ax(2).YAxisLocation = "right";
    linkaxes(ax,'x');

    figure();
    histogram(ampMain(vI),logspace(floor(min(log10(ampMain(vI)))),ceil(max(log10(ampMain(vI)))),201));
    ax = gca; ax.XScale = 'log'; grid on; zoom on;
    ylabel("frequency");
    xlabel("amplitude [$Pa$]");

    [rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = ...
        t2r(tMain(vI),hours(nHours),ampMain(vI),causalFlag);
    figure();
    axx(1) = subplot(311);
    plot(tMain(vI),rate/nHours,'.'); grid on; hold on;
    axx(2) = subplot(312);
    semilogy(tMain(vI),meanMagsFixedTimeWin,'.'); grid on; hold on;
    axx(3) = subplot(313);
    %semilogy(tMain(vI),medianMagsFixedTimeWin.*rate/nHours,'.');
    semilogy(tMain(vI),sumEnergyFixedTimeWin,'.');
    hold on;
    semilogy(tMain(vI),(meanMagsFixedTimeWin.^2).*rate/nHours,'.');
    semilogy(tMain(vI),(medianMagsFixedTimeWin.^2).*rate/nHours,'.');
    zoom on; grid on; hold on;
    linkaxes(axx,'x');

    figure();
    histogram(bazMain(vI),linspace(0,360,361));
    ax = gca;  grid on; zoom on; ax.YScale = 'log';
    ylabel("frequency");
    xlabel("back azimuth [$^\circ$]");

    figure();
    histogram(velMain(vI),logspace(floor(min(log10(velMain(vI)))),ceil(max(log10(velMain(vI)))),201));
    ax = gca; ax.XScale = 'log'; grid on; zoom on;
    ylabel("frequency");
    xlabel("velocity [$m \cdot s^{-1}$]");

    if exist("Gcolumns",'var')
        if Gcolumns > 2
            figure();
            axx(1) = subplot(211);
            plot(axx(1),tMain(vI),bazMain(vI),'.'); zoom on; grid on; hold on;
            plot(axx(1),tMain(~vI),bazMain(~vI),'.');
            axx(2) = subplot(212);
            plot(axx(2),tMain(vI),belMain(vI),'.'); zoom on; grid on; hold on;
            plot(axx(2),tMain(~vI),belMain(~vI),'.');
            linkaxes(axx,'x');
        end
    end

    % tMain2 = tMain(vI);
    % vI2 = tMain2 >= datetime(2023,08,01);
    % figure(); SS = scatter(rate(vI2)/nDays,meanMagsFixedTimeWin(vI2),[],datenum(tMain2(vI2)),'filled'); zoom on; grid on; colorbar;
    % ax = gca; SS.MarkerEdgeColor = 'k';
    % ax.XScale = 'log'; ax.YScale = 'log'; SS.MarkerFaceAlpha = 0.5; SS.MarkerEdgeAlpha = 0.5;
    %
    % tMain2 = tMain(vI);
    % vI2 = tMain2 < datetime(2022,08,01);
    % figure(); SS = scatter(rate(vI2)/nDays,meanMagsFixedTimeWin(vI2),[],datenum(tMain2(vI2)),'filled'); zoom on; grid on; colorbar;
    % ax = gca; SS.MarkerEdgeColor = 'k';
    % ax.XScale = 'log'; ax.YScale = 'log'; SS.MarkerFaceAlpha = 0.5; SS.MarkerEdgeAlpha = 0.5;
    %
    % tMain2 = tMain(vI);
    % vI2 = tMain2 > datetime(2022,08,01) & tMain2 <= datetime(2023,05,01);
    % figure(); SS = scatter(rate(vI2)/nDays,meanMagsFixedTimeWin(vI2),[],datenum(tMain2(vI2)),'filled'); zoom on; grid on; colorbar;
    % ax = gca; SS.MarkerEdgeColor = 'k';
    % ax.XScale = 'log'; ax.YScale = 'log'; SS.MarkerFaceAlpha = 0.5; SS.MarkerEdgeAlpha = 0.5;
end
toc;

%%
% ampMain = ampMain(vI);
% bazMain = bazMain(vI);
% meanCCsMain = meanCCsMain(vI);
% medCCsMain = medCCsMain(vI);
% nUsed = nUsed(vI);
% tMain = tMain(vI);
% velMain = velMain(vI);
% clearvars -except ampMain bazMain meanCCsMain medCCsMain nUsed tMain velMain
% [rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = ...
% t2r(tMain,days(1),ampMain,true);
%%save('SAG1_InfrasoundCatalog_23OCT2023');