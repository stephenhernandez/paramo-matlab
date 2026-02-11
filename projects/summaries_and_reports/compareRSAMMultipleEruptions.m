%compareFernandinaEruptions
% the idea here is to plot rsam at both stations for all these eruptions
% (2018,2019,2020) and align them according to their eruption time and see
% if there are any patterns

cd ~/research/now/fernandina/
% datetime(2013,07,14,11,46,00); % Tungu (RETU)
% datetime(2013,10,18,14,26,00); % Tungu (RETU)
% datetime(2014,02,01,22,39,00); % Tungu (RETU)
% datetime(2015,08,14,09,02,00); % Coto (BREF)
% datetime(2016,02,26,17,11,00); % Tungu (RETU)
% datetime(2017,09,04,18,22,00); % Fer (FER2)
% datetime(2018,06,16,17,05,00); % Fer (FER1)
% datetime(2018,06,26,19,40,00); % SN (VCH1)
% datetime(2019,05,07,20,40,00); % Sangay (SAGA)
% datetime(2020,01,13,00,00,00); % Fer (FER1)

close all;
eruptionStartTimes = [datetime(2013,07,14,11,46,00); ...    % Tungu (RETU)
    datetime(2013,10,18,14,26,00);...                       % Tungu (RETU)
    datetime(2014,02,01,22,39,00);...                       % Tungu (RETU)
    datetime(2015,08,14,09,02,00);...                       % Coto (BREF)
    datetime(2016,02,26,17,11,00);...                       % Tungu (RETU)
    datetime(2017,09,04,18,22,00);...                       % Fer (FER2)
    datetime(2018,06,16,17,05,00);...                       % Fer (FER1)
    datetime(2018,06,26,19,40,00);...                       % SN (VCH1)
    datetime(2019,05,07,20,40,00);...                       % Sangay (SAGA)
    datetime(2020,01,13,00,00,00)];                         % Fer (FER1)
disp(eruptionStartTimes)

% eruptionStartTimes = [datetime(2017,09,04,18,22,00); ...
%     datetime(2018,06,16,17,05,00); ...
%     datetime(2020,01,12,23,59,60)];
% lfc = 0.25;
% hfc = 4;
% diffFlag = 1;
% dur = 15;
% S = rmsGather(datetime(2017,09,03),3,dur,lfc,hfc,"FER2","HHZ",0,1,"EC","",diffFlag);
% S(2) = rmsGather(datetime(2018,06,15),3,dur,lfc,hfc,"FER1","BHZ",0,1,"EC","",diffFlag);
% S(3) = rmsGather(datetime(2020,01,11),3,dur,lfc,hfc,"FER1","BHZ",0,1,"EC","",diffFlag);
% % S = rmsGather(datetime(2020,01,11),3,60,2,6,"FER1","BHZ",false,false,"EC","");

%%
close all;
nHours = 96;
method1 = true;
%Ssmooth = medfiltWaveforms(S,40,0,true);

%figure('units','normalized','outerposition',[0 0 1 1]);
lS = length(eruptionStartTimes);
for i = 1:lS
    if ismember(i,[1 2 3 5]) % RETU
        if method1
            load RETU_RSAM_MEDFILT.mat
            S = RETU_rsam_medfilt;
            clear RETU_rsam_medfilt;
        else
            S = extractWaveforms(tabs-seconds(5),seconds(65),"PINO","SHZ");
        end
    elseif i == 4 % BREF
        if method1
            load BREF_RSAM_MEDFILT.mat
            S = BREF_rsam_medfilt;
            clear BREF_rsam_medfilt;
        end
    elseif i == 6 % FER2
        if method1
            load FER2_RSAM_MEDFILT.mat
            S = FER2_rsam_medfilt;
            clear FER2_rsam_medfilt;
        end
    elseif i == 7 || i == 10 % FER1
        if method1
            load FER1_RSAM_MEDFILT.mat
            S = FER1_rsam_medfilt;
            clear FER1_rsam_medfilt;
        end
    elseif i == 8 % VCH1
        if method1
            load VCH1_RSAM_MEDFILT.mat
            S = VCH1_rsam_medfilt;
            clear VCH1_rsam_medfilt;
        end
    elseif i == 9 % SAGA
        if method1
            load SAGA_RSAM_MEDFILT.mat
            S = SAGA_rsam_medfilt;
            clear SAGA_rsam_medfilt;
        end
    end
    Ssmooth = medfiltWaveforms(S,40,0,true);
    Scut(i) = cutWaveforms(Ssmooth,eruptionStartTimes(i)-hours(nHours),0,(nHours*2)*3600);
    %     t = getTimeVec(Scut(i));
    %     t = t-hours(6);
    %     d = Scut(i).d;
    %     tI = t >= datetime(2000,01,12);
    %     a(i) = subplot(lS,1,i);
    %     hh = plot(t(tI),d(tI),'linewidth',1);
    %     axis tight
    %     ax = gca;
    %     ax.YScale = 'log';
    %     xlabel('Tiempo Galapagos');
    %     ylabel('Amplitud [cuentas]');
    %     %legend('FER1.BHZ.EC')
    %     hh.LineWidth = 2;
    %     %grid on;
    %     %ax.YLim = [2 500];
    %     zoom on;
end

%%
%yyyy = ["2017";"2018";"2020"];
ax = plotWaveforms(Scut,[],[],[],[],true);
for i = 1:length(ax)
    ax(i).YScale = 'log';
    ax(i).XLim = [0 (2*nHours)]-nHours;
    ylims = ax(i).YLim;
    hold(ax(i),'on');
    plot(ax(i),[0 0],ylims,'k--','linewidth',1);
    grid(ax(i),'on');
    %ax(i).YLabel.String = yyyy(i);
    ax(i).FontSize = 14;
end
