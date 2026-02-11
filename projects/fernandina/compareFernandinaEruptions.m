%compareFernandinaEruptions
% the idea here is to plot rsam at both stations for all these eruptions
% (2018,2019,2020) and align them according to their eruption time and see
% if there are any patterns

close all;
eruptionStartTimes = [datetime(2017,09,04,18,22,00); ...
    datetime(2018,06,16,17,05,00); ...
    datetime(2020,01,12,23,59,60)];
lfc = 0.25;
hfc = 4;
diffFlag = 1;
dur = 15;
S = rmsGather(datetime(2017,09,03),3,dur,lfc,hfc,"FER2","HHZ",0,1,"EC","",diffFlag);
S(2) = rmsGather(datetime(2018,06,15),3,dur,lfc,hfc,"FER1","BHZ",0,1,"EC","",diffFlag);
S(3) = rmsGather(datetime(2020,01,11),3,dur,lfc,hfc,"FER1","BHZ",0,1,"EC","",diffFlag);
% S = rmsGather(datetime(2020,01,11),3,60,2,6,"FER1","BHZ",false,false,"EC","");

%%
close all;
Ssmooth = medfiltWaveforms(S,40,0,true);

figure('units','normalized','outerposition',[0 0 1 1]);
lS = length(S);
for i = 1:lS
    Scut(i) = cutWaveforms(Ssmooth(i),eruptionStartTimes(i)-hours(12),0,24*3600);
    t = getTimeVec(Scut(i));
    t = t-hours(6);
    d = Scut(i).d;
    tI = t >= datetime(2000,01,12);
    a(i) = subplot(lS,1,i);
    hh = plot(t(tI),d(tI),'linewidth',1);
    axis tight
    ax = gca;
    ax.YScale = 'log';
    xlabel('Tiempo Galapagos');
    ylabel('Amplitud [cuentas]');
    %legend('FER1.BHZ.EC')
    hh.LineWidth = 2;
    %grid on;
    %ax.YLim = [2 500];
    zoom on;
end

%%
yyyy = ["2017";"2018";"2020"];
ax = plotWaveforms(Scut,[],[],[],[],true);
for i = 1:length(ax)
    ax(i).YScale = 'log';
    ax(i).XLim = [0 24]-12;
    ylims = ax(i).YLim;
    hold(ax(i),'on');
    plot(ax(i),[0 0],ylims,'k--','linewidth',1);
    grid(ax(i),'on');
    ax(i).YLabel.String = yyyy(i);
    ax(i).FontSize = 18;
end
