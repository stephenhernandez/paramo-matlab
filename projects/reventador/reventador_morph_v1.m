clear; close all;
cd ~/Desktop/;
T = readtable('new_temperaturas.xlsx');

%%
tImage = T.FechaUTC_new_HoraEnArchivo;
T_south = T.VS__C_;
T_north = T.VN__C_;
T_nw = T.VNW__C_;
T_central = T.VC;
notas = string(T.nota);

%%
Temps = [T_south T_central T_north T_nw];

%%
colNames = ["Vento Sur";"Vento Central";"Vento Norte";"Vento Noroccidental"];
lineWidth = 2;
ncols = size(Temps,2); 
maxHeight = 20;
maxCumHeight = 1100;
figure();
for col = 1:ncols-1
    ventI = Temps(:,col);
    ventI = isfinite(ventI);
    tVent = tImage(ventI);
    ax(col) = subplot(ncols,1,col);
    [N,edges_] = histcounts(tVent,...
        (dateshift(min(tVent),'start','day'):dateshift(max(tVent),'end','day'))');
    edges_ = edges_(1:end-1);
    N = N';
    nI = N>0;
    plot(edges_(nI),N(nI),'o','LineWidth',lineWidth); hold on; grid on;
    ax(col).YLim = [0 maxHeight];
    yyaxis right;
    plot(tVent,(1:sum(ventI))','k.','LineWidth',lineWidth+2);
    %ax(col).YLim = [0 maxCumHeight];
    ax(col).YAxis(2).Color = 'k';
    ax(col).YAxis(2).Label.String = 'Num. Acum.';
    ylim_ = ax(col).YLim;
    ax(col).YLim = [min(ylim_) max(ylim_)+1];
    ax(col).YAxis(1).Label.String = 'Num. Diario';
    titStr = sprintf("%s",colNames(col)); 
    title(ax(col),titStr);
end
linkaxes(ax,'x');

%%
close all; 
figure(); 
load('tGood_tmp','tGood');
histogram(tGood,(dateshift(min(tGood),'start','day'):dateshift(max(tGood),'end','day'))'); zoom on; grid on;
ax = gca;
xlim(ax,[datetime(2018,04,01) datetime(2019,12,01)]);
ylim([0 100])
ax.YAxis(1).Label.String = 'Num. Diario';
cd ~/Desktop/
goodI = tGood>=datetime(2018,04,01) & tGood<=datetime(2019,12,01);
yyaxis right; plot(ax,tGood(goodI),0:sum(goodI)-1,'-','linewidth',4);
ax.YAxis(2).Label.String = 'Num. Acum.';