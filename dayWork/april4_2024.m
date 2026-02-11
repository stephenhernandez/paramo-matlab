clear; close all; 
[tsheet,type,Tsp,amp,magnitude,codaDuration,period,refStation] = loadVolcanicEvents("cotopaxi",7);
clearvars -except tsheet amp type Tsp magnitude codaDuration period refStation

[tsheetF,ampF,typeF,TspF,magnitudeF,codaDurationF,periodF,refStationF] = ...
    filterCatalog(tsheet,amp,10,type,Tsp,magnitude,codaDuration,period,refStation);

goodI = tsheetF >= datetime(2019,01,01);
axL1 = linkedPlot(tsheetF(goodI),ampF(goodI),t2r(tsheetF(goodI),days(1),[],false));
axL1(1).YScale = 'log';
axL1(1).YLabel.String = 'amplitude';
axL1(2).YLabel.String = 'daily rate'; axis tight;

goodI = tsheetF >= datetime(2019,01,01) & ampF >= 3e2;
axL2 = linkedPlot(tsheetF(goodI),ampF(goodI),t2r(tsheetF(goodI),days(1),[],false));
axL2(1).YScale = 'log';
axL2(1).YLabel.String = 'amplitude';
axL2(2).YLabel.String = 'daily rate'; axis tight;

[N,edges] = histcounts(tsheetF(goodI),dateshift(min(tsheetF(goodI)),'start','day'):dateshift(max(tsheetF(goodI)),'end','day'));
N = N';
edges = edges(1:end-1)';

T = readtable('~/igdata/GasCotopaxi_03APR2024.xlsx');

[lia,locb] = ismember(T.Var1,edges);
tGas = T.Var1;
correctedFlux = T.Var9;
tGas = tGas(lia);
correctedFlux = correctedFlux(lia);
 
postI = edges >= datetime(2023,08,01);
eruptionI = edges >= datetime(2022,08,01) & edges < datetime(2023,08,01);
close all; figure(); ax(1) = subplot(211); SS = scatter(N(postI),correctedFlux(postI),[],datenum(edges(postI)),'filled'); zoom on; grid on;
SS.MarkerEdgeColor = 'k'; SS.MarkerFaceAlpha = 0.7; SS.MarkerEdgeAlpha = 0.7; cbar1 = colorbar; clim([datenum(2022,08,01) datenum(2024,04,05)]); cbar1.TickLabels = datestr(cbar1.Ticks);
hold on; ax(2) = subplot(212);
SS = scatter(N(eruptionI),correctedFlux(eruptionI),[],datenum(edges(eruptionI)),'v','filled'); zoom on; grid on;
SS.MarkerEdgeColor = 'k'; SS.MarkerFaceAlpha = 0.7; SS.MarkerEdgeAlpha = 0.7;
cbar2 = colorbar; clim([datenum(2022,08,01) datenum(2024,04,05)]); cbar2.TickLabels = datestr(cbar2.Ticks);
%
ax(1).XScale = 'log'; ax(1).YScale = 'log'; ax(2).XScale = 'log'; ax(2).YScale = 'log';
linkaxes(ax,'xy');
ax(1).XLabel.String = 'Numero de Sismos/Dia (amp > 300 cuentas)';
ax(2).XLabel.String = 'Numero de Sismos/Dia (amp > 300 cuentas)';
ax(1).YLabel.String = 'Flujo SO2';
ax(2).YLabel.String = 'Flujo SO2';
ax(1).Title.String = {'Pos-Erupcion';'t > 2023-08-01'};
ax(2).Title.String = {'Syn-Erupcion';'2022-08-01 < t < 2023-08-01'};
%
nMeasurements = T.Var3;
nMeasurements = nMeasurements(lia);
axL3 = linkedPlot(edges,N,correctedFlux,nMeasurements);
axL3(1).Title.String = 'Numero de Sismos/Dia (amp >= 300 cuentas)';
axL3(2).Title.String = 'Flujo SO2';
axL3(3).Title.String = 'Numero de Medidas Validas';