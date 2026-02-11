cd ~/research/now/sangay/
clear; close all;
load Puyo_tmp.mat

%%
pointsPerDay = 1440;
t = getTimeVec(S);
d = S.d; 
d(d == 0) = NaN;
plottingWeeks = 400;

tday = dateshift(t,'start','day');
tm = minutes(t - tday);
nWeeks = floor(length(tm)/pointsPerDay/7);
tm2 = reshape(tm(1:pointsPerDay*7*nWeeks),[pointsPerDay*7 nWeeks]);
d2 = reshape(d(1:pointsPerDay*7*nWeeks),[pointsPerDay*7 nWeeks]);
tthisweek = t(pointsPerDay*7*(nWeeks-plottingWeeks)+1:end);

%%
dmed = nanmedian(d2,2);
dmean = nanmean(d2,2);
bgBounds = repmat(dmed,1,2);
for i = 1:size(d2,1)
    Y = prctile(d2(i,:),[10 90]);
    bgBounds(i,:) = Y;
end

%%
dmed = [repmat(dmed,nWeeks,1); dmed(1:length(d)-pointsPerDay*7*nWeeks)];

figure('units','normalized','outerposition',[0 0 1 1]);
dcut = d(pointsPerDay*7*(nWeeks-plottingWeeks)+1:end);
pp = plot(tthisweek,dcut,'-','linewidth',4); 
pp.Color(4) = 0.5;
zoom on; grid on; 

hold on;
ax = gca;
ax.YScale = 'log';

dmedcut = dmed(1:length(tthisweek));
plot(tthisweek,dmedcut,'-','linewidth',3); zoom on; grid on;

bgBounds2 = [repmat(bgBounds(:,1),nWeeks,1); bgBounds(1:length(d)-pointsPerDay*7*nWeeks,1)];
bgBounds2(:,2) = [repmat(bgBounds(:,2),nWeeks,1); bgBounds(1:length(d)-pointsPerDay*7*nWeeks,2)];

bgCut = bgBounds2(1:length(tthisweek),:);
ppB(1) = plot(tthisweek,bgCut(:,1),'-','linewidth',2,'Color',[0.5 0.5 0.5]); zoom on; grid on;
ppB(2) = plot(tthisweek,bgCut(:,2),'-','linewidth',2,'Color',[0.5 0.5 0.5]); zoom on; grid on;

%%
ppB(1).Color(4) = 0.8;
ppB(2).Color(4) = 0.8;
xlim([datetime(2020,09,20,08,00,00) datetime(2020,09,20,12,00,00)]);

%%
outFile = '~/amplitude_puyo_20SEP2020.txt';
tI = tthisweek >= datetime(2020,09,20) & tthisweek < datetime(2020,09,21);

tprint = datenum(tthisweek(tI))-693960; %convert to 1900 excel format

formatSpec = '%f,%f,%f,%f';
str = compose(formatSpec,tprint,dcut(tI),bgCut(tI,1),bgCut(tI,2));
str = string(str);
fileID = fopen(outFile,'w');
fprintf(fileID,'%s\n',str);
fclose(fileID);
%!open ~/ampltude_puyo_20SEP2020.txt;

figure(1);
cd ~/research/now/sangay/
fname = 'puyo_amplitudes_20SEP2020';
print('-djpeg',fname);
print('-depsc',fname);

%%
% % clear; close all; clc;
% % cd ~/products/rsam/;
% %
% % %
% % load EC.BULB..BHZ_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat
% % W = S;
% % %load EC.PUYO..HHZ_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat
% % load ~/research/now/sangay/Puyo_tmp.mat
% % %S.ref = S.ref - minutes(8.5);
% % W(2,1) = S;
% % load EC.TAMH..HHZ_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat
% % W(3,1) = S;
% % load EC.TAIS..HHZ_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat
% % W(4,1) = S;
% % W = padWaveforms(W);
% % close all; [~,ax] = plotWaveforms(W);
% % for i = 1:length(ax)
% %     ax(i).YScale = 'log';
% % end
% %
% % %%
% % snclList = ["BULB" "BHZ" "EC" "";...
% %     "PUYO" "HHZ" "EC" "";...
% %     "TAMH" "HHZ" "EC" "";...
% %     "TAIS" "HHZ" "EC" ""];
% %
% % %%
% % kstnms = snclList(:,1);
% % refEllipse = referenceEllipsoid('wgs84');
% % [stlat,stlon,stelev] = metaDataFromStationList(kstnms);
% % [d_,az_] = distance(-2.00535,-78.341294,stlat,stlon,refEllipse);
% % d_ = d_*1e-3;
% %
% % sensitivities = [4.872110e+08,...
% %     3.141950e+08,...
% %     3.141950e+08,...
% %     3.141950e+08,...
% %     3.141950e+08];
% %
% %
% % S_09JUN2020 = cutWaveforms(W,datetime(2020,06,08),0,minutes(1440*3),true);
% % S_20SEP2020 = cutWaveforms(W,datetime(2020,09,19),0,minutes(1440*3),true);
% % S_23JAN2021 = cutWaveforms(W,datetime(2021,01,23),0,minutes(1440*3),true);
% %
% % %% nanometers/sec
% % for i = 1:length(kstnms)
% %     S_09JUN2020(i) = scaleWaveforms(S_09JUN2020(i),1e9/sensitivities(i));
% %     S_20SEP2020(i) = scaleWaveforms(S_20SEP2020(i),1e9/sensitivities(i));
% %     S_23JAN2021(i) = scaleWaveforms(S_23JAN2021(i),1e9/sensitivities(i));
% % end
% %
% % %%
% % [~,ax1] = plotWaveforms(S_09JUN2020);
% % [~,ax2] = plotWaveforms(S_20SEP2020);
% % [~,ax3] = plotWaveforms(S_23JAN2021);
% %
% % for i = 1:length(kstnms)
% %     ax1(i).YScale = 'log';
% %     ax2(i).YScale = 'log';
% %     ax3(i).YScale = 'log';
% % end
% %
% %
% % %%
% % % ampFact = -1e7;
% % %
% % % figure();
% % % hold on;
% % % for i = 1:length(S_20SEP)
% % %     t_ = getTimeVec(S_20SEP(i));
% % %     pp(i) = plot(t_,d_(i)+ampFact*d_20SEP(:,i)/sensitivities(i),'linewidth',3);
% % %     pp(i).Color(4) = 0.8;
% % % end
% % % ax = gca;
% % % ylim([50 105]);
% % % ax.YDir = 'reverse';
% % % zoom on;
% % % legend(kstnms,'location','northwest');

