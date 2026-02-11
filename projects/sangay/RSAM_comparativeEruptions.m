clear; close all; %clc;
exponent = 1;
maxHours1 = 6.3;
medianFlag = false;
secDur = 60;
stride = 60*1/4;
rmsFlag = true;
newFs = 25;
units = 'vel';
direction = true;
scalar = 1e9;
tw = 0.002;
derivedRate = false;
verboseFlag = false;
lfc = 1/100;
sagalfc = 3/8;
sagahfc = 3;
reglfc = 0.6;
reghfc = 1.2;
noiseWinHours = 0.4;
maxHours = maxHours1 - noiseWinHours;

tStarts = [datetime(2020,06,08,23,12,00);...
    datetime(2020,09,20,07,35,00);...
    datetime(2021,03,06,02,45,00);...
    datetime(2021,03,11,08,18,00);...
    datetime(2021,04,12,23,25,00);...
    datetime(2021,05,07,13,35,00);...
%     datetime(2021,05,21,13,45,00);...
    datetime(2021,05,30,21,50,00)];

%%
tic;
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 09 jun 2020\n');
S09 = loadWaveforms(datetime(2020,06,08),2,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
S09 = detrendWaveforms(differentiateWaveforms(S09));
S09 = scaleWaveforms(transferWaveforms(S09,lfc,-inf,4,newFs,units,direction),scalar);
S09(1) = filterWaveforms(S09(1),sagalfc,sagahfc);
S09(2:end) = filterWaveforms(S09(2:end),reglfc,reghfc);matlab
S09 = intWaveforms(taperWaveforms(S09,tw));
S09 = interpolateWaveforms(S09);
S09env = S09;
S09env = powWaveforms(S09env,2);
if medianFlag
    S09env = medfiltWaveforms(S09env,secDur*newFs+1);
else
    S09env = convWaveforms(S09env,secDur*newFs);
end
S09env = powWaveforms(S09env,1/2);
S09env = resampleWaveforms(S09env,1/stride);
S09env = cutWaveforms(S09env,tStarts(1),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 09 jun 2020\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 20 sep 2020\n');
S20 = loadWaveforms(datetime(2020,09,19),2,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
S20 = detrendWaveforms(differentiateWaveforms(S20));
S20 = scaleWaveforms(transferWaveforms(S20,lfc,-inf,4,newFs,units,direction),scalar);
S20(1) = filterWaveforms(S20(1),sagalfc,sagahfc);
S20(2:end) = filterWaveforms(S20(2:end),reglfc,reghfc);
S20 = intWaveforms(taperWaveforms(S20,tw));
S20env = S20;
S20env = powWaveforms(S20env,2);
if medianFlag
    S20env = medfiltWaveforms(S20env,secDur*newFs+1);
else
    S20env = convWaveforms(S20env,secDur*newFs);
end
S20env = powWaveforms(S20env,1/2);
S20env = resampleWaveforms(S20env,1/stride);
S20env = cutWaveforms(S20env,tStarts(2),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 20 sep 2020\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 23 jan 2021\n');
% S23 = loadWaveforms(datetime(2021,01,23),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
% S23 = detrendWaveforms(differentiateWaveforms(S23));
% S23 = scaleWaveforms(transferWaveforms(S23,lfc,-inf,4,newFs,units,direction),scalar);
% S23(1) = filterWaveforms(S23(1),sagalfc,sagahfc);
% S23(2:end) = filterWaveforms(S23(2:end),reglfc,reghfc);
% S23 = intWaveforms(taperWaveforms(S23,tw));
% S23env = S23;
% S23env = powWaveforms(S23env,2);
% S23env = convWaveforms(S23env,secDur*newFs);
% S23env = powWaveforms(S23env,1/2);
% S23env = resampleWaveforms(S23env,1/stride);
% S23env = cutWaveforms(S23env,dateshift(S23env(1).ref,'start','day')+hours(13)+minutes(40),0,hours(maxHours));
% toc;
% fprintf(1,'done with 23 jan 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 06 mar 2021\n');
S06 = loadWaveforms(datetime(2021,03,05),2,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
S06 = detrendWaveforms(differentiateWaveforms(S06));
S06 = scaleWaveforms(transferWaveforms(S06,lfc,-inf,4,newFs,units,direction),scalar);
S06(1) = filterWaveforms(S06(1),sagalfc,sagahfc);
S06(2:end) = filterWaveforms(S06(2:end),reglfc,reghfc);
S06 = intWaveforms(taperWaveforms(S06,tw));
S06env = S06;
S06env = powWaveforms(S06env,2);
if medianFlag
    S06env = medfiltWaveforms(S06env,secDur*newFs+1);
else
    S06env = convWaveforms(S06env,secDur*newFs);
end
S06env = powWaveforms(S06env,1/2);
S06env = resampleWaveforms(S06env,1/stride);
S06env = cutWaveforms(S06env,tStarts(3),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 06 mar 2021\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 11 mar 2021\n');
S11 = loadWaveforms(datetime(2021,03,11),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
S11 = detrendWaveforms(differentiateWaveforms(S11));
S11 = scaleWaveforms(transferWaveforms(S11,lfc,-inf,4,newFs,units,direction),scalar);
S11(1:end) = filterWaveforms(S11(1:end),reglfc,reghfc);
S11 = intWaveforms(taperWaveforms(S11,tw));
S11env = S11;
S11env = powWaveforms(S11env,2);
if medianFlag
    S11env = medfiltWaveforms(S11env,secDur*newFs+1);
else
    S11env = convWaveforms(S11env,secDur*newFs);
end
S11env = powWaveforms(S11env,1/2);
S11env = resampleWaveforms(S11env,1/stride);
S11env = cutWaveforms(S11env,tStarts(4),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 11 mar 2021\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 04 apr 2021\n');
% S04 = loadWaveforms(datetime(2021,04,04),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
% S04 = detrendWaveforms(differentiateWaveforms(S04));
% S04 = scaleWaveforms(transferWaveforms(S04,lfc,-inf,4,newFs,units,direction),scalar);
% S04(1:end) = filterWaveforms(S04(1:end),reglfc,reghfc);
% S04 = intWaveforms(taperWaveforms(S04,tw));
% S04env = S04;
% S04env = powWaveforms(S04env,2);
% S04env = convWaveforms(S04env,secDur*newFs);
% S04env = powWaveforms(S04env,1/2);
% S04env = resampleWaveforms(S04env,1/stride);
% S04env = cutWaveforms(S04env,dateshift(S04env(1).ref,'start','day')+hours(14)+minutes(00),0,hours(maxHours));
% toc;
% fprintf(1,'done with 04 apr 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 12 apr 2021\n');
S12 = loadWaveforms(datetime(2021,04,12),2,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
S12 = detrendWaveforms(differentiateWaveforms(S12));
S12 = scaleWaveforms(transferWaveforms(S12,lfc,-inf,4,newFs,units,direction),scalar);
S12(1:end) = filterWaveforms(S12(1:end),reglfc,reghfc);
S12 = intWaveforms(taperWaveforms(S12,tw));
S12env = S12;
S12env = powWaveforms(S12env,2);
if medianFlag
    S12env = medfiltWaveforms(S12env,secDur*newFs+1);
else
    S12env = convWaveforms(S12env,secDur*newFs);
end
S12env = powWaveforms(S12env,1/2);
S12env = resampleWaveforms(S12env,1/stride);
S12env = cutWaveforms(S12env,tStarts(5),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 12 apr 2021\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 28 apr 2021\n');
% S28 = loadWaveforms(datetime(2021,04,28),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"],"EC","",derivedRate,verboseFlag);
% S28 = detrendWaveforms(differentiateWaveforms(S28));
% S28 = scaleWaveforms(transferWaveforms(S28,lfc,-inf,4,newFs,units,direction),scalar);
% S28(1) = filterWaveforms(S28(1),sagalfc,sagahfc);
% S28(2:end) = filterWaveforms(S28(2:end),reglfc,reghfc);
% S28 = intWaveforms(taperWaveforms(S28,tw));
% S28env = S28;
% S28env = powWaveforms(S28env,2);
% S28env = convWaveforms(S28env,secDur*newFs);
% S28env = powWaveforms(S28env,1/2);
% S28env = resampleWaveforms(S28env,1/stride);
% S28env = cutWaveforms(S28env,dateshift(S28env(1).ref,'start','day')+hours(14)+minutes(00),0,hours(maxHours));
% toc;
% fprintf(1,'done with 28 apr 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 07 may 2021\n');
S07 = loadWaveforms(datetime(2021,05,07),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"]);
S07 = detrendWaveforms(differentiateWaveforms(S07));
S07 = scaleWaveforms(transferWaveforms(S07,lfc,-inf,4,newFs,units,direction),scalar);
S07(1) = filterWaveforms(S07(1),sagalfc,sagahfc);
S07(2:end) = filterWaveforms(S07(2:end),reglfc,reghfc);
S07 = intWaveforms(taperWaveforms(S07,tw));
S07env = S07;
S07env = powWaveforms(S07env,2);
if medianFlag
    S07env = medfiltWaveforms(S07env,secDur*newFs+1);
else
    S07env = convWaveforms(S07env,secDur*newFs);
end
S07env = powWaveforms(S07env,1/2);
S07env = resampleWaveforms(S07env,1/stride);
S07env = cutWaveforms(S07env,tStarts(6),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 07 may 2021\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 16 may 2021\n');
% S16 = loadWaveforms(datetime(2021,05,16),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"]);
% S16 = detrendWaveforms(differentiateWaveforms(S16));
% S16 = scaleWaveforms(transferWaveforms(S16,lfc,-inf,4,newFs,units,direction),scalar);
% S16(1:end) = filterWaveforms(S16(1:end),reglfc,reghfc);
% S16 = intWaveforms(taperWaveforms(S16,tw));
% S16env = S16;
% S16env = powWaveforms(S16env,2);
% S16env = convWaveforms(S16env,secDur*newFs);
% S16env = powWaveforms(S16env,1/2);
% S16env = resampleWaveforms(S16env,1/stride);
% S16env = cutWaveforms(S16env,dateshift(S16env(1).ref,'start','day')+hours(12)+minutes(25),0,hours(maxHours));
% toc;
% fprintf(1,'done with 16 may 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');
%
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 18 may 2021\n');
% S18 = loadWaveforms(datetime(2021,05,18),1,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"]);
% S18 = detrendWaveforms(differentiateWaveforms(S18));
% S18 = scaleWaveforms(transferWaveforms(S18,lfc,-inf,4,newFs,units,direction),scalar);
% S18(1) = filterWaveforms(S18(1),sagalfc,sagahfc);
% S18(2:end) = filterWaveforms(S18(2:end),reglfc,reghfc);
% S18 = intWaveforms(taperWaveforms(S18,tw));
% S18env = S18;
% S18env = powWaveforms(S18env,2);
% S18env = convWaveforms(S18env,secDur*newFs);
% S18env = powWaveforms(S18env,1/2);
% S18env = resampleWaveforms(S18env,1/stride);
% S18env = cutWaveforms(S18env,dateshift(S18env(1).ref,'start','day')+hours(00)+minutes(00),0,hours(maxHours));
% toc;
% fprintf(1,'done with 18 may 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');
% 
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'starting 21 may 2021\n');
% S21 = loadWaveforms(datetime(2021,05,21),2,["PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"]);
% S21 = detrendWaveforms(differentiateWaveforms(S21));
% S21 = scaleWaveforms(transferWaveforms(S21,lfc,-inf,4,newFs,units,direction),scalar);
% S21(1:end) = filterWaveforms(S21(1:end),reglfc,reghfc);
% S21 = intWaveforms(taperWaveforms(S21,tw));
% S21env = S21;
% S21env = powWaveforms(S21env,2);
% if medianFlag
%     S21env = medfiltWaveforms(S21env,secDur*newFs+1);
% else
%     S21env = convWaveforms(S21env,secDur*newFs);
% end
% S21env = powWaveforms(S21env,1/2);
% S21env = resampleWaveforms(S21env,1/stride);
% S21env = cutWaveforms(S21env,tStarts(7),-hours(noiseWinHours),hours(14));
% toc;
% fprintf(1,'done with 21 may 2021\n');
% fprintf(1,'-----------------------------------------------------------\n');
% fprintf(1,'\n');
% fprintf(1,'\n');

fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'starting 30 may 2021\n');
S30 = loadWaveforms(datetime(2021,05,30),2,["SAGA";"PUYO";"TAIS";"PKYU";"BULB";"TAMH";"PORT";"BPAT";"BMAS";"BRUN"],["HHZ";"BHZ"]);
S30 = detrendWaveforms(differentiateWaveforms(S30));
S30 = scaleWaveforms(transferWaveforms(S30,lfc,-inf,4,newFs,units,direction),scalar);
S30(1) = filterWaveforms(S30(1),sagalfc,sagahfc);
S30(2:end) = filterWaveforms(S30(2:end),reglfc,reghfc);
S30 = intWaveforms(taperWaveforms(S30,tw));
S30env = S30;
S30env = powWaveforms(S30env,2);
if medianFlag
    S30env = medfiltWaveforms(S30env,secDur*newFs+1);
else
    S30env = convWaveforms(S30env,secDur*newFs);
end
S30env = powWaveforms(S30env,1/2);
S30env = resampleWaveforms(S30env,1/stride);
S30env = cutWaveforms(S30env,tStarts(7),-hours(noiseWinHours),hours(maxHours));
toc;
fprintf(1,'done with 30 may 2021\n');
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'\n');
fprintf(1,'\n');

%%
close all;
clear sumvelsquared durations stdamp meanamp maxamp medamp hax axPUYO madamp;
TAIS = [S09env(3); S20env(3); S06env(3); S11env(2); S12env(1); S07env(2); S30env(2)];

[masterFig,axPUYO] = plotWaveforms(TAIS,[],[],[],[],true);

% figure('units','normalized','outerposition',[0.72 0 0.28 1]);
% lTais = length(TAIS);
% for i = 1:lTais
%     axPUYO(i) = subplot(lTais,1,i);
%     t = getTimeVec(TAIS(i));
%     d = TAIS(i).d;
%
%     tStart = tStarts(i);
%     tEnd = tEnds(i);
%
%     tI = t >= tStart & t <= tEnd;
%     d_ = d(tI);
%     t_ = t(tI);
%     plot(axPUYO(i),hours(t-min(t)),d,'k-','linewidth',1);
%     hold on; zoom on; %grid on;
%     axPUYO(i).ColorOrderIndex = 1;
%     aa = area(axPUYO(i),hours(t-min(t)),d,'FaceAlpha',0.5,'EdgeColor','none');
%     if i < lTais
%         axPUYO(i).XTickLabel = '';
%     end
%     grid on; %grid minor;
% end
linkaxes(axPUYO,'xy');

masterFig.OuterPosition = [0.72 0 0.28 1];
tits = flipud(["30 may 2021"; "07 may 2021";"13 abr 2021";"11 mar 2021";"06 mar 2021";"20 sep 2020";"09 jun 2020"]);
sup_title_ = "TAIS.HHZ, 0.6 - 1.2 Hz.";
if ~medianFlag
    sup_title_ = strcat(sup_title_,", root mean sq., ",num2str(secDur)," sec. win., $nm \cdot s^{-1}$");
else
    sup_title_ = strcat(sup_title_,", root med. sq.",num2str(secDur)," sec. win., $nm \cdot s^{-1}$");
end
title(axPUYO(1),sup_title_);
xlim([0 maxHours1]);

%sensorFloor = [5; 4; 5; 5; 6; 5; 7]; %nm/s sensor floors (RMedS)
sensorFloor = [7; 6; 7; 8; 9; 8; 11]; %nm/s sensor floors (RMeanS)

exponent = 1; %should always be 1
dMaster = [];
timeMaster = [];

% tStarts = [datetime(2020,06,08,23,12,00);...
%     datetime(2020,09,20,07,35,00);...
%     datetime(2021,03,06,02,45,00);...
%     datetime(2021,03,11,08,18,00);...
%     datetime(2021,04,12,23,25,00);...
%     datetime(2021,05,07,13,35,00);...
%     datetime(2021,05,30,21,50,00)];

tEnds = [datetime(2020,06,09,04,15,00);...
    datetime(2020,09,20,12,08,00);...
    datetime(2021,03,06,08,35,00);...
    datetime(2021,03,11,12,10,00);...
    datetime(2021,04,13,02,16,00);...
    datetime(2021,05,07,18,35,00);...
%     datetime(2021,05,22,01,40,00);...
    datetime(2021,05,31,03,30,00)];

Titles = ["Erupcion 08/09 Junio:";...
    "Erupcion 20 Septiembre:";...
    "Erupcion 05/06 Marzo:";...
    "Erupcion 11 Marzo:";...
    "Erupcion 12/13 Abril:";...
    "Erupcion 07 Mayo:";...
%     "Erupcion 21/22 Mayo:";...
    "Erupcion 30/31 Mayo:"];

tt = getTimeVec(TAIS(end));
fprintf('Ultima actualizacion: %s\n',datestr(tt(end)));

zoomFlag = true;
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(tStarts)
    sensorFloor_ = sensorFloor(i);
    t = getTimeVec(TAIS(i));
    d = TAIS(i).d;

    tStart = tStarts(i);
    tEnd = tEnds(i);

    tI = t >= tStart & t <= tEnd;
    d_ = d(tI);
    t_ = t(tI);

    mainLength = 10*3600/stride;
    ltt = length(t_);
    if ltt >= mainLength
        timeMaster = [timeMaster t(1:mainLength)];
        dMaster = [dMaster d(1:mainLength)];
    else
        timeMaster = [timeMaster NaT(mainLength,1)];
        dMaster = [dMaster NaN(mainLength,1)];
        timeMaster(1:ltt,end) = t_;
        dMaster(1:ltt,end) = d_;
    end

    fprintf('%s\n',Titles(i));
    durations(i,1) = stride*sum(d_>=sensorFloor_ & t_ >= tStart & t_ <= tEnd)/3600;
    sumvelsquared(i,1) = (nansum(abs(d_(d_>=sensorFloor_)).^2))*(stride/secDur);
    meanamp(i,1) = nanmean(abs(d_(d_>=sensorFloor_) - 0).^exponent);
    maxamp(i,1) = max(abs(d_(d_>=sensorFloor_) - 0).^exponent);
    stdamp(i,1) = nanstd(abs(d_(d_>=sensorFloor_) - 0).^exponent);
    medamp(i,1) = nanmedian(abs(d_(d_>=sensorFloor_) - 0).^exponent);
    madamp(i,1) = mad(abs(d_(d_>=sensorFloor_) - 0).^exponent,1);
    fprintf('Max Amp.: <strong>%5.2f</strong>; Mean Amp: <strong>%5.2f</strong>; Std. Amp: <strong>%5.2f</strong>; Med. Amp: <strong>%5.2f</strong>\n',...
        maxamp(i,1),meanamp(i,1),stdamp(i,1),medamp(i,1));


    fprintf('Duration: %f horas\n',durations(i));
    fprintf('Sum Vel: %u\n',round(nansum(abs(d_(d_>=sensorFloor_) - 0).^1)));
    fprintf('Sum Vel^2: %u\n',sumvelsquared(i));
    fprintf('Estimated sensor floor: %f\n',sensorFloor_);

    fprintf('\n');
    axPUYO(i).ColorOrderIndex = 2;
    hold(axPUYO(i),'on');
    %aa = area(axPUYO(i),hours(t-min(t)),d,'FaceAlpha',0.5,'EdgeColor','none');
    %aa.BaseLine.Visible = 'off';
    %axPUYO(i).Legend.Visible = 'off';
    aa = area(axPUYO(i),[0 maxHours1],sensorFloor_.*[1 1],'FaceAlpha',0.35,'EdgeColor','none','FaceColor','k');
    %axPUYO(i).ColorOrderIndex = axPUYO(i).ColorOrderIndex + 1;

    grid(axPUYO(i),'on');
    if zoomFlag
        cOrder = axPUYO(i).ColorOrder;
        scatter(axPUYO(i),hours(tStart-min(t)),interp1(hours(t-min(t)),d,hours(tStart-min(t))),100,cOrder(2,:),'o','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
        scatter(axPUYO(i),hours(tEnd-min(t)),2+interp1(hours(t-min(t)),d,hours(tEnd-min(t))),100,cOrder(3,:),'v','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
        text(axPUYO(i),0.08,45,Titles(i),'FontSize',14);
        %axPUYO(i).YLabel.String = '$nm \cdot s^{-1}$';
        %axPUYO(i).YLim = [0 50];
        axPUYO(i).YScale = 'log';
        axPUYO(i).YLim = [4 100];
    else
        text(axPUYO(i),0.1,500,tits(i),'FontSize',14);
        cOrder = axPUYO(i).ColorOrder;
        scatter(axPUYO(i),hours(tStart-min(t)),interp1(hours(t-min(t)),d,hours(tStart-min(t))),60,cOrder(2,:),'o','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
        scatter(axPUYO(i),hours(tEnd-min(t)),2+interp1(hours(t-min(t)),d,hours(tEnd-min(t))),60,cOrder(3,:),'v','filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k');
    end

    %
    figure(2);
    hax(i,1) = subplot(2,4,i);
    H = histogram(abs(d_(d_>=sensorFloor_)),logspace(log10(1),log10(700),101));
    binEdges = H.BinEdges;

    H.EdgeColor = 'none';
    zoom on;
    hold on;
    grid on; grid minor;
    lwidth = 1;
    plot(hax(i),maxamp(i,1)*ones(2,1),[0 0+interp1(binEdges(1:end-1),H.BinCounts,maxamp(i,1))],'-','linewidth',lwidth,'Color',cOrder(2,:));
    plot(hax(i),medamp(i,1)*ones(2,1),[0 0+interp1(binEdges(1:end-1),H.BinCounts,medamp(i,1))],'-','linewidth',lwidth,'Color',cOrder(3,:));
    plot(hax(i),meanamp(i,1)*ones(2,1),[0 3+interp1(binEdges(1:end-1),H.BinCounts,meanamp(i,1))],'-','linewidth',lwidth,'Color',cOrder(4,:));

    scatter(hax(i),maxamp(i,1),2,110,cOrder(2,:),'v','filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.8);
    scatter(hax(i),medamp(i,1),2,110,cOrder(3,:),'^','filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.8);
    scatter(hax(i),meanamp(i,1),2,110,cOrder(4,:),'d','filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.8);
    %ylabel(hax(i),'N');
    %hax(i).YTickLabel = '';
    hax(i).YTick = [];
    xlabel(hax(i),'amp. [$nm \cdot s^{-1}$]','FontSize',14);
end
linkaxes(hax,'xy');
%axis(hax,'tight');
pause(2);
maxylim = 0;

clear htitle;
for i = 1:length(hax)
    htitle(i) = title(hax(i),{Titles(i);...
        ['$max$=\bf',num2str(maxamp(i,1)),', $med$=\bf',num2str(medamp(i,1))];...
        ['$\mu$=\bf',num2str(meanamp(i,1)),', $\sigma$=\bf',num2str(stdamp(i,1))]},...
        'FontSize',15);
    %hax(i).YScale = 'lin';
    hax(i).XScale = 'log';
    plot(hax(i),maxamp(i,1)*ones(2,1),[0 max(hax(i).YLim)],'-','linewidth',lwidth,'Color',cOrder(2,:));
    plot(hax(i),medamp(i,1)*ones(2,1),[0 max(hax(i).YLim)],'-','linewidth',lwidth,'Color',cOrder(3,:));
    plot(hax(i),meanamp(i,1)*ones(2,1),[0 max(hax(i).YLim)],'-','linewidth',lwidth,'Color',cOrder(4,:));
end

legax = subplot(2,4,8);
legax.Box = 'off';
hold(legax,'on');
scatter(legax,0.25,0.75,110,cOrder(2,:),'v','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
text(0.3,0.75,'max','FontSize',18);
scatter(legax,0.25,0.50,110,cOrder(3,:),'^','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
text(0.3,0.50,'median','FontSize',18);
scatter(legax,0.25,0.25,110,cOrder(4,:),'d','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
text(0.3,0.25,'mean','FontSize',18);
xlim([0 1]); ylim([0 1]);
legax.XAxis.Visible = 'off';
legax.YAxis.Visible = 'off';
sgtitle(sup_title_);

fontsize = 18;
clear pax;
figure('units','normalized','outerposition',[0 0 1 1]);
pax(1,1) = subplot(2,3,1);
scatter(pax(1),tStarts,durations,110,cOrder(1,:),'s','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('Hours','FontSize',fontsize); grid on; zoom on; hold on;
ylim([2 12]); title('Duration','FontSize',fontsize);

pax(2,1) = subplot(2,3,2);
scatter(pax(2),tStarts,sumvelsquared,110,cOrder(5,:),'^','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('$\propto$ Energy','FontSize',fontsize); grid on; zoom on; hold on;
pax(2).YScale = 'log';
ylim([2 40]*1e5);
title('Sum Of Velocity Squared','FontSize',fontsize);

pax(3,1) = subplot(2,3,3);
scatter(pax(3),tStarts,sumvelsquared./durations,110,cOrder(6,:),'o','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('$\propto$ Energy/Hour','FontSize',fontsize); grid on; zoom on; hold on;
pax(3).YScale = 'log';
ylim([0.4 10]*1e5);
title('Energy Density','FontSize',fontsize);

pax(4,1) = subplot(2,3,4);
scatter(pax(4),tStarts,meanamp,110,cOrder(4,:),'d','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('amp [$nm \cdot s^{-1}$]','FontSize',fontsize);
grid on; zoom on; hold on;
scatter(pax(4),tStarts,medamp,110,cOrder(3,:),'^','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
pax(4).YScale = 'log';
ylim([10 100]); title('Mean/Median RMS Amplitude','FontSize',fontsize);
legend(pax(4),'Mean','Median','Location','Northwest');

pax(5,1) = subplot(2,3,5);
scatter(pax(5),tStarts,maxamp,110,cOrder(2,:),'v','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('amp [$nm \cdot s^{-1}$]','FontSize',fontsize);
grid on; zoom on; hold on;
pax(5).YScale = 'log';
ylim([80 1200]); title('Max RMS Amplitude','FontSize',fontsize);

pax(6,1) = subplot(2,3,6);
scatter(pax(6),tStarts,stdamp,110,cOrder(7,:),'p','filled','MarkerFaceAlpha',0.85,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.85);
ylabel('amp [$nm \cdot s^{-1}$]','FontSize',fontsize);
grid on; zoom on; hold on;
pax(6).YScale = 'log';
ylim([10 100]); title('Standard Deviation of Amplitudes','FontSize',fontsize);

linkaxes(pax,'x');
sgtitle("Assorted Source Parameters of Sangay Eruptions",'FontSize',fontsize+2);

figure();
plot(tStarts,meanamp./stdamp,'.','MarkerSize',20); zoom on; grid on; hold on; %ylim([0.5 2]);
hold on;
plot(tStarts,medamp./madamp,'s','MarkerSize',20); zoom on; grid on; hold on; %ylim([0.5 2]);
legend('$\mu/\sigma$','med/MAD','Location','Best');
ylim([0.25 3]);

%%
% C = who('*env');
% %close all;
% for i = 1:length(C)
%     disp(i)
%     C_ = string(C(i));
%     disp(C_);
%     eval(strcat("[~,ax] = plotWaveforms(",C_,"); for j = 1:length(ax); ax(j).YScale = 'log'; end;"));
% end
