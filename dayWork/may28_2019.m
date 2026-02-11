%may28_2019

tStart = datetime(2019,05,20);
tEnd = datetime(2019,05,31);
referenceStartTime = datetime(2019,05,26);

newFs = 32;
secDur = 2^12;
maxLag = 32;
dailyFlag = false;
lfc = 2;
hfc = 4;

stnm1 = "BRTU"; chan1 = "HHN"; net1 = "EC"; locid1 = "";
stnm2 = stnm1; chan2 = "HHN"; net2 = net1; locid2 = locid1;
SNCL1 = [stnm1,chan1,net1,locid1];
SNCL2 = [stnm2,chan2,net2,locid2];
[dayStack,t,caus,acaus,symStack,lags] = greenFromNoise(tStart,tEnd,SNCL1,SNCL2,lfc,hfc,newFs,secDur,maxLag,dailyFlag);

%%
close all;
secStart = 2; %4/lfc;
cwiDur = 10; %/lfc;
secEnd = secStart + cwiDur;

tic;
refTrace = nanmean(caus(:,t <= referenceStartTime),2);
testTrace = caus;
[prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
figure('units','normalized','outerposition',[0 0 1 0.75]);
ax(1) = subplot(311);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([0 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
toc;
disp('finished first');

refTrace = nanmean(acaus(:,t <= referenceStartTime),2);
testTrace = acaus;
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    [prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
end
ax(2) = subplot(312);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([0 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
toc;
disp('finished second');

refTrace = nanmean(symStack(:,t <= referenceStartTime),2);
testTrace = symStack;
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    [prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
end
ax(3) = subplot(313);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([0 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
linkaxes(ax,'x');
toc;
disp('finished third');

ax(1).YGrid = 'on';
ax(2).YGrid = 'on';
ax(3).YGrid = 'on';

lags1 = (0:size(symStack,1)-1)'/newFs;
figure('units','normalized','outerposition',[0 0 1 0.75]);
imagesc(lags1,(1:length(t))',symStack'); colorbar; zoom on;
caxis(0.5*0.5*[-0.1 0.1]);
figure('units','normalized','outerposition',[0 0 1 0.75]);
imagesc(lags,(1:length(t))',dayStack'); colorbar; zoom on;
caxis(0.5*0.5*[-0.1 0.1]);
