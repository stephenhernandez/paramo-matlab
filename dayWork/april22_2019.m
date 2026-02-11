%april22_2019
clear;
close all;
clc;

cd ~;
%cd '~/research/now/sn_eruption/noiseCorrelations'
%cd '~/cr/2012/INDI/'

dayStart = datetime(2016,11,25); %001;
dayEnd = datetime(2016,12,05); %010;
referenceStartTime = datetime(2016,12,01);

days = dayStart:dayEnd;
lDays = length(days);
lfc = 2;
hfc = 4;
npoles = 8;
tw = 0.001;
plotFlag = true;
dailyFlag = 0;

newFs = 32;
secDur = 2^8;
totN = secDur*newFs;
mxl = totN - 1;
stackN = 2*mxl+1;
nOverlap = 0;

stnm1 = "VTCE"; chan1 = "HHN"; net1 = "OV"; locid1 = "";
stnm2 = stnm1; chan2 = "HHN"; net2 = net1; locid2 = locid1;
SNCL1 = [stnm1,chan1,net1,locid1];
SNCL2 = [stnm2,chan2,net2,locid2];

zeroPhaseFlag = 1;
ampFact = 20;
normFlag = false;
smoother = 2./lfc;
lags = (-mxl:mxl)';
lags = lags/newFs;
maxLag = 128;
mI = abs(lags) <= maxLag;
sum_mI = sum(mI);
maxWindows = floor(86400/secDur);

%%
if dailyFlag
    dayStack = NaN(stackN,lDays);
    parfor i = 1:lDays
        dayStack_ = getCF(days(i),SNCL1,SNCL2,newFs,lfc,hfc,npoles,tw,zeroPhaseFlag,smoother,totN,nOverlap,normFlag,maxWindows,stackN,true);
        dayStack(:,i) = dayStack_;
    end
    t = days;
    dayStack = dayStack(mI,:);
    goodI = find(sum(~isnan(dayStack)));
    dayStack = dayStack(:,goodI);
    t = t(goodI);
    
    Nsmooth = 1;
    box = ones(Nsmooth,1)/Nsmooth;
    dayStack = flipud(convn(flipud(dayStack'),box));
    dayStack = dayStack(Nsmooth:end,:);
    dayStack = normalize_traces(dayStack');
else
    totSubWindows = sum_mI * maxWindows;
    dayData = NaN(sum_mI,lDays*maxWindows);
    dayData2 = NaN(totSubWindows,lDays);
    t = NaT(maxWindows,lDays);
    si = (1:maxWindows:totSubWindows)';
    ei = flipud((totSubWindows:-maxWindows:maxWindows)');
    parfor i = 1:lDays
        disp(datestr(days(i)));
        [dayData_,t_] = getCF(days(i),SNCL1,SNCL2,newFs,lfc,hfc,npoles,tw,zeroPhaseFlag,smoother,totN,nOverlap,normFlag,maxWindows,stackN,false);
        dayData_ = dayData_(mI,:);
        dayData_ = dayData_(:);
        dayData2(:,i) = dayData_;
        t(:,i) = t_;
    end
    
    %%
    for i = 1:lDays
        dayData_ = dayData2(:,i);
        dayData(:,si(i):ei(i)) = reshape(dayData_,sum_mI,maxWindows);
    end
    t = t(:);
    
    %%
    Nsmooth = 16;
    box = ones(Nsmooth,1)/Nsmooth;
    dayData = flipud(convn(dayData',box));
    dayData = flipud(dayData(Nsmooth:end,:));
    dayStack = normalize_traces(dayData');
    clear dayData*
end

%%
lags = lags(mI);
dayStack = normalize_traces(dayStack);

%%
[lDayStack,lDays] = size(dayStack);

if mod(lDayStack,2) == 0 %if even
    disp('hmmm... its possible something went wrong...');
    return
else
    newMXL = (lDayStack-1)/2;
    midpoint = newMXL+1;
    startInd1 = 1; endInd1 = midpoint;
    startInd2 = midpoint; endInd2 = lDayStack;
    
    totalStack = nanmean(dayStack,2);
    [maxPos,maxPosI] = max(totalStack);
    [maxNeg,maxNegI] = max(-totalStack);
    
    if maxNeg > maxPos
        trueMidPoint = maxNegI;
    else
        trueMidPoint = maxPosI;
    end
    
    if trueMidPoint ~= midpoint
        disp('the traces arent exactly symmetric...')
        diffMidPoints = sign(trueMidPoint - midpoint);
        absOMP = abs(trueMidPoint - midpoint);
        midpoint = trueMidPoint;
        
        if diffMidPoints+1
            startInd1 = 2*absOMP+1;             % acausal side
            endInd1 = midpoint;
            startInd2 = midpoint;               % causal side
            endInd2 = lDayStack;
        else
            startInd1 = 1;                      % acausal side
            endInd1 = trueMidPoint;
            startInd2 = trueMidPoint;           % causal side
            endInd2 = 2*trueMidPoint-1;
        end
    end
    
    disp(['midpoint = ',num2str(midpoint)]);
    acaus = flipud(dayStack(startInd1:endInd1,:));
    caus = dayStack(startInd2:endInd2,:);
    
    acaus = normalize_traces(acaus);
    caus = normalize_traces(caus);
end

symStack = normalize_traces(caus(1:size(acaus,1),:)+acaus);

%symStack = downsample(symStack',Nsmooth)';
%t = downsample(t,Nsmooth);

%%
close all;
secStart = 4/lfc;
cwiDur = 16/lfc;
secEnd = secStart + cwiDur;

tic;
refTrace = nanmean(caus(:,t <= referenceStartTime),2);
testTrace = caus;
[prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
figure('units','normalized','outerposition',[0 0 1 0.75]);
ax(1) = subplot(311);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([-1 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
toc;
disp('finished first');

refTrace = nanmean(acaus(:,t <= referenceStartTime),2);
testTrace = acaus;
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    [prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
end
ax(2) = subplot(312);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([-1 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
toc;
disp('finished second');

refTrace = nanmean(symStack(:,t <= referenceStartTime),2);
testTrace = symStack;
if ~strcmp(chan1,chan2) || ~strcmp(stnm1,stnm2)
    [prcntVector,peakCC] = cwiShifts(refTrace,testTrace,secStart,secEnd,newFs);
end
ax(3) = subplot(313);
Sh = scatter(t,prcntVector,[],peakCC,'filled'); colorbar; caxis([-1 1]); zoom on;
Sh.MarkerEdgeColor = 'k';
linkaxes(ax,'x');
toc;
disp('finished third');

figure('units','normalized','outerposition',[0 0 1 0.75]);
imagesc((0:size(symStack,1)-1)'/newFs,(1:lDays)',symStack'); colorbar; zoom on;
caxis(0.5*0.5*[-0.1 0.1]);
figure('units','normalized','outerposition',[0 0 1 0.75]);
imagesc(lags,(1:lDays)',dayStack'); colorbar; zoom on;
caxis(0.5*0.5*[-0.1 0.1]);

% [pxx,fxx] = pmtm(dayStack,[],[],32);
% [~,maxPxxI] = max(pxx);
% maxPxx = fxx(maxPxxI);
% figure(); plot(t,maxPxx,'o'); zoom on;

%%
newMax = 20*newFs;
newRef = testTrace(1:newMax,:); %preallocate
shifts_ = base + base*prcntVector/100;
base = 10000;
parfor i = 1:length(prcntVector)
    disp(i);
    shift_ = shifts_(i);
    newTest = resample(testTrace(:,i),shift_,base);
    newRef(:,i) = newTest(1:newMax);
end
disp('i am done');
toc;

%
[prcntVector,peakCC] = cwiShifts(normalize_traces(nanmean(newRef(:,t <= referenceStartTime),2)),testTrace,secStart,secEnd,newFs);
prcntVector = prcntVector(goodI); peakCC = peakCC(goodI);
figure(); plot(prcntVector,peakCC,'o'); zoom on;
figure(); plot(peakCC,'o'); zoom on;
figure(); plot(t,prcntVector,'o'); zoom on;
toc;
