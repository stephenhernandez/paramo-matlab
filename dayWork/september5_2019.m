clear; clc; close all;
cd ~/research/now/tungurahua/dv

%retu_2011_2018_3Hz6Hz_4s36s
%fnames(1).name = 'retu_2008_2018_1p5Hz6Hz_5s20s.mat';

%cd retu_2011_2018_1Hz2Hz_6s38s
%cd retu_2011_2018_2Hz4Hz_4s20s
%cd retu_2009_2018_1Hz4Hz_4s36s
%cd retu_2009_2018_8s2Hz_10s80s
%cd retu_2009_2018_1s4s_10s80s
cd retu_2009_2018_1Hz4Hz_10s80s
fnames = dir('ictpRun_*.mat');

tF = [];
causF = [];
acausF = [];
symStackF = [];
dayStackF = [];
Nsmooth = 11;

%%
for i = 1:length(fnames)
    load(fnames(i).name)
    dayData = dayStack;
    box = ones(Nsmooth,1)/Nsmooth;
    dayData = flipud(convn(dayData',box));
    dayData = flipud(dayData(Nsmooth:end,:));
    dayStack = normalizeWaveforms(dayData');
    lDayStack = size(dayStack,1);
    dayStackF = [dayStackF dayStack];
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
        
        disp(['midpoint = ',num2str(midpoint)]);
        acaus = flipud(dayStack(startInd1:endInd1,:));
        caus = dayStack(startInd2:endInd2,:);
        
        acaus = normalizeWaveforms(acaus);
        caus = normalizeWaveforms(caus);
    end
    symStack = normalizeWaveforms(caus(1:size(acaus,1),:)+acaus);
    tF = [tF; t];
    causF = [causF caus];
    acausF = [acausF acaus];
    symStackF = [symStackF symStack];
end

%%
clc;
STNM1 = ["RETU","SHZ","EC",""];
STNM2 = STNM1;
referenceStartTime = datetime(2020,04,17);
cwiStart = 10;
cwiEnd = 100;
newFs = 32;
[dVcaus_,dVacaus_,dVsymmetric_,ccCaus_,ccAcaus_,ccSymmetric_] = ...
    getDV(tF,causF,acausF,symStackF,referenceStartTime,cwiStart,cwiEnd,newFs,STNM1(1),STNM2(1),STNM1(2),STNM2(2));

%%
close all;
Nmed = 41;
figure('units','normalized','outerposition',[0 0 1 0.75]);
scatter(tF,dVcaus_,[],ccCaus_,'o'); zoom on;
hold on;
y = medfiltSH(dVcaus_,Nmed,1);
plot(tF,y,'linewidth',2);
title('Causal');
colorbar;
caxis([0.4 0.6]);

figure('units','normalized','outerposition',[0 0 1 0.75]);
scatter(tF,dVacaus_,[],ccAcaus_,'o'); zoom on;
hold on;
y = medfiltSH(dVacaus_,Nmed,1);
plot(tF,y,'linewidth',2);
title('Acausal');
colorbar;
caxis([0.4 0.6]);

figure('units','normalized','outerposition',[0 0 1 0.75]);
scatter(tF,dVsymmetric_,[],ccSymmetric_,'o'); zoom on;
hold on;
y = medfiltSH(dVsymmetric_,Nmed,1);
plot(tF,y,'linewidth',2);
title('Symmetric');
colorbar;
caxis([0.4 0.6]);

figure('units','normalized','outerposition',[0 0 1 0.75]);
imagesc(lags,(1:length(tF))',dayStackF'); colorbar; zoom on;
caxis(0.5*0.5*[-0.001 0.001]);

%%
stack = nanmedian(causF,2);
stack = stack / norm(stack);
N = size(causF,2);
NFFT = 4096*2;
cc = max(ifft(fft(causF(:,1:N),NFFT).*fft(stack,NFFT),'symmetric'));

figure();
plot(tF(1:N),cc,'o');
zoom on;

figure();
plot(cc,dVcaus_,'o');
zoom on;
