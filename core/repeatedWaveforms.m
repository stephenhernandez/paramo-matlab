function [Sf,ccnorm,pks,locs,amps,templateStack] = ...
    repeatedWaveforms(S,varargin)
nVarargin = length(varargin);
functionDefaults = {...
    datetime(2022,04,20,07,43,26),...
    24,...
    0.2,...
    0.8,...
    16,...
    6.5,...
    0.44};

optsToUse = functionDefaults;
optsToUse(1:nVarargin) = varargin;
[templateStart,templateEnd,lfc,hfc,newFs,thresh,minThresh] = deal(optsToUse{:});

%%
Sf = filterWaveforms(S,lfc,hfc);
Sf = resampleWaveforms(Sf,newFs);
Sf = syncWaveforms(Sf);

if isstruct(templateStart)
    Scut = templateStart;
    clear templateStart;
    templates = double(pull(Scut));
    templates = templates(1:1+newFs*templateEnd,:);
    [winlen,lT] = size(templates);
else
    Scut = cutWaveforms(Sf,templateStart,0,templateEnd);
    templates = double(pull(Scut));
    [winlen,lT] = size(templates);
end

lS = length(S);
if lS ~= lT
    fprintf(2,'Number of templates not the same as number of inputs\n');
    return;
end

templates = templates(1:winlen,:);
templates = normalizeWaveforms(detrend(templates));

t = getTimeVec(Sf);
npts = Sf(1).npts;
ccnorm = NaN(npts,lT);
for i = 1:lT
    d = Sf(i).d;
    p = templates(:,i);
    p = p./norm(p);
    p = flipud(p);
    cc = fftfilt(p,d);
    winlen = length(p);
    norms = fftfilt(ones(winlen,1),d.^2);
    norms = sqrt(abs(norms));
    ccnorm(:,i) = cc./norms;
end

%thresh = 10;
%minThresh = 0.45;
minAmp = 5e2;

ccnorm = mean(ccnorm,2,'omitnan');
mad_ = mad(ccnorm,1);
[pks,locs] = findpeaks(ccnorm,'MINPEAKDISTANCE',0.25*winlen,'MINPEAKHEIGHT',thresh*mad_);

badLI = locs <= winlen | locs > npts - winlen + 1;

locs(badLI) = [];
pks(badLI) = [];
locs = locs - winlen + 1;
tabs = t(locs);

Sdup1 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(1),tabs+seconds(0),0,minutes(1)))));
Sdup2 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(2),tabs+seconds(0),0,minutes(1)))));
Sdup3 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(3),tabs+seconds(0),0,minutes(1)))));

amps = median([peak2peak(double(pull(Sdup1)))' peak2peak(double(pull(Sdup2)))' peak2peak(double(pull(Sdup3)))']./2,2);
goodI = amps >= minAmp & pks >= minThresh; 
amps = amps(goodI); 
tabs = tabs(goodI); 
pks = pks(goodI); 
locs = locs(goodI);

Sdup1 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(1),tabs+seconds(0),0,minutes(1)))));
Sdup2 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(2),tabs+seconds(0),0,minutes(1)))));
Sdup3 = detrendWaveforms((detrendWaveforms(cutWaveforms(Sf(3),tabs+seconds(0),0,minutes(1)))));


[P,tpratio,~,~] = pratio(tabs,51);
close all; 
figure(); 
plot(t(1:end-winlen+1),ccnorm(winlen:end)); zoom on; grid on;
hold on;
plot(tabs,pks,'.'); zoom on; grid on;

figure('units','normalized','outerposition',[0 0 1 1]); 
plot(tpratio,P,'.'); zoom on; title('p-ratio'); grid on;

ddup1 = double(pull(Sdup1));
ddup2 = double(pull(Sdup2));
ddup3 = double(pull(Sdup3));
normers = sqrt(max([doAutoCorrFreqDom(detrend(ddup1)); ...
    doAutoCorrFreqDom(detrend(ddup2));...
    doAutoCorrFreqDom(detrend(ddup3))]))';

figure('units','normalized','outerposition',[0 0 1 1]); 
d2 = ddup1; 
imagesc((0:size(d2,1)-1)'/newFs,(1:size(d2,2))',sign(d2)'); zoom on; colorbar; colormap gray;

figure('units','normalized','outerposition',[0 0 1 1]); 
d2 = ddup2; 
imagesc((0:size(d2,1)-1)'/newFs,(1:size(d2,2))',sign(d2)'); zoom on; colorbar; colormap gray;


figure('units','normalized','outerposition',[0 0 1 1]); 
d2 = ddup3; 
imagesc((0:size(d2,1)-1)'/newFs,(1:size(d2,2))',sign(d2)'); zoom on; colorbar; colormap gray;
clear d2

figure(); plot(tabs,pks,'p'); zoom on; axis tight; grid on;
[P,~,meanP,stdP] = pratio(log10(amps(1:end-1)),51);

figure(); 
plot(tpratio,meanP,'.'); zoom on; grid on;

figure(); 
plot(tpratio,stdP,'.'); zoom on; grid on;

figure(); 
plot(tpratio,P,'.'); zoom on; grid on;

%%
if ~sum(~isfinite(amps))
    nHours = 2; 
    [rate,meanAmps] = t2r(tabs,hours(nHours),log10(amps));

    figure('units','normalized','outerposition',[0 0 1 1]); 
    aXX(1) = subplot(211); semilogy(tabs,amps,'.'); zoom on; grid on; hold on; 
    semilogy(tabs,10.^meanAmps,'linewidth',4); title('Average Amplitude'); zoom on; grid on; 
    
    aXX(2) = subplot(212); 
    semilogy(tabs,rate/nHours,'.'); grid on; title('Average Hourly Rate'); linkaxes(aXX,'x');

    figure('units','normalized','outerposition',[0 0 1 1]); 
    semilogy(tabs,rate.*(10.^meanAmps)/nHours,'.'); zoom on; grid on; title('$\propto$ Energy');
    
    figure('units','normalized','outerposition',[0 0 1 1]); 
    ax(1) = subplot(311); plot(templates(:,1)); 
    ax(2) = subplot(312); plot(templates(:,2)); 
    ax(3) = subplot(313); plot(templates(:,3)); 
    zoom on; grid on; axis tight; linkaxes(ax,'x');
end

%%
ddup1 = double(pull(Sdup1)); 
ddup2 = double(pull(Sdup2)); 
ddup3 = double(pull(Sdup3));

ddup1 = ddup1 ./ normers'; 
ddup2 = ddup2 ./ normers'; 
ddup3 = ddup3 ./ normers';

%whos
templatesNew = [mean(ddup1,2) mean(ddup2,2) mean(ddup3,2)];
figure('units','normalized','outerposition',[0 0 3/5 1]);
ax(1) = subplot(4,3,[1 2 3]); hold on;
plot(templatesNew(:,1),'linewidth',2);
ax(2) = subplot(4,3,[4,5,6]); hold on; plot(templatesNew(:,2),'linewidth',2);
ax(3) = subplot(4,3,[7,8,9]); hold on; plot(templatesNew(:,3),'linewidth',2);
zoom on; grid on; axis tight; linkaxes(ax(1:3),'x'); figure(13);
ax(4) = subplot(4,3,10); hold on; plot(templatesNew(:,3),templatesNew(:,2),'linewidth',2); axis equal; grid on;
ax(5) = subplot(4,3,11); hold on; plot(templatesNew(:,2),templatesNew(:,1),'linewidth',2); axis equal; grid on;
ax(6) = subplot(4,3,12); hold on; plot(templatesNew(:,3),templatesNew(:,1),'linewidth',2); axis equal; grid on;

%%
templateStack = NaN(winlen,lT);
for i = 1:lT
    Sdup1 = detrendWaveforms(cutWaveforms(Sf(i),tabs+seconds(0),0,minutes(1)));
    ddup1 = pws(double(pull(Sdup1)));
    templateStack(:,i) = normalizeWaveforms(ddup1(1:winlen,:));
end

