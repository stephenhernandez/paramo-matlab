function [detection_times,overtoneFreq,harmStrength,fundPow,...
    t,dcut,overtonePow,interharmPow,overtoneI,powerOrig,pxx,fxx,fundPowOrig] = ...
    tremometer_control(S,nOvertones,MINMAX,MINMED,MINMIN,MINFREQ,MINPOW,MINPSD,plotFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Tremometer 1.1                                                                                    %
%   Heavily modified by Stephen Hernandez, November 2025                                              %
%                                                                                                     %
%   Matlab code to automatically detect and characterize harmonic tremor                              %
%   in continuous seismic data using a pitch-detection algorithm                                      %
%                                                                                                     %
%   Diana C. Roman                                                                                    %
%   Last modified July 31, 2017                                                                       %
%                                                                                                     %
%   Please cite: Roman, D.C. (2017), Automated detection and characterization of harmonic tremor      %
%   in continuous seismic data. Geophys. Res. Lett.,44,doi: 10.1002/2017GL073715.                     %
%                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if nargin < 3
    MINMAX = 10;
end

if nargin < 4
    MINMED = 4;
end

if nargin < 5
    MINMIN = 2;
end

if nargin < 6
    MINFREQ = 0.35;
end

if nargin < 7
    MINPOW = 1e-18;
end

if nargin < 8
    MINPSD = 1e-6;
end

if nargin < 9
    plotFlag = false;
end

% MINMAX = 10;
% MINMED = 3;
% MINMIN = 2;
% MINFREQ = 0.35;
% MINPOW = 1e-17;
% MINPSD = 1e-5;
% plotFlag = false;

%%
S = S(1);
delta = S.delta;
secDur = 120;
nOverlap = 0.5;
d = S.d;
npts = S.npts;
baseVec = delta*(0:npts-1)';
tStart = S.ref;
e = S.e;
tEnd = tStart + e;

Fs = 1/delta;
winlen = secDur*Fs;
detrendFlag = true;

tStartNew = dateshift(tStart,"start","minute");
tEndNew = dateshift(tEnd,"start","minute");
localStart = seconds(tStart - tStartNew);
new_npts = seconds(tEndNew-tStartNew)/delta;
localI = ceil(localStart/delta);
new_ref = delta*localI;
mainI = 1+localI;
tOrig = baseVec + localStart; %with local sampling rate

tQuery = baseVec + new_ref; %with local sampling rate
dinterp = interp1(tOrig,d,tQuery,"linear","extrap");

%%
d = zeros(new_npts+1,1);
d(mainI:mainI + npts - 1) = dinterp;
[dcut,startI] = cutWindows(d,winlen,nOverlap,detrendFlag);
t = tStartNew + seconds(delta*(0:new_npts)');
t = t(startI);

%%
[overtoneFreq,harmStrength,fundPow,...
    overtonePow,interharmPow,overtoneI,powerOrig,pxx,fxx] = tremometer(dcut,Fs,nOvertones);
harmStrength = sort(harmStrength,2,"descend");
harmStrength = harmStrength(:,1:3);
maxStrength = harmStrength(:,1); %max(harmStrength,[],2,"omitnan");
medStrength = harmStrength(:,2); %median(harmStrength,2,"omitnan");
minStrength = harmStrength(:,3); %min(harmStrength,[],2,"omitnan");
maxStrength = mean([medStrength maxStrength],2,"omitnan");

goodI = maxStrength >= MINMAX & ...
    medStrength >= MINMED & ...
    minStrength >= MINMIN & ...
    overtoneFreq(:,1) >= MINFREQ & ...
    sum(overtonePow >= MINPSD,2) > nOvertones & ...
    fundPow >= MINPOW & ...
    isfinite(maxStrength);

%% Make a table of harmonic tremor detections
detection_times = t(goodI);
fund_freq = overtoneFreq(goodI,1);
harm_freq_1 = overtoneFreq(goodI,2);
harm_freq_2 = overtoneFreq(goodI,3);
harm_freq_3 = overtoneFreq(goodI,4);
HSI_1 = harmStrength(goodI,1);
HSI_2 = harmStrength(goodI,2);
HSI_3 = harmStrength(goodI,3);
fundPow_ = fundPow(goodI);
fundPowOrig = fundPow;

%%
detections = table((1:sum(goodI))',detection_times,fund_freq,...
    harm_freq_1,harm_freq_2,harm_freq_3,HSI_1,HSI_2,HSI_3,fundPow_);
disp(detections);

%%
if ~plotFlag
    return;
end

figure('units','normalized','outerposition',[0 0.1 0.5 0.9]);
nSubplots = 3;
tiledlayout(nSubplots,1,"Padding","compact"); %,"TileSpacing","none");
ax = gobjects(nSubplots,1);

ax(1) = nexttile;
semilogy(ax(1),t,overtoneFreq(:,1),'p'); zoom on; grid on; hold on;
semilogy(ax(1),t,overtoneFreq(:,2:end),'.');
semilogy(ax(1),t(goodI),overtoneFreq(goodI,1),'ko'); zoom on; grid on; hold on;

ax(2) = nexttile;
semilogy(ax(2),t,harmStrength,'.'); zoom on; grid on; hold on;
semilogy(ax(2),t(goodI),maxStrength(goodI),'ko'); zoom on; grid on; hold on;

ax(3) = nexttile;
semilogy(ax(3),t,fundPow,'.'); zoom on; grid on; hold on;
linkaxes(ax,"x");

overtoneFreq(~goodI,:) = [];
harmStrength(~goodI,:) = [];
fundPow(~goodI) = [];