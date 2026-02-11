clear; close all; clc;

%% TODO:
% -do binary classifier based on amplitude ratios (of the 12 components, compare
% each component with every other component)
%
% -use max instead of median for attenuation inversion and magnitude
% calulation
%
% -experiment with models to synthesize some amplitudes based on other
% stations
%
% -for subspace detectors in general, in the future you should be saving
% the weights from the SVD inversion since they can be used in a
% classification tree binary detector
%
% -use narrower band filter (0.6-1.2 Hz.)
%
% -far in the future, re-run detector but with response removed before
% processing
%

tic;
load ~/igdata/ReventadorSubspaceDetectorResults_v10.mat
%load ~/research/now/reventador/ReventadorSubspaceDetectorResults_v10.mat

[stla,stlo] = metaDataFromStationList(["CASC";"BONI";"ANTS";"ANTG"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse)*1e-3;

maxFlag = true;
if maxFlag
    CASC = 1e9*max(z2p(:,1:3),[],2,"omitnan")/3.141950e+08;

    bI = tabs >= datetime(2020,01,338);
    BONI = 1e9*max(z2p(:,4:6),[],2,"omitnan");
    BONI(bI) = BONI(bI)/5.037350e+08;
    BONI(~bI) = BONI(~bI)/2.014940e+09;

    bI = tabs >= datetime(2014,01,226);
    ANTS = 1e9*max(z2p(:,7:9),[],2,"omitnan");
    ANTS(bI) = ANTS(bI)/5.037350e+08;
    ANTS(~bI) = ANTS(~bI)/3.141950e+08;

    bI = tabs >= datetime(2017,01,105);
    ANTG = 1e9*max(z2p(:,10:12),[],2,"omitnan");
    ANTG(bI) = ANTG(bI)/3.141950e+08;
    ANTG(~bI) = ANTG(~bI)/3.141950e+08;
else
    CASC = 1e9*median(z2p(:,1:3),2,"omitnan")/3.141950e+08;

    bI = tabs >= datetime(2020,01,338);
    BONI = 1e9*median(z2p(:,4:6),2,"omitnan");
    BONI(bI) = BONI(bI)/5.037350e+08;
    BONI(~bI) = BONI(~bI)/2.014940e+09;

    bI = tabs >= datetime(2014,01,226);
    ANTS = 1e9*median(z2p(:,7:9),2,"omitnan");
    ANTS(bI) = ANTS(bI)/5.037350e+08;
    ANTS(~bI) = ANTS(~bI)/3.141950e+08;

    bI = tabs >= datetime(2017,01,105);
    ANTG = 1e9*median(z2p(:,10:12),2,"omitnan");
    ANTG(bI) = ANTG(bI)/3.141950e+08;
    ANTG(~bI) = ANTG(~bI)/3.141950e+08;
end

tOrig = tabs;

CASC(CASC == 0) = NaN;
BONI(BONI == 0) = NaN;
ANTS(ANTS == 0) = NaN;
ANTG(ANTG == 0) = NaN;

amps = [CASC BONI ANTS ANTG];
fourI = sum(isfinite(amps),2) == 4;
fourI = fourI & CASC >= 100;
fourI = fourI & BONI >= 40;
fourI = fourI & ANTS >= 80;
fourI = fourI & ANTG >= 40;

fourI = find(fourI);
tabs2 = tabs(fourI);

CASC4 = CASC(fourI);
BONI4 = BONI(fourI);
ANTS4 = ANTS(fourI);
ANTG4 = ANTG(fourI);

[tabs3,z2p3] = filterCatalog(tabs2,CASC4,30);
t1_commonI = synchronizeCatalog(tabs2,tabs3,3,true,false);
CASC4 = CASC4(t1_commonI);
BONI4 = BONI4(t1_commonI);
ANTS4 = ANTS4(t1_commonI);
ANTG4 = ANTG4(t1_commonI);

load('~/igdata/CASC_BONI_ANTS_ANTG_REVS_ampsCorrected.mat','tsample','ampsCorrected');
%load ~/research/now/reventador/CASC_BONI_ANTS_ANTG_REVS_ampsCorrected.mat
[t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(tabs3,tsample,3,true,false);

%
close all;
figure(); loglog(CASC4(t1_commonI),ampsCorrected(:,1),'.'); zoom on; grid on;
bCASC = flipud(robustfit(log10(CASC4(t1_commonI)),log10(ampsCorrected(:,1))));
yq = polyval(bCASC,sort(log10(CASC4(t1_commonI))));
hold on;
ll = loglog(10.^sort(log10(CASC4(t1_commonI))),10.^yq,'k','linewidth',4);
ll.Color(4) = 0.5;

figure(); loglog(BONI4(t1_commonI),ampsCorrected(:,3),'.'); zoom on; grid on;
bBONI = flipud(robustfit(log10(BONI4(t1_commonI)),log10(ampsCorrected(:,3))));
yq = polyval(bBONI,sort(log10(BONI4(t1_commonI))));
hold on;
ll = loglog(10.^sort(log10(BONI4(t1_commonI))),10.^yq,'k','linewidth',4);
ll.Color(4) = 0.5;

figure(); loglog(ANTS4(t1_commonI),ampsCorrected(:,5),'.'); zoom on; grid on;
bANTS = flipud(robustfit(log10(ANTS4(t1_commonI)),log10(ampsCorrected(:,5))));
yq = polyval(bANTS,sort(log10(ANTS4(t1_commonI))));
hold on;
ll = loglog(10.^sort(log10(ANTS4(t1_commonI))),10.^yq,'k','linewidth',4);
ll.Color(4) = 0.5;

figure(); loglog(ANTG4(t1_commonI),ampsCorrected(:,7),'.'); zoom on; grid on;
bANTG = flipud(robustfit(log10(ANTG4(t1_commonI)),log10(ampsCorrected(:,7))));
yq = polyval(bANTG,sort(log10(ANTG4(t1_commonI))));
hold on;
ll = loglog(10.^sort(log10(ANTG4(t1_commonI))),10.^yq,'k','linewidth',4);
ll.Color(4) = 0.5;

%
CASC_WA = 10.^polyval(bCASC,log10(CASC4));
BONI_WA = 10.^polyval(bBONI,log10(BONI4));
ANTS_WA = 10.^polyval(bANTS,log10(ANTS4));
ANTG_WA = 10.^polyval(bANTG,log10(ANTG4));

% convert everything to Wood-Anderson
amps(:,1) = 10.^polyval(bCASC,log10(amps(:,1)));
amps(:,2) = 10.^polyval(bBONI,log10(amps(:,2)));
amps(:,3) = 10.^polyval(bANTS,log10(amps(:,3)));
amps(:,4) = 10.^polyval(bANTG,log10(amps(:,4)));
toc;

%%
close all;
nStations = 4;
G_ = full(Gvdcc(nStations));
G_ = G_(1:end-1,:);

d = [CASC_WA BONI_WA ANTS_WA ANTG_WA];
maxNumberOfAmpsForInversion = 5e5;
findI = (1:length(CASC_WA))';
sampleIndex = sort(randsample(length(CASC_WA),min([maxNumberOfAmpsForInversion length(CASC_WA)])));
ampsScrambled = d(sampleIndex,:);

mbest = [];
rdum = logspace(log10(5),log10(700),501)';

fixFlag = true;
if fixFlag
    %n = 1/2; % median mad: 0.0176244, Q = 30, Vs = 2.3;
    n = 1; % median mad: 0.0176993, Q=60, Vs = 2.3;
    G_ = [getDD(d_) G_];
else
    G_ = [getDD(log10(d_)) getDD(d_) G_];
end

tic;
lKstnms = nStations;
ampsScrambled2 = ampsScrambled;
sizeSamples = size(ampsScrambled,1);
lCombos = 0.5*lKstnms*(lKstnms-1);
nn = 1;
nnMax = sizeSamples*lCombos;
dd = zeros(nnMax,1);
G = NaN(nnMax,lKstnms+1);
for i = 1:sizeSamples
    disp(i);
    amps_ = ampsScrambled(i,:);
    amps_ = amps_./median(amps_);

    G(nn:nn+lCombos-1,:) = G_;
    if fixFlag
        %dd = [dd; -getDD(log10(amps_')) - n*getDD(log10(d_))];
        dd(nn:nn+lCombos-1) = -getDD(log10(amps_')) - n*getDD(log10(d_));
    else
        dd = [dd; -getDD(log10(amps_'))];
    end

    %
    ampsScrambled2(i,:) = amps_;
    nn = nn+lCombos;
end
toc;

if fixFlag
    G = [G; 0 ones(1,nStations)];
else
    G = [G; 0 0 ones(1,nStations)];
end
dd = [dd; 0];

mbest_ = pinv(G)*dd;
mbest = [mbest mbest_];

format long g
disp(mbest);

if fixFlag
    att1 = n*log10(rdum) + mbest_(1)*rdum;              %hernandez
else
    att1 = mbest_(1)*log10(rdum) + mbest_(2)*rdum;      %hernandez
end

att2 = 1.11*log10(rdum) + 0.00189*rdum + 0.591;     %uhrhammer
gamma1 = -interp1(rdum,att1,1);
gamma17 = -interp1(rdum,att1,17) + 2;
gamma100 = -interp1(rdum,att1,100) + 3;
attGamma1 = att1 + gamma17;
attGamma100 = att1 + gamma100;

%%
plotFlag = true;
if plotFlag
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogx(rdum,att1,'.'); zoom on; grid on; hold on;
    semilogx(rdum,att2,'.');
    ylim([0 6]);

    hold on;
    semilogx(rdum,attGamma1,'.','color',[0.5 0.5 0.5]); grid on;
    hold on;
    semilogx(rdum,attGamma100,'.','color',[0.3 0.3 0.3]); grid on;
    legend('new','uhrhammer','fix17','fix100','location','northwest');
end

%%
close all;
McorrOrig = amps;
McorrOrig(~isfinite(McorrOrig) | McorrOrig == 0) = NaN;
Mcorr = McorrOrig;
for i = 1:4
    if fixFlag
        Mcorr(:,i) = log10(Mcorr(:,i)) + (n*log10(d_(i)) + mbest(1)*d_(i) + gamma100 + mbest(i+1));
    else
        Mcorr(:,i) = log10(Mcorr(:,i)) + (mbest(1)*log10(d_(i)) + mbest(2)*d_(i) + gamma100 + mbest(i+2));
    end
end
fprintf('mean mad: %g\n',mean(mad(Mcorr,1,2),"omitnan"));

figure('units','normalized','outerposition',[0 0 1 1]);
M1 = median(Mcorr,2,"omitnan");
t = tOrig;
magErr = mad(Mcorr,1,2);
minMag = 0.2;
maxMag = 10;
mI = M1 >= minMag & magErr <= 1 & M1 <= maxMag;
plot(t(mI),M1(mI),'.'); zoom on; grid on; hold on;
title(sprintf("Reventador Estimated Pseudo-Magnitudes, Min. Mag: %f, Max. Mag.: %f\n",minMag,maxMag));

%
nDays = 14;
zeroPhaseFlag = ~true;
Nmed = 1;
plot(t(mI),medfiltSH(M1(mI),Nmed,zeroPhaseFlag),'.'); zoom on;
legend('orig',sprintf("%d-point running median\n",Nmed),'location','best');

[rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = ...
    t2r(t(mI),days(nDays),medfiltSH(M1(mI),Nmed,zeroPhaseFlag),~true);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(t(mI),rate/nDays,'.'); zoom on; grid on;
title(sprintf("Rate Per Day, Smoothed over: %d day(s)\n",nDays));

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t(mI),meanMagsFixedTimeWin,'.'); zoom on; grid on;
title("Mean Mag (Fixed Time Win)");

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t(mI),sumEnergyFixedTimeWin/nDays,'.'); zoom on; grid on;
title("Average Daily Energy (Fixed Time Win)");
ylabel('$\propto$ Energy');

figure('units','normalized','outerposition',[0 0 1 1]);
plot(t(mI),medianMagsFixedTimeWin,'.'); zoom on; grid on;
title("Median Mag (Fixed Time Win)");

figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(t(mI),(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); zoom on; grid on;
title("Projected Daily Energy (From Median Mag.)");
ylabel('$\propto$ Energy');
