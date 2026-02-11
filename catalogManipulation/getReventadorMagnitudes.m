function [t,M,magErr,Mcorr,amps] = getReventadorMagnitudes()
load('~/igdata/ReventadorSubspaceDetectorResults_v10.mat','z2p','tabs');
%load ~/research/now/reventador/ReventadorSubspaceDetectorResults_v10.mat

[stla,stlo] = metaDataFromStationList(["CASC";"BONI";"ANTS";"ANTG"]);
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(stla,stlo,-0.080850,-77.657995,refEllipse)*1e-3;

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

[tabs3,~] = filterCatalog(tabs2,CASC4,30);
t1_commonI = synchronizeCatalog(tabs2,tabs3,3,true,false);
CASC4 = CASC4(t1_commonI);
BONI4 = BONI4(t1_commonI);
ANTS4 = ANTS4(t1_commonI);
ANTG4 = ANTG4(t1_commonI);

load('~/igdata/CASC_BONI_ANTS_ANTG_REVS_ampsCorrected.mat','tsample','ampsCorrected');
%load ~/research/now/reventador/CASC_BONI_ANTS_ANTG_REVS_ampsCorrected.mat
[t1_commonI,~,~,~] = synchronizeCatalog(tabs3,tsample,3,true,false);

%
%close all;
%figure(); loglog(CASC4(t1_commonI),ampsCorrected(:,1),'.'); zoom on; grid on;
bCASC = flipud(robustfit(log10(CASC4(t1_commonI)),log10(ampsCorrected(:,1))));
yq = polyval(bCASC,sort(log10(CASC4(t1_commonI))));
%hold on;
%ll = loglog(10.^sort(log10(CASC4(t1_commonI))),10.^yq,'k','linewidth',4);
%ll.Color(4) = 0.5;

%figure(); loglog(BONI4(t1_commonI),ampsCorrected(:,3),'.'); zoom on; grid on;
bBONI = flipud(robustfit(log10(BONI4(t1_commonI)),log10(ampsCorrected(:,3))));
yq = polyval(bBONI,sort(log10(BONI4(t1_commonI))));
%hold on;
%ll = loglog(10.^sort(log10(BONI4(t1_commonI))),10.^yq,'k','linewidth',4);
%ll.Color(4) = 0.5;

%figure(); loglog(ANTS4(t1_commonI),ampsCorrected(:,5),'.'); zoom on; grid on;
bANTS = flipud(robustfit(log10(ANTS4(t1_commonI)),log10(ampsCorrected(:,5))));
yq = polyval(bANTS,sort(log10(ANTS4(t1_commonI))));
%hold on;
%ll = loglog(10.^sort(log10(ANTS4(t1_commonI))),10.^yq,'k','linewidth',4);
%ll.Color(4) = 0.5;

%figure(); loglog(ANTG4(t1_commonI),ampsCorrected(:,7),'.'); zoom on; grid on;
bANTG = flipud(robustfit(log10(ANTG4(t1_commonI)),log10(ampsCorrected(:,7))));
yq = polyval(bANTG,sort(log10(ANTG4(t1_commonI))));
%hold on;
%ll = loglog(10.^sort(log10(ANTG4(t1_commonI))),10.^yq,'k','linewidth',4);
%ll.Color(4) = 0.5;


% convert everything to Wood-Anderson
amps(:,1) = 10.^polyval(bCASC,log10(amps(:,1)));
amps(:,2) = 10.^polyval(bBONI,log10(amps(:,2)));
amps(:,3) = 10.^polyval(bANTS,log10(amps(:,3)));
amps(:,4) = 10.^polyval(bANTG,log10(amps(:,4)));

%%
load('~/research/now/reventador/ReventadorRegionalAttenuationModel_v02AUG2023.mat',...
    'mbest','d_','gamma100','fixFlag','n');
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

M = median(Mcorr,2,"omitnan");
t = tabs;
magErr = mad(Mcorr,1,2);
