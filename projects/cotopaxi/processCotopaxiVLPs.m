%function processCotopaxiVLPs
clear; close all;
tMain = [];
ccMain = tMain;
ampMain = tMain;
dMag = tMain;
magMain = tMain;
templateNumber = tMain;
madMain = tMain;
evidMain = tMain;
nUsedMain = tMain;

for i = 1:43
    customPrefix = sprintf("CotopaxiVLP_%02d",i);
    [tMain_,ccMain_,ampMain_,dMag_,magMain_,templateNumber_,madMain_,evidMain_,nUsedMain_] = ...
        loadRepeaterCatalog(fullfile(customPrefix));
    tMain = [tMain; tMain_];
    ccMain = [ccMain; ccMain_];
    ampMain = [ampMain; ampMain_];
    dMag = [dMag; dMag_];
    magMain = [magMain; magMain_];
    templateNumber = [templateNumber; templateNumber_*i];
    madMain = [madMain; madMain_];
    evidMain = [evidMain; evidMain_];
    nUsedMain = [nUsedMain; nUsedMain_];
end

[tMain,sortI] = sort(tMain);
ccMain = ccMain(sortI);
ampMain = ampMain(sortI);
dMag = dMag(sortI);
magMain = magMain(sortI);
templateNumber = templateNumber(sortI);
madMain = madMain(sortI);
evidMain = evidMain(sortI);
nUsedMain = nUsedMain(sortI);

%%
close all; 
[tMain,ccMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain] = ...
    filterCatalog(tMain,ccMain,45,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain);

ccI = ccMain >= 0.45 & madMain >= 10 & ampMain >= 2e1 & nUsedMain >= 6 & ...
    tMain >= datetime(2016,09,01); % & templateNumber == 42;
figure(); semilogy(tMain(ccI),ampMain(ccI),'.'); zoom on; grid on;
figure(); semilogy(tMain(ccI),ccMain(ccI),'.'); zoom on; grid on;
figure(); semilogy(tMain(ccI),madMain(ccI),'.'); zoom on; grid on;
figure(); plot(tMain(ccI),(1:sum(ccI))','.'); zoom on; grid on;
cd ~/research/now/cotopaxi/;
