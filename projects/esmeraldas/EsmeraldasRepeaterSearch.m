clear; close all;
%cd ~/research/now/esmeraldas/
cd ~/igdata

%load("EsmeraldasTemplates.mat");
load("EsmeraldasTemplates_2.mat");
Torig = T;

%%
lfc = 2;
hfc = 12; %-inf;
newFs = 40; %20;
thresh = 8; %10;

%dayStart = datetime(2022,05,06);
%dayEnd = datetime(2022,05,06);
dayStart = datetime(2021,12,01);
dayEnd = datetime(2022,01,01);
%dayStart = datetime(2022,01,01);
%dayEnd = datetime(2022,05,12);
dayInc = 1;

kstnms = ["HB15";"HB18";"HB12";"HB17";"HB19";"PTGL";...
    "HB15";"HB18";"HB12";"HB17";"HB19";"PTGL";...
    "HB15";"HB18";"HB12";"HB17";"HB19";"PTGL"];

knetwk_ = [repmat("XF",5,1);"EC"];
knetwks = repmat(knetwk_,3,1);
kcmpnms = [repmat("HHZ",6,1); repmat("HHN",6,1); repmat("HHE",6,1)];
kholes = repmat("",18,1);
lKstnms = length(kstnms);

dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

%%
writeFlag = true;
plotFlag = false;
verboseFlag = true;

%%
for i = 1:lDays
    day_ = dayVec(i);
    S = populateWaveforms(lKstnms);
    for j = 1:lKstnms
        kstnm_ = kstnms(j);
        kcmpnm_ = kcmpnms(j);
        khole_ = kholes(j);
        knetwk_ = knetwks(j);

        S(j) = loadWaveforms(day_,dayInc,...
            kstnm_,kcmpnm_,knetwk_,khole_);
    end

    badI = isnat(pull(S,'ref'));
    S(badI) = [];
    lS = length(S);
    if ~lS
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end

    Sf = padWaveforms(filterWaveforms(detrendWaveforms((S)),lfc,hfc));
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    customPrefix = "esmeraldas";
    Tsncls = strcat(pull(T,'knetwk'),pull(T,'kstnm'),pull(T,'khole'),pull(T,'kcmpnm'));
    Ssncls = strcat(pull(Sf,'knetwk'),pull(Sf,'kstnm'),pull(Sf,'khole'),pull(Sf,'kcmpnm'));
    [lia,locb] = ismember(Ssncls,Tsncls);
    if sum(lia)
        [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain] = ...
            templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,verboseFlag,thresh,customPrefix);
    end
end
