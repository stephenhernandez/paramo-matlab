%function esmeraldas_2025_detector(dayStart,dayEnd,dayInc)
clear; close all;
dayStart = datetime(2025,04,20);
dayEnd = datetime(2025,04,30);
dayInc = 1;

template_function_directory = fullfile("~","masa","template_search","template_functions");
load(fullfile(template_function_directory,"esmeraldas_2025_v1"),"T");

%%
Torig = T;

lfc = 2;
hfc = 8;
newFs = 50;
thresh = 7;

dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

writeFlag = true;
plotFlag = false;
verboseFlag = true;
minNpts = max(pull(Torig(:),'npts'));

%%
for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        ["AES1";"SNLR";"ESM1";"PAC1"],...
        ["HNZ";"HHZ";"BHZ";"BLZ";"HNN";"HHN";"BHN";"BLN";"HNE";"HHE";"BHE";"BLE"],"EC");
    kstnms = pull(S,'kstnm');
    kcmpnms = pull(S,'kcmpnm');

    %%
    badI = isnat(pull(S,'ref')) | max(pull(S,'npts')) < minNpts;
    S(badI) = [];
    lS = length(S);
    if ~lS %min number of traces is 1
        fprintf("not enough data for day: %s\n",string(day_));
        continue;
    end

    %%
    Sf = padWaveforms(filterWaveforms(detrendWaveforms(nanGapWaveforms(detrendWaveforms(S),0)),lfc,hfc));
    Sf = resampleWaveforms(Sf,newFs);

    customPrefix = "esmeraldas2025";
    Tsncls = strcat(pull(T(:,1),'knetwk'),pull(T(:,1),'kstnm'),pull(T(:,1),'khole'),pull(T(:,1),'kcmpnm'));
    Ssncls = strcat(pull(Sf,'knetwk'),pull(Sf,'kstnm'),pull(Sf,'khole'),pull(Sf,'kcmpnm'));

    [lia,locb] = ismember(Ssncls,Tsncls);
    if sum(lia)
        [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain,ccnorm,tLong] = ...
            templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,verboseFlag,thresh,customPrefix);
    end
end
