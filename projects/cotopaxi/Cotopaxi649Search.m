clear; close all;
cd ~/igdata

lfc = 0.7;
hfc = 7;
newFs = 50;
thresh = 8;

dayStart = datetime(2024,03,12);
dayEnd = datetime(2024,04,09);
dayInc = 1;

if dayEnd < dayStart
    fprintf('dayEnd should be dayStart, fix!\n');
    return;
end
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

writeFlag = true;
plotFlag = false;
verboseFlag = true;

load('Cotopaxi649Events.mat');
Torig = C';
T = Torig;
clear C

%%
minNpts = max(pull(Torig,'npts'));
for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        ["CO1V";"BREF";"BTAM";"BVC2"],...
        ["HHZ";"BHZ"],"EC");

    %%
    badI = isnat(pull(S,'ref')) | max(pull(S,'npts') < minNpts);
    S(badI) = [];
    lS = length(S);
    if ~lS
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end

    Sf = padWaveforms(filterWaveforms(detrendWaveforms((S)),lfc,hfc));
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    T = Torig(:,end);
    Tchan = char(pull(T,'kcmpnm'));
    Tchan = string(Tchan(:,2:3));
    Schan = char(pull(S,'kcmpnm'));
    Schan = string(Schan(:,2:3));
    Tsncls = strcat(pull(T,'knetwk'),pull(T,'kstnm'),pull(T,'khole'),Tchan);
    Ssncls = strcat(pull(Sf,'knetwk'),pull(Sf,'kstnm'),pull(Sf,'khole'),Schan);

    [lia,locb] = ismember(Ssncls,Tsncls);
    if sum(lia)
        %this section only works when the template lengths in npts is the
        %same from multiplet template to multiplet template...
        [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain] = ...
            templateSearch(Sf(lia),Torig(locb(lia),:),writeFlag,plotFlag,...
            verboseFlag,thresh,'Cotopaxi649DiscreteEventsMarch2024');
    end
end
