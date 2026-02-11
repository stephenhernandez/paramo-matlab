clear; close all;
cd ~/research/now/fernandina/
load FER1_templates.mat
T2 = table2struct(T2);
T = reshape(T2,sizeT);
clear T2

lfc = 4;
hfc = 16;
newFs = 50;

dayStart = datetime(2022,11,20);
dayEnd = datetime(2022,11,22);

dayVec = (dayEnd:-1:dayStart)';
lDays = length(dayVec);

for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,1,"FER1",["BHZ";"BHN";"BHE"],"EC");
    if isnat(S(1).ref)
        continue;
    end

    Sf = syncWaveforms(filterWaveforms(detrendWaveforms((S)),lfc,hfc));
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain] = ...
        templateSearch(Sf,T);
end