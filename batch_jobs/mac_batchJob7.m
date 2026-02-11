clear; close all;
cd ~/igdata/
load BREF_templates.mat
T2 = table2struct(T2);
T = reshape(T2,sizeT);
clear T2

lfc = 1;
hfc = 8;
newFs = 50;

dayStart = datetime(2022,10,25);
dayEnd = datetime(2022,10,25);

dayVec = (dayEnd:-1:dayStart)';
lDays = length(dayVec);

for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,1,"BREF",["BHZ";"BHN";"BHE"],"EC");
    if isnat(S(1).ref)
        continue;
    end

    Sf = syncWaveforms(filterWaveforms(detrendWaveforms((S)),lfc,hfc));
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain] = ...
        templateSearch(Sf,T,false);

    [tMain,sI] = sort(tMain);
    dMag = dMag(sI);
    ampMain = ampMain(sI);
    ccMain = ccMain(sI);
    evidMain = evidMain(sI);
    templateNumber = templateNumber(sI);
    magMain = magMain(sI);
    madMain = madMain(sI);

    fname = strcat('~/Desktop/dayTemplateSearch_',datestr(dateshift(datetime(day_),'start','day'),'yyyy.mm.dd'),'.mat');
    save(fname,"madMain","templateNumber","magMain","evidMain","ccMain","ampMain","dMag","tMain");
end
