function tPotential = sangayRegionalUniqueEvents()
load('~/research/now/sangay/SangayRegionalAnalysis_v1','tabs','pks','templateIndex','maxAmpRMS','pksOrig','Neff');

%%
[tabs,sI] = sort(tabs);
maxAmpRMS = maxAmpRMS(sI,:);
Neff = Neff(sI);
pks = pks(sI);
templateIndex = templateIndex(sI);
pksOrig = pksOrig(sI);

%%
[tabs,pksOrig,removeIndices] = removeRepeatedMatches(tabs,pksOrig,60);
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    maxAmpRMS(rI,:) = [];
    Neff(rI) = [];
    templateIndex(rI) = [];
    pks(rI) = [];
end

%%
[tabs,sI] = sort(tabs);
maxAmpRMS = maxAmpRMS(sI,:);
Neff = Neff(sI);
pks = pks(sI);
templateIndex = templateIndex(sI);
pksOrig = pksOrig(sI);

diffSeconds = seconds(diff(tabs));
secThresh = 60;
badI = false(size(sI));
badI2 = find(diffSeconds <= secThresh);
badI3 = badI2+1;
badPKS1 = pksOrig(badI2);
badPKS2 = pksOrig(badI3);
secondGood = badPKS2 >= badPKS1;

badI2(secondGood) = badI3(secondGood);
badI(badI2) = true;

%%
badI = badI | Neff == 0 | (Neff == 1 & pksOrig <= 13) | ...
    (Neff == 2 & pksOrig <= 12.5) | ...
    (Neff == 3 & pksOrig <= 12) | ...
    (Neff == 4 & pksOrig <= 11.5) | ...
    (Neff == 5 & pksOrig <= 11) | ...
    (Neff == 6 & pksOrig <= 10.5) | ...
    (Neff == 7 & pksOrig <= 10) | ...
    (Neff == 8 & pksOrig <= 9.5) | ...
    (Neff == 9 & pksOrig <= 9) | ...
    (Neff == 10 & pksOrig <= 8.5) |...
    (Neff >= 11 & pksOrig <= 8);

%%
tPotential = sort(tabs(~badI));
