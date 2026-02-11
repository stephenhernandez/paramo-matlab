function [tPotential,z2p,NCC,Neff,templateIndex,scaledCC,kurt] = sangayRegionalUniqueEvents_v4()
load('~/research/now/sangay/SangayRegionalAnalysis_v4','tabs','scaledCC','templateIndex','z2p','NCC','Neff','kurt');

%%
[tabs,sI] = sort(tabs);
z2p = z2p(sI,:);
Neff = Neff(sI);
scaledCC = scaledCC(sI);
templateIndex = templateIndex(sI);
NCC = NCC(sI);
kurt = kurt(sI);

%%
[tabs,NCC,removeIndices] = removeRepeatedMatches(tabs,NCC,60);
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    z2p(rI,:) = [];
    Neff(rI) = [];
    templateIndex(rI) = [];
    scaledCC(rI) = [];
    kurt(rI) = [];
end

%%
[tabs,sI] = sort(tabs);
z2p = z2p(sI,:);
Neff = Neff(sI);
scaledCC = scaledCC(sI);
templateIndex = templateIndex(sI);
NCC = NCC(sI);
kurt = kurt(sI);

%%
diffSeconds = seconds(diff(tabs));
secThresh = 15;
badI = false(size(sI));
badI2 = find(diffSeconds <= secThresh);
badI3 = badI2+1;
badPKS1 = NCC(badI2);
badPKS2 = NCC(badI3);
secondGood = badPKS2 >= badPKS1;

%%
badI2(secondGood) = badI3(secondGood);
badI(badI2) = true;

%%
lTemplates = length(unique(templateIndex));
lN = length(unique(Neff));
pivotPoint = 17; %16 is too "loose" and 18 is probably too strict
for i = 1:lTemplates
    for j = 1:lN
        badI = badI | (templateIndex == i & Neff == j & scaledCC < (pivotPoint + (i - 1)/4 - (j - 1)));
    end
end

%%
goodI = sort(find(~badI));
tPotential = tabs(goodI);
z2p = z2p(goodI,:);
Neff = Neff(goodI);
scaledCC = scaledCC(goodI);
templateIndex = templateIndex(goodI);
NCC = NCC(goodI);
kurt = kurt(goodI);
