function [tPotential,z2p,NCC,Neff,p2rms,kurt] = filterUniqueEvents(longtermCatalogName,MPD)
%,pivotPoint)
if nargin < 1
    longtermCatalogName = []; %'~/research/now/sangay/SangayRegionalAnalysis_v5';
end

if nargin < 2
    MPD = [];
end

% if nargin < 3
%     pivotPoint = 15; %15 is too "loose" and 18 is probably too strict
% end

%%
if isempty(longtermCatalogName)
    longtermCatalogName = '~/research/now/sangay/SangayRegionalAnalysis_v6';
end

if isempty(MPD)
    MPD = 30;
end

load(longtermCatalogName,'tabs','p2rms','z2p','NCC','Neff','kurt');

%%
[tabs,sI] = sort(tabs);
z2p = z2p(sI,:);
Neff = Neff(sI);
p2rms = p2rms(sI,:);
%templateIndex = templateIndex(sI);
NCC = NCC(sI);
kurt = kurt(sI,:);

%%
[tabs,NCC,removeIndices] = removeRepeatedMatches(tabs,NCC,MPD);
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    z2p(rI,:) = [];
    Neff(rI) = [];
    %templateIndex(rI) = [];
    p2rms(rI,:) = [];
    kurt(rI,:) = [];
end

%%
[tabs,sI] = sort(tabs);
z2p = z2p(sI,:);
Neff = Neff(sI);
p2rms = p2rms(sI,:);
%templateIndex = templateIndex(sI);
NCC = NCC(sI);
kurt = kurt(sI,:);

%%
diffSeconds = seconds(diff(tabs));
secThresh = MPD; %/2;
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
% lTemplates = length(unique(templateIndex));
% lN = length(unique(Neff));
% for i = 1:lTemplates
%     for j = 1:lN
%         badI = badI | (templateIndex == i & Neff == j & p2rms < (pivotPoint - (i - 1)/4 - (j - 1)));
%     end
% end

%%
goodI = sort(find(~badI));
tPotential = tabs(goodI);
z2p = z2p(goodI,:);
Neff = Neff(goodI);
p2rms = p2rms(goodI,:);
%templateIndex = templateIndex(goodI);
NCC = NCC(goodI);
kurt = kurt(goodI,:);
