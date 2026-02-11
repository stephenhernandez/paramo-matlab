function [tPotential,maxAmpRMS,pksOrig,Neff,templateIndex,pks] = sangayRegionalUniqueEvents_v3()
load('~/research/now/sangay/SangayRegionalAnalysis_v3','tabs','pks','templateIndex','maxAmpRMS','pksOrig','Neff');

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

%% loop through and apply correction to magnitudes
wd = pwd;
cd ~/response/
load('~/research/now/sangay/nineTemplatesEighteenSensorsSangay','kstnm','chan');

%%
[stla,stlo] = metaDataFromStationList(kstnm);
refEllipse = referenceEllipsoid('wgs84');
r = distance(stla,stlo,-2.0051,-78.3415,refEllipse)*1e-3;

%%
for i = 1:length(kstnm)
    chan_ = char(chan(i));
    kstnm_ = char(kstnm(i));
    
    dirName = strcat(kstnm_,'_EC');
    if exist(dirName,'dir')
        cd(dirName);
        searchPattern = strcat('SAC_PZs_EC_',kstnm_,'_',chan_,'*');
        pzFiles = dir(searchPattern);
        lPZs = length(pzFiles);
        disp([kstnm_,chan_,':',num2str(lPZs)]);
        for k = 1:lPZs
            pz_ = pzFiles(k).name;
            [~,~,~,~,~,~,~,tStart_,tEnd_,~,~,~,~,...
                ~,~,~,~,~,sensitivity_] = read_sac_pole_zero(pz_);
            if sensitivity_ > 1
                tI = tabs >= tStart_ & tabs < tEnd_;
                if sum(tI)
                    amps_ = maxAmpRMS(:,i);
                    tI = tI & isfinite(amps_) & amps_ > 0;
                    %amps_(tI) = log10(amps_(tI)/sensitivity_);
                    amps_(tI) = log10(amps_(tI)/sensitivity_) + 1.11*log10(r(i)) + 0.00189*r(i) + 4.5;
                    maxAmpRMS(:,i) = amps_;
                end
            end
        end
    end
    cd ..;
end
cd(wd);

%%
%ampMED = nanmedian(maxAmpRMS,2);
%ampSTD = mad(maxAmpRMS,1,2);

%%
diffSeconds = seconds(diff(tabs));
secThresh = 30;
badI = false(size(sI));
badI2 = find(diffSeconds <= secThresh);
badI3 = badI2+1;
badPKS1 = pksOrig(badI2);
badPKS2 = pksOrig(badI3);
secondGood = badPKS2 >= badPKS1;

badI2(secondGood) = badI3(secondGood);
badI(badI2) = true;

%%
uniqueTemplates = unique(templateIndex);
% for i = 1:length(uniqueTemplates)
badI = badI | Neff == 0 | (Neff == 1 & pksOrig < 10) | ...
    (Neff == 2 & pksOrig < 9.5) | ...
    (Neff == 3 & pksOrig < 9) | ...
    (Neff == 4 & pksOrig < 8.5) | ...
    (Neff == 5 & pksOrig < 8) |...
    (Neff == 6 & pksOrig < 7.6) | ...
    (Neff == 7 & pksOrig < 7) | ...
    (Neff == 8 & pksOrig < 6.5) | ...
    (Neff == 9 & pksOrig < 6) | ...
    (Neff >= 10 & pksOrig < 5.5);% | ...
%end
%(Neff >= 11 & pksOrig < 4);

%%
goodI = sort(find(~badI));
tPotential = tabs(goodI);
maxAmpRMS = maxAmpRMS(goodI,:);
Neff = Neff(goodI);
pks = pks(goodI);
templateIndex = templateIndex(goodI);
pksOrig = pksOrig(goodI);