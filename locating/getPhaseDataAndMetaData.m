function [observedTT,stla,stlo,stel,SFlag] = getPhaseDataAndMetaData(E,SFlag)
Pphases = E.Pphases;

%%
stnm = pull(Pphases,'stnm');
[stla,stlo,stel] = metaDataFromStationList(stnm);
goodstnm = isfinite(stla);
nPphases = sum(goodstnm);
stla = stla(goodstnm);
stlo = stlo(goodstnm);
stel = stel(goodstnm);
stnm = stnm(goodstnm);

%%
observedTT = pull(Pphases(goodstnm),'t');

%% get only stations that i have metadata for
if ~SFlag
    return;
end

%%
Sphases = E.Sphases;
nSphases = E.nSphases;

if nSphases < 2
    SFlag = false;
    return;
end

nS = 0;
stnmS = pull(Sphases,'stnm');
observedTT = [observedTT NaT(nPphases,1)];
for i = 1:nSphases
    stnmS_ = stnmS(i);
    [lia,locb] = ismember(stnmS_,stnm);
    if ~lia
        continue;
    end
    locb = locb(lia);
    observedTT(locb,2) = pull(Sphases(i),'t');
    nS = nS+1;
end

%%
if nSphases < 2
    SFlag = false;
    observedTT = observedTT(:,1);
    return;
end
