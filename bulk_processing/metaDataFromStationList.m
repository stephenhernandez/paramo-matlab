function [stlat,stlon,stelev] = metaDataFromStationList(stnm,ntwk,chan,locID)
if nargin < 2; ntwk = []; end
if nargin < 3; chan = []; end
if nargin < 4; locID = []; end

%%
try
    load('~/igdata/ecuadorSensorDataTable10','kstnm','knetwk','kcmpnm','khole','stla','stlo','stel');
catch
    try
        load('~/Desktop/igdata/ecuadorSensorDataTable10','kstnm','knetwk','kcmpnm','khole','stla','stlo','stel');
    catch ME
        rethrow(ME)
    end
end

lS = length(stnm);
lN = length(ntwk);
lC = length(chan);
lL = length(locID);

%%
stlat = NaN(lS,1);
stlon = stlat;
stelev = stlat;

%%
allFour = lS == lN && lS == lC && lS == lL;
if ~allFour
    for i = 1:lS
        stnm_ = stnm(i);
        lia = ismember(kstnm,stnm_);
        if sum(lia)
            locs = find(lia);
            stlat(i,1) = stla(locs(end));
            stlon(i,1) = stlo(locs(end));
            stelev(i,1) = stel(locs(end));
        else
            fprintf('Station %s not found.\n',stnm_);
        end
    end
    return;
end

%%
mySNCLs = strcat(ntwk,stnm,locID,chan);
allSNCLs = strcat(knetwk,kstnm,khole,kcmpnm);

[lia,locb] = ismember(mySNCLs,allSNCLs);
if ~sum(lia)
    return;
end
locb = locb(lia);
stlat(lia) = stla(locb);
stlon(lia) = stlo(locb);
stelev(lia) = stel(locb);
