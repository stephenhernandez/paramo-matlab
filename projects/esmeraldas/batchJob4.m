clear; close all; clc;
[~,~,t,eqlat,eqlon,eqdepth,eqmag,id] = readCat1(datetime(2022,01,01),datetime(2022,05,12),1,6.5);

%
fzqpI = find(eqmag == max(eqmag));
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(eqlat,eqlon,eqlat(fzqpI),eqlon(fzqpI),refEllipse);
d_ = d_ * 1e-3;
lI = d_ <= 50;
ids = id(lI);
for i = 1:sum(lI)
    E(i,1) = readSCBulletin(ids(i));
    disp(ids(i));
end

%
n = 0;
tnew = NaT(length(E),1);
tMain = tnew;
idNew = repmat("",length(E),1);
for i = 1:length(E)
    Pphases = E(i).Pphases;
    kstnms = pull(Pphases,'stnm');
    [lia,locb] = ismember("HB15",kstnms);
    locb = locb(lia);
    if lia
        disp(locb);
        n = n+1;
        tMain(n) = E(i).t;
        tnew(n) = Pphases(locb).t;
        idNew(n) = ids(i);
    end
end
tnew = tnew(1:n);
idNew = idNew(1:n);
tMain = tMain(1:n);

[tnew,sortI] = sort(tnew);
idNew = idNew(sortI);
tMain = tMain(sortI);

%
lNew = length(idNew);
E = populateSeisCompStructure(lNew);
for i = 1:lNew
    E(i,1) = readSCBulletin(idNew(i));
    disp(i);
end

%
%HB19.HHN is damaged
%HB17.HHZ and HB17.HHN are damaged
kstnms = ["HB15";"HB18";"HB12";"HB17";"HB19";"PTGL";...
    "HB15";"HB18";"HB12";"HB17";"HB19";"PTGL";...
    "HB15";"HB18";"HB12";"HB17";"HB19";"PTGL"];

knetwk_ = [repmat("XF",5,1);"EC"];
knetwks = repmat(knetwk_,3,1);
kcmpnms = [repmat("HHZ",6,1); repmat("HHN",6,1); repmat("HHE",6,1)];
kholes = repmat("",18,1);

%%
newFs = 40;
lfc = 2;
hfc = 12;
lKstnms = length(kstnms);
allTemplates = [];
for i = 1:lKstnms
    kstnm_ = kstnms(i);
    kcmpnm_ = kcmpnms(i);
    khole_ = kholes(i);
    knetwk_ = knetwks(i);

    D = extractWaveforms(tMain(1:lNew)-seconds(0),seconds(30),...
        kstnm_,kcmpnm_,knetwk_,khole_,true,true,[],1,true,[lfc,hfc,false,false,false,newFs]);
    D = resampleWaveforms(D,newFs);
    allTemplates = [allTemplates; ...
        pull(D)];
end

%
clear D

%%
T = populateWaveforms([lKstnms lNew]);
for i = 1:lNew
    disp(i);
    at_ = allTemplates(:,i);
    at_ = reshape(at_,[1201 length(at_)/1201]);
    T_ = populateWaveforms(lKstnms);
    for j = 1:size(at_,2)
        T_(j) = dealHeader(T_(j),at_(:,j),newFs,tMain(i));
        T_(j).kstnm = kstnms(j);
        T_(j).kcmpnm = kcmpnms(j);
        T_(j).knetwk = knetwks(j);
        T_(j).khole = kholes(j);
        T_(j).eqmag = E(i).mag;
        T_(j).evla = E(i).lat;
        T_(j).evlo = E(i).lon;
        T_(j).evdp = E(i).depth;
        T_(j).evid = E(i).id;
    end
    T(:,i) = T_;
end

%%
% close all;
% plotWaveforms(T(:,4));

%%
% lfc = 0.75;
% hfc = 6;
% diffFlag = false;
% plotFlag = false;
% writeFlag = true;
% minNStations = 9;
% diasFlag = false;
% maxDist = 100;
% 
% for i = 13:length(idNew)
%     close all;
%     %tnew,newMags,pks,matches,E,pTime,ccnorm,T,stack,templates,S,t,locs_
%     [tnew,newMags,pks,matches,E,pTime,ccnorm,T,stack,templates,S,t,locs_] = ...
%         findDuplicates(idNew(i),(datetime(2022,03,20):datetime(2022,04,20))',...
%         100,100,2,18,0.2,lfc,hfc,diffFlag,plotFlag,writeFlag,minNStations,diasFlag,maxDist);
% end
