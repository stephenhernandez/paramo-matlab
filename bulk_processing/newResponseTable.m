% generate new table of response info
% run updateAll.sh in the response directory
clear; close all; clc;

%%
cd ~/masa/response/;

files = dir('*/SAC_PZs_*');
ll = length(files);
ZP = populatePoleZeroStruct(ll);
for i = 1:ll
    disp(i);
    fname = fullfile(files(i).folder,files(i).name);

    [zz,pp,constant,knetwk,kstnm,khole,kcmpnm,tStart,tEnd,stla,stlo,stel,stde,...
        stdi,staz,Fs,instType,instGain,sensitivity,a0] = read_sac_pole_zero(fname);
    ZP(i).zeros = zz;
    ZP(i).poles = pp;
    ZP(i).constant = constant;
    ZP(i).knetwk = knetwk;
    ZP(i).kstnm = kstnm;
    ZP(i).khole = khole;
    ZP(i).kcmpnm = kcmpnm;
    ZP(i).tStart = tStart;
    ZP(i).tEnd = tEnd;
    ZP(i).stla = stla;
    ZP(i).stlo = stlo;
    ZP(i).stel = stel;
    ZP(i).stde = stde;
    ZP(i).stdi = stdi;
    ZP(i).staz = staz;
    ZP(i).Fs = Fs;
    ZP(i).instType = instType;
    ZP(i).instGain = instGain;
    ZP(i).sensitivity = sensitivity;
    ZP(i).a0 = a0;
end

%%
clearvars -except ZP
kstnm2 = pull(ZP,'kstnm');
knetwk2 = pull(ZP,'knetwk');
kcmpnm2 = pull(ZP,'kcmpnm');
khole2 = pull(ZP,'khole');
stla2 = pull(ZP,'stla');
stlo2 = pull(ZP,'stlo');
stde2 = pull(ZP,'stde');
Fs2 = pull(ZP,'Fs');
sensitivity2 = pull(ZP,'sensitivity');
stel2 = pull(ZP,'stel');
a02 = pull(ZP,'a0');
constant2 = pull(ZP,'constant');
tStart2 = pull(ZP,'tStart');
tEnd2 = pull(ZP,'tEnd');
allSNCLs2 = strcat(knetwk2,kstnm2,khole2,kcmpnm2);

%%
ZPT2 = struct2table(ZP);
clear ZP
load ~/igdata/ecuadorSensorDataTable10.mat

%
lia = ismember(allSNCLs,allSNCLs2);
kstnm(lia) = [];
knetwk(lia) = [];
kcmpnm(lia) = [];
khole(lia) = [];
stla(lia) = [];
stlo(lia) = [];
stde(lia) = [];
Fs(lia) = [];
sensitivity(lia) = [];
stel(lia) = [];
a0(lia) = [];
constant(lia) = [];
tStart(lia) = [];
tEnd(lia) = [];
allSNCLs(lia) = [];
ZPT(lia,:) = [];

kstnm = [kstnm; kstnm2];
knetwk = [knetwk; knetwk2];
kcmpnm = [kcmpnm; kcmpnm2];
khole = [khole; khole2];
stla = [stla; stla2];
stlo = [stlo; stlo2];
stde = [stde; stde2];
Fs = [Fs; Fs2];
sensitivity = [sensitivity; sensitivity2];
stel = [stel; stel2];
a0 = [a0; a02];
constant = [constant; constant2];
tStart = [tStart; tStart2];

tEnd = [tEnd; tEnd2];
allSNCLs = [allSNCLs; allSNCLs2];
ZPT = [ZPT; ZPT2];

%%
clearvars *2
[~,ia] = unique([allSNCLs datestr(tStart) datestr(tEnd)],'rows');
kstnm = kstnm(ia);
knetwk = knetwk(ia);
kcmpnm = kcmpnm(ia);
khole = khole(ia);
stla = stla(ia);
stlo = stlo(ia);
stde = stde(ia);
Fs = Fs(ia);
sensitivity = sensitivity(ia);
stel = stel(ia);
a0 = a0(ia);
constant = constant(ia);
tStart = tStart(ia);
tEnd = tEnd(ia);
allSNCLs = allSNCLs(ia);
ZPT = ZPT(ia,:);

%% experimental
nExperiment = 0;
while nExperiment < 5
    uniqSNCLs = unique(allSNCLs);
    lAll = length(uniqSNCLs);
    removeI = [];
    for i = 1:lAll
        disp(i);
        thisSNCL = uniqSNCLs(i);
        lia = find(ismember(allSNCLs,thisSNCL));
        tStart_ = ZPT.tStart(lia);
        tEnd_ = ZPT.tEnd(lia);
        [tStart_,sortlia] = sort(tStart_);
        lia = lia(sortlia);
        tEnd_ = tEnd_(sortlia);
        llia = length(lia);
        if (llia < 1)
            disp('continue');
        end
        tDiff = seconds(diff(tStart_));
        repeatI = tDiff == 0;
        if sum(repeatI)
            ri1 = find(repeatI);
            ri1 = ri1(1);
            ri2 = ri1+1;
            tEnd1 = tEnd_(ri1);
            tEnd2 = tEnd_(ri2);
            if tEnd1 < tEnd2
                deletei_ = lia(ri2);
            else
                deletei_ = lia(ri1);
            end
            removeI = [removeI; deletei_];
        end
    end

    kstnm(removeI) = [];
    knetwk(removeI) = [];
    kcmpnm(removeI) = [];
    khole(removeI) = [];
    stla(removeI) = [];
    stlo(removeI) = [];
    stde(removeI) = [];
    Fs(removeI) = [];
    sensitivity(removeI) = [];
    stel(removeI) = [];
    a0(removeI) = [];
    constant(removeI) = [];
    tStart(removeI) = [];
    tEnd(removeI) = [];
    allSNCLs(removeI) = [];
    ZPT(removeI,:) = [];
    clear deletei_ i ia lAll lia llia removeI repeatI ri1 ri2 sortlia tDiff tEnd1 tEnd2 tEnd_ tStart_ thisSNCL uniqSNCLs
    nExperiment = nExperiment + 1;
end
clear nExperiment;

%%
instType = ZPT.instType;
sensitivity = ZPT.sensitivity;
accI = (contains(instType,"CMG5T") | contains(instType,"CMG-5T")) & a0 < 1;
table(unique(instType(accI)),groupcounts(instType(accI)))

t5lI = instType == "CMG5T,serial number: T5L...,guralp CMG5T" & a0 > 1;
goodT5 = ZPT(t5lI,:);
goodT5 = goodT5(1,:);
t5lI = ((instType == "CMG5T,serial number: T5L...,guralp CMG5T" | instType == "CMG5T,serial number: T5R...,guralp CMG5T") ...
    & a0 < 1) | (instType == "CMG5T,serial" & sensitivity < 1) | ...
    (instType == "CMG5T,serial number: T5R...,guralp CMG5T" & sensitivity < 1);
badT5 = ZPT(t5lI,:);
for i = 1:sum(t5lI)
    bad_ = badT5(i,:);
    bad_.zeros = goodT5.zeros;
    bad_.poles = goodT5.poles;
    bad_.a0 = goodT5.a0;
    instGainTmp = bad_.instGain;
    bad_.sensitivity = 312500*instGainTmp;
    bad_.constant = bad_.sensitivity.*bad_.a0;
    badT5(i,:) = bad_;
end
ZPT(t5lI,:) = badT5;

sensitivity = ZPT.sensitivity;
a0 = ZPT.a0;
constant = ZPT.constant;

clear accI goodT5 bad_ badT5 instGainTmp i t5lI ans

%% 5TC
instType = ZPT.instType;
accI = contains(instType,"5TC") & a0 > 1;
goodT5 = ZPT(accI,:);
goodT5 = goodT5(1,:);
t5lI = contains(instType,"5TC") & a0 < 1;
badT5 = ZPT(t5lI,:);
%1.02            427819    3022009000000
for i = 1:sum(t5lI)
    bad_ = badT5(i,:);
    bad_.zeros = goodT5.zeros;
    bad_.poles = goodT5.poles;
    bad_.a0 = 3022009000000; %goodT5.a0;
    bad_.sensitivity = 427819;
    bad_.constant = bad_.sensitivity.*bad_.a0;
    badT5(i,:) = bad_;
end

ZPT(t5lI,:) = badT5;
sensitivity = ZPT.sensitivity;
a0 = ZPT.a0;
constant = ZPT.constant;

clear accI goodT5 bad_ badT5 instGainTmp i t5lI ans