%function mergeQ330(stnm,searchComponent)
%if nargin < 1
%this script is superceded by "dataselect" which is a passoft command line
%tool
clear; clc;
stnm = "FER1";

searchComponent = "BHZ";
mseedDirName = sprintf("~/data/galapagos/merged/%s.D/",searchComponent);
if ~exist(mseedDirName,"dir")
    mkdir(mseedDirName);
end
%end

%cd ~/data/VCH1/June2023/DISK1/data
%cd ~/data/FER1/data/;
%cd(['~/data/galapagos/',char(string(stnm)),'_Nov2023/data/']);
cd(['~//data/galapagos/RawIGEPNQ330Data/',char(string(stnm))]);

files = dir('EC-*');
lFiles = length(files);
%fmtStr = "NN.SSSSS.LC.CCC.T";

quality = "D";

i = 1;
fName = files(i).name;
fprintf("%s\n",fName);
S = readMiniSeed(fName);
kcmpnms = pull(S,'kcmpnm');
[lia,locb] = ismember(searchComponent,kcmpnms);
if sum(lia)
    Z(1) = S(locb);
    knetwk = Z.knetwk;
    kstnm = Z.kstnm;
    khole = Z.khole;
    kcmpnm = Z.kcmpnm;
    mseedFileName = sprintf("%s.%s.%s.%s",knetwk,kstnm,khole,kcmpnm);
end

newFs = 1./Z.delta;
tStart = Z.ref; %dateshift(Z.ref,'start','day');
tEnd = dateshift(tStart,'end','day');
day1 = dateshift(Z.ref,'start','day');
day2 = dateshift(Z.ref+Z.e,'start','day');
if day2 > day1
    mseedFullFileName = fullfile(mseedDirName,mseedFileName);
    [Z,tStart,tEnd,Zorig] = writeZ(Z,tStart,mseedFullFileName);
end

%%
for i = 2:lFiles
    fName = files(i).name;
    %fprintf("%s\n",fName);
    S = readMiniSeed(fName);
    kcmpnms = pull(S,'kcmpnm');
    [lia,locb] = ismember(searchComponent,kcmpnms);
    if sum(lia)
        S = S(locb);
    end

    sStart = S.ref;
    if sStart >= tEnd
        tStart = sStart;
        tEnd = dateshift(tStart,'end','day');
        Z = S;
    else
        Z = mergeWaveforms([Z;S]);
    end

    zEndTime = Z.ref + Z.e;
    if zEndTime < tEnd% && i < lFiles
        %fprintf('do not write yet, keep adding stations\n');
        continue;
    end

    %%
    mseedFullFileName = fullfile(mseedDirName,mseedFileName);
    [Z,tStart,tEnd,Zorig] = writeZ(Z,tStart,mseedFullFileName);
end

function [Z,tStart,tEnd,Zorig] = writeZ(Z,tStart,mseedFullFileName)
Zorig = Z;
Z = cutWaveforms(Zorig,tStart,0,dateshift(tStart,'end','day')-tStart); %,true,true);

%
d = single(Z.d);
t0 = datenum(Z.ref);
newFs = 1./Z.delta;
mkmseed(char(mseedFullFileName),d,t0,newFs,11);

zEndTime = Z.ref + Z.e;
tStartOrig = tStart;
tStart = zEndTime + seconds(1/newFs);
tEnd = dateshift(tStart,'end','day');
newEnd = Zorig.ref+Zorig.e;
if newEnd > tStart
    Z = cutWaveforms(Zorig,tStart,0,newEnd-tStart); %,true,true);
end
end
