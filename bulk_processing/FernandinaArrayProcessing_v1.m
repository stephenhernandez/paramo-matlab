clear;
close all;
clc;
cd ~/data/Galapagos/FERNANDINA_2023/SERVICE2/
nUsed = [];
tMain = [];
ampMain = [];
bazMain = [];
velMain = [];
meanCCsMain = [];
medCCsMain = [];
dirFiles = dir('202*');
for j = 1:length(dirFiles)
    dayName = fullfile('~/data/Galapagos/FERNANDINA_2023/SERVICE2/',dirFiles(j).name);
    cd(dayName)
    files = dir("*_FAR*_HHZ.m");
    lfiles = length(files);
    if lfiles < 5
        continue;
    end
    
    clear S;
    for i = 1:lfiles
        S(i,1) = readMiniSeed(files(i).name);
    end
    disp(dayName);
    [~,nUsed_,tMain_,ampMain_,bazMain_,velMain_,meanCCsMain_,medCCsMain_] = process_fernandina_array_v1(S);
    nUsed = [nUsed; nUsed_];
    tMain = [tMain; tMain_];
    ampMain = [ampMain; ampMain_];
    bazMain = [bazMain; bazMain_];
    velMain = [velMain; velMain_];
    meanCCsMain = [meanCCsMain; meanCCsMain_];
    medCCsMain = [medCCsMain; medCCsMain_];
end