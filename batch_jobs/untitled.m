clear; close all;
cd ~/data/FER1/data/;
files = dir('EC-*');
lFiles = length(files);
Z = populateWaveforms(lFiles); n = 0;
for i = 1:3%lFiles
    fName = files(i).name;
    fprintf("%s\n",fName);
    S = readMiniSeed(fName);
    kcmpnms = pull(S,'kcmpnm');
    [lia,locb] = ismember("BHZ",kcmpnms);
    if sum(lia)
        n = n+1;
        Z(n,1) = S(locb);
    end
end