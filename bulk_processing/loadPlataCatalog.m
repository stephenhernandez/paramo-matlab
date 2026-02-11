clear; close all; clc;
cd ~/templateSearch/

dirs = dir('Family*');
dataMain = [];
for i = 1:length(dirs)
    dirName = dirs(i).name;
    cd(dirName);
    fprintf("processing: %s\n",dirName);
    files = dir('dayTemplateSearch*.txt');
    lFiles = length(files);

    maxN = 1e5;
    data_ = NaN(maxN,13);
    n = 0;
    for j = 1:lFiles
        data = importdata(files(j).name);
        nrows = size(data,1);
        if nrows
            data_(n+1:n+nrows,:) = data;
            n = n+nrows;
        end
    end
    data_ = data_(1:n,:);
    data_(:,10) = i;

    dataMain = [dataMain; data_];
    cd ..
end

%%
tMain = datetime(dataMain(:,1:6));
[tMain,sI] = sort(tMain);
dataMain = dataMain(sI,:);
dataMainOrig = dataMain;
clearvars -except dataMainOrig tMain
cd ~/research/now/plata;
%save("plata_catalog_10JAN2023");

%%
close all; 
tMain = datetime(dataMainOrig(:,1:6));
[tMain,sI] = sort(tMain);
dataMain = dataMainOrig(sI,:);
[t,cc,removeIndices,keepIndices] = removeRepeatedMatches(tMain,dataMainOrig(:,11),90);
for i = 1:length(removeIndices)
rI = removeIndices{i};
dataMain(rI,:) = [];
end
cc1 = dataMain(:,9);
cc2 = dataMain(:,11);
nUsed = dataMain(:,13);
dataMain(1,:)
ccI = cc1 >= 0.06 & cc1 <= 1 & cc2 >= 10 & nUsed >= 9;% & t >= datetime(2020,01,01);
sum(ccI)

close all; figure(); ax(1) = subplot(311); semilogy(t(ccI),cc1(ccI),'.'); zoom on; grid on; hold on; semilogy(t(~ccI),cc1(~ccI),'.'); zoom on; grid on;
ax(2) = subplot(312); semilogy(t(ccI),cc2(ccI),'.'); zoom on; grid on; hold on; semilogy(t(~ccI),cc2(~ccI),'.');
ax(3) = subplot(313); semilogy(t(ccI),dataMain(ccI,9),'.'); zoom on; grid on; hold on; semilogy(tMain(~ccI),dataMain(~ccI,9),'.'); linkaxes(ax,'x');
figure(); plot(t(ccI),1:sum(ccI),'.'); zoom on; grid on;
figure(); semilogy(t(ccI),dataMain(ccI,end),'.'); zoom on; grid on;
