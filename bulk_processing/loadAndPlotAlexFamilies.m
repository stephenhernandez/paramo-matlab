close all;
dataMain = []; 
load ~/igdata/ec_boundaries.mat;
cd ~/research/now/plata/laplataregiontemplatesfamiliesfromgrowclust20112016/;
files = dir('group_*.txt');
lFiles = length(files);
for i = 1:lFiles
    tmpFile = files(i).name;
    if length(tmpFile) == 18
        famNumber = str2double(tmpFile(13:14));
    elseif length(tmpFile) == 17
        famNumber = str2double(tmpFile(13));
    end
    groupNumber = str2double(tmpFile(8));
    data = importdata(tmpFile);
    data = [data repmat(groupNumber,size(data,1),1) repmat(famNumber,size(data,1),1)];
    dataMain = [dataMain; data];
end

tMain = datetime(dataMain(:,1:6));
[tMain,sI] = sort(tMain);
dataMain = dataMain(sI,:);
lon = dataMain(:,7);
lat = dataMain(:,8);
depth = dataMain(:,9);
groupNumber = dataMain(:,10);
familyNumber = dataMain(:,11);

badI = seconds(diff(tMain)) == 0;
dataMain(badI,:) = [];
tMain = datetime(dataMain(:,1:6));
lon = dataMain(:,7);
lat = dataMain(:,8);
depth = dataMain(:,9);
groupNumber = dataMain(:,10);
familyNumber = dataMain(:,11);

%%
for i = 1:6
    gI = groupNumber == i;
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(8,1,[1 2]); 
    plot(tMain(gI),1:sum(gI),'.'); zoom on; title(strcat("Group: ",num2str(i)));
    grid on; ax = gca; ax.LineWidth = 1.5;
    ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';

    grid on
    uniqFams = unique(familyNumber(gI));
    nFams = 0;
    legStr = [];
    for j = 1:length(uniqFams)
        nFams = nFams + 1;
        fI = gI & familyNumber == uniqFams(j);
        if nFams <= 7
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'p','linewidth',1); zoom on;
        elseif nFams <= 14
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'h','linewidth',1); zoom on;
        elseif nFams <= 21
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'o','linewidth',1); zoom on;
        elseif nFams <= 28
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'d','linewidth',1); zoom on;
        elseif nFams <= 35
            subplot(8,1,[3 4 5 6 7 8]); plot(lon(fI),lat(fI),'*','linewidth',1); zoom on;
        end

        %title(strcat("Group: ",num2str(i)));
        hold on; grid on; axis equal;

        legStr = [legStr; ...
            strcat("Family ",num2str(uniqFams(j)),"(",num2str(sum(fI)),")")];
        axis([min(lon)-0.1 max(lon)+0.1 min(lat)-0.1 max(lat)+0.1]);
    end
    hold on; plot(lonEC,latEC,'k-','linewidth',2);
    legend(legStr,'Location','Best');
    grid on; ax = gca; ax.LineWidth = 1.5;
    ax.GridAlpha = 0.2; zoom on; ax.Box = 'on';
end
