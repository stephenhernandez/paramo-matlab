clear; close all; clc;

%%
cd ~/research/now/reventador/
[tabs,z2p,NCC,Neff,p2rms,kurt] = ...
    filterUniqueEvents('~/research/now/reventador/ReventadorSubspaceDetectorResults_v6',10); %,10);

%%
minKurt = min(kurt(:,1:2),[],2,'omitnan');
maxKurt = max(kurt(:,1:2),[],2,'omitnan');
kurtRatio = maxKurt./minKurt;

%%
z2p = z2p(:,1); %max(z2p,[],2);
kurt = kurt(:,1); %max(kurt,[],2);

%%
minAmp = 150;
winlen = 25;
nboot = 1e3;
zeroPhaseFlag = false;

%%
%goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 1) & kurt < 25;
goodI = z2p >= minAmp & ~(Neff == 2 & NCC < 0.04) & ~(Neff == 0) & kurt < 25;

%%
[N,edges] = histcounts(tabs(goodI),dateshift(min(tabs),'start','day'):dateshift(max(tabs),'end','day'));
N = N';
edges = edges(1:end-1)';

figure();
stairs(edges,N,'.-');
zoom on; grid on;

%%
figure();
semilogy(tabs(goodI),z2p(goodI),'.'); zoom on; grid on;

%%
figure();
plot(tabs(goodI),1:sum(goodI),'.'); zoom on; grid on;

tgood = tabs(goodI);
z2pGood = z2p(goodI);

%%
tStart = [datetime(2014,04,01); ...
    datetime(2016,11,10);...
    datetime(2017,09,23);...
    datetime(2018,04,02);...
    datetime(2018,07,11);...
    datetime(2018,09,19);...
    datetime(2018,11,24)];

tEnd = [datetime(2014,07,23);
    datetime(2017,09,22);...
    datetime(2018,04,01);...
    datetime(2018,07,10);...
    datetime(2018,09,18);...
    datetime(2018,11,23);...
    datetime(2019,03,27)];

%%
% 0.518167539
% 0.7063
% 0.007142857
% 0.03
% 0.104193548

%%
prcntJuv = [91.45;...
    70.10;...
    81.33;...
    85.04;...
    86.48;...
    77.53];

%%
figure();
ax_(1) = subplot(311); hold on; zoom on; grid on;
ax_(2) = subplot(312); hold on; zoom on; grid on;
ax_(3) = subplot(313); hold on; zoom on; grid on;
%ax_(4) = subplot(414); hold on; zoom on; grid on;

lT = length(tStart);
for i = 2:lT
    ts_ = tStart(i);
    te_ = tEnd(i);
    I = tgood >= ts_ & tgood < te_ & z2pGood >= 150;
    sumi = sum(I);
    
    %subplot(411);
    %plot(te_,sumi,'.','markersize',20);
    
    subplot(311);
    plot(te_,sumi/days(te_-ts_),'.','markersize',35);
    ylabel('$\#$ events / day');
    
    subplot(312);
    plot(te_,median((z2pGood(I)))/days(te_-ts_),'.','markersize',35);
    ylabel('median amp. / day');
    
    subplot(313);
    %plot(prcntJuv(i-1),sumi/days(te_-ts_),'.','markersize',20);
    plot(prcntJuv(i-1),median((z2pGood(I)))/days(te_-ts_),'.','markersize',35);
    xlabel('juvenile content');
    ylabel('median amp. / day');
    fprintf('Ndays: %f\n',days(te_-ts_));
end
zoom on; grid on;
linkaxes(ax_(1:2),'x');

%%
clear; close all; clc;
cd ~/research/now/reventador/ash_shape_analysis/
files = dir('201*.xlsx');
files = files(2:end);
lFiles = length(files);

markerSize = 16/2;
symbolList = ['h';'p';'v';'s';'o';'^';'d'];
labelNames = ["solidity";"convexity";"form factor";"axial ratio"];
legStr = [];
fn = [];
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:lFiles
    fn_ = strsplit(files(i).name,'.');
    fn_ = string(fn_{1});
    fn_ = replace(fn_,"_","/");
    
    %legStr = [legStr; string(files(i).name)];
    n = 0;
    disp(files(i).name);
    [num,txt] = xlsread(files(i).name);
    %num = (num' ./ rssq(num'))';
    [lGrains,ncols] = size(num);
    legStr = [legStr; strcat(fn_,"; $N_{g}=$",string(num2str(lGrains)))];
    for j = 1:3
        x = num(:,j);
        for k = j+1:size(num,2)
            y = num(:,k);
            y = sort(x./y);
            n = n+1;
            
            plotPos = 3*(j-1)+k-1;
            if plotPos == 9
                plotPos = [8 9];
            end
            ax_(n) = subplot(3,3,plotPos);
            zoom on; grid on;
            hold on;
            
            hold on;
            plot(y,(0:lGrains-1)'/lGrains,symbolList(i),'markersize',markerSize,'linewidth',0.2); zoom on; grid on;
            %xlim([0 1]); ylim([0 1]);
            %xlabel(labelNames(j)); ylabel(labelNames(k));
            lab1 = char(labelNames(j));
            lab2 = char(labelNames(k));
            lab1(isspace(lab1)) = [];
            lab2(isspace(lab2)) = [];
            fname(n) = strcat(string(lab1),"_",string(lab2));
            title(replace(fname(n),"_","/"))
            disp(fname);
        end
    end
end
linkaxes(ax_,'y');
subplot(3,3,[8 9])
%axis off
legend(legStr,'Location','WestOutside');
% for i = 1:n
%     figure(i);
% %     ax = gca;
% %     ax.YScale = 'log';
% %     ax.XScale = 'log';
%     legend(legStr,'Location','Best');
%     figure(i);
%     cd ~/Desktop;
%     print('-djpeg',fname(i));
% end

%
cd ~/research/now/reventador/ash_shape_analysis/;
files = dir('201*.xlsx');
files = files(2:end);
lFiles = length(files);
markerSize = 15;
symbolList = ['h';'p';'v';'s';'o';'^';'d'];
labelNames = ["solidity";"convexity";"form factor";"axial ratio"];
legStr = [];
fn = [];
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
zoom on; grid on;
for i = 1:lFiles
    fn_ = strsplit(files(i).name,'.');
    fn_ = string(fn_{1});
    fn_ = replace(fn_,"_","/");
    
    [num,txt] = xlsread(files(i).name);
    [lGrains,ncol] = size(num);
    legStr = [legStr; strcat(fn_,"; N=",string(num2str(lGrains)))];
    
    for j = 1:ncol
        subplot(2,2,j); hold on;
        x = sort(num(:,j));
        plot(x,(1:lGrains)'/lGrains,symbolList(i),'markersize',8,'linewidth',0.1);
    end
end

for i = 1:ncol
    figure(2);
    subplot(2,2,i);
    ax = gca;
    %ax.YScale = 'log';
    %ax.XScale = 'log';
    grid on;
    legend(legStr,'Location','Best');
    title(labelNames(i))
end
