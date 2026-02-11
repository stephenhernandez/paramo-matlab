clear; %vars -except S;
close all;
tic;
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    loadRepeaterCatalog("vlps_cotopaxi");

[tMain,sortI] = sort(tMain);
ccMain = ccMain(sortI);
ampMain = ampMain(sortI);
dMag = dMag(sortI);
magMain = magMain(sortI);
templateNumber = templateNumber(sortI);
madMain = madMain(sortI);
nUsedMain = nUsedMain(sortI);

[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain] = ...
    filterCatalog(tMain,ccMain,10,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain);

close all;
sortI = (ccMain >= 0.20 | madMain >= 12) & ampMain >= 1e1 & nUsedMain >= 6 & ...
    tMain >= datetime(2000,01,01);

% sortI = (ccMain >= 0.10 | madMain >= 10) & ampMain >= 1e1 & nUsedMain >= 6 & ...
%     tMain >= datetime(2000,01,01);

tMain = tMain(sortI);
ccMain = ccMain(sortI);
ampMain = ampMain(sortI);
dMag = dMag(sortI);
magMain = magMain(sortI);
templateNumber = templateNumber(sortI);
madMain = madMain(sortI);
nUsedMain = nUsedMain(sortI);
clear sortI;

ax1 = linkedPlot(tMain,[ampMain ccMain madMain nUsedMain],".","compact");
linkaxes(ax1,'x');
toc;

%%
close all;
cd("~/masa/template_search/vlps_cotopaxi/");
toc;

transferFlag = true;
diffFlag = true;
waveI = true(size(tMain));
waveformKstnms = ["BREF";"CO1V";"BVC2";"BTAM"];
if exist("CotopaxiVLPCatalog_v3.mat","file")
    if ~exist("S","var")
        fprintf("loading waveform file\n");
        load("CotopaxiVLPWaveforms","S");
    end
    refs = pull(S,"ref");
    refs1 = max(refs,[],2);
    badI = ~isfinite(refs1);
    refs1(badI) = [];
    S(badI,:) = [];

    [t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(refs1,tMain,10,true,false);
    if ~isempty(just_t1I)
        S(just_t1I,:) = [];
        refs = pull(S,"ref");
        refs1 = max(refs,[],2);
        badI = ~isfinite(refs1);
        refs1(badI) = [];
        S(badI,:) = [];
        clear badI;
        [t1_commonI,t2_commonI,just_t1I,just_t2I] = synchronizeCatalog(refs1,tMain,10,true,false);
    end

    if length(t1_commonI) ~= length(t2_commonI)
        fprintf("length(t1_commonI) ~= length(t2_commonI)\n");
        fprintf("what should I do??\n");
        %do stuff
    end

    if ~isempty(just_t2I)
        tDownload = tMain(just_t2I);
        S_ = extractWaveforms(tDownload,seconds(60),waveformKstnms,["BHZ";"HHZ"],"EC","",true,true,...
            1,true,[1/20,-inf,false,transferFlag,diffFlag,20]);
        S = [S; S_(:,1:3)];

        %
        refs = pull(S,"ref");
        refs1 = max(refs,[],2);
        badI = ~isfinite(refs1);
        refs1(badI) = [];
        S(badI,:) = [];
        %
        [~,sortI] = sort(refs1);
        S = S(sortI,:);

        %%
        fprintf("saving waveform file...\n");
        save('CotopaxiVLPWaveforms',"S",'-v7.3');
        fprintf("done saving waveform file\n");
        toc;
    else
        fprintf("they are both synced, no need to download new events...\n");
    end
else
    if sum(waveI)
        S = extractWaveforms(tMain(waveI),seconds(60),waveformKstnms,"BHZ","EC","",true,true,...
            [],1,true,[1/20,-inf,false,transferFlag,diffFlag,20]);
        fprintf("saving waveform file...\n");
        save('CotopaxiVLPWaveformsOrig_03JUN2025',"S",'-v7.3');
        %save('CotopaxiVLPWaveformsOrig_17OCT2023',"S",'-v7.3');
        fprintf("done saving waveform file...\n");
    end
end

%%
Sorig = S;
stnms = pull(S,"kstnm");
refs = pull(S,"ref");
refs1 = max(refs,[],2);
[SR,SC] = size(S);

d = pull(S);
[dR,dC] = size(d);
if SC > 1
    d = reshape(d,[dR SR SC]);
    d = permute(d,[1 3 2]);
    d = reshape(d,[],SR,1);
end

d(~isfinite(d)) = 0;
badI = ~isfinite(refs1) | mean(abs(d))' < 1;
fprintf("number of not-existing events: %d\n",sum(badI));

%%
ccMain(badI) = [];
tMain(badI) = [];
ampMain(badI) = [];
dMag(badI) = [];
magMain(badI) = [];
templateNumber(badI) = [];
madMain(badI) = [];
nUsedMain(badI) = [];
S(badI,:) = [];

ccMainOrig = ccMain;
tMainOrig = tMain;
ampMainOrig = ampMain;
dMagOrig = dMag;
magMainOrig = magMain;
templateNumberOrig = templateNumber;
madMainOrig = madMain;
nUsedMainOrig = nUsedMain;

%%
lfc = 1/5;
hfc = 1/2;
Slow = filterWaveforms(detrendWaveforms(S),lfc,hfc,4,false,true);
Shigh = filterWaveforms(detrendWaveforms(S),1,8,4,false,true);

lowAmpsOrig = 0.5*(pull(Slow,'depmax') - pull(Slow,'depmin'));
highAmpsOrig = 0.5*(pull(Shigh,'depmax') - pull(Shigh,'depmin'));
fIndexOrig = lowAmpsOrig./highAmpsOrig;
fIndexOrig = median(fIndexOrig,2,"omitnan");
% [fIndexOrig,maxFII] = max(fIndexOrig,[],2,'linear');
% lowAmpsOrig = lowAmpsOrig(maxFII);
% highAmpsOrig = highAmpsOrig(maxFII);
lowAmpsOrig = median(lowAmpsOrig,2,"omitnan");
highAmpsOrig = median(highAmpsOrig,2,"omitnan");

%%
close all;
ratioThreshold = 1e-2; %0.1;
minLowAmp = 5e1;
ccThresh = 0.4;
madThresh = 14;
try
    vlpI = fIndexOrig >= ratioThreshold & ...
        (ccMainOrig >= ccThresh & madMainOrig >= madThresh) & ...
        lowAmpsOrig >= minLowAmp & tMainOrig >= datetime(2000,01,01);
catch
    whos
end
vlpIOrig = vlpI;

close all;
axL = linkedPlot(tMainOrig(vlpI),[(0:sum(vlpI)-1)' ccMainOrig(vlpI) ...
    ampMainOrig(vlpI) madMainOrig(vlpI) lowAmpsOrig(vlpI) ...
    templateNumberOrig(vlpI) nUsedMainOrig(vlpI)]); zoom on;
axL(3).YScale = "log";
axL(5).YScale = "log";

%Figure 1
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(tMainOrig(vlpI),fIndexOrig(vlpI),'.'); zoom on; grid on;
hold on;
semilogy(tMainOrig(~vlpI),fIndexOrig(~vlpI),'.','Linewidth',2); zoom on; grid on;
ylabel("Frequency Index [2 - 5 seconds / 1 - 8 Hz]");
pause(2);
legend("Kept","Rejected");

%Figure 3
figure('units','normalized','outerposition',[0 0 3/4 1]);
clear ax;
ax = subplot(211);
semilogy(tMainOrig(vlpI),lowAmpsOrig(vlpI),'.');
zoom on; grid on;
ax(2) = subplot(212);
ax(2).ColorOrderIndex = 2;
semilogy(tMainOrig(~vlpI),lowAmpsOrig(~vlpI),'.');
zoom on; grid on;
linkaxes(ax,'xy');
sgtitle("Amplitude, low-frequency band");
ylabel(ax(1),"Kept Events");
ylabel(ax(2),"Rejected Events");
pause(1);

%Figure 4
figure('units','normalized','outerposition',[0 0 3/4 1]);
clear ax;
ax = subplot(211);
semilogy(tMainOrig(vlpI),highAmpsOrig(vlpI),'.');
zoom on; grid on;
ax(2) = subplot(212);
ax(2).ColorOrderIndex = 2;
semilogy(tMainOrig(~vlpI),highAmpsOrig(~vlpI),'.');
zoom on; grid on;
linkaxes(ax,'xy');
sgtitle("Amplitude, high-frequency band");
ylabel(ax(1),"Kept Events");
ylabel(ax(2),"Rejected Events");
pause(1);

%Figure 5
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(tMainOrig(vlpI),ccMainOrig(vlpI),'.'); zoom on; grid on; hold on;
semilogy(tMainOrig(~vlpI),ccMainOrig(~vlpI),'.');
title("CC")
legend("Kept","Rejected");

%Figure 6
figure('units','normalized','outerposition',[0 0 3/4 1]);
semilogy(tMainOrig(vlpI),madMainOrig(vlpI),'.'); zoom on; grid on; hold on;
semilogy(tMainOrig(~vlpI),madMainOrig(~vlpI),'.');
title("MAD");
legend("Kept","Rejected");

%
saveFlag = true;
if saveFlag
    cd ~/masa/template_search/vlps_cotopaxi/
    save("CotopaxiVLPCatalog_v3.mat","ampMainOrig","ccMainOrig","dMagOrig",...
        "fIndexOrig","highAmpsOrig","lowAmpsOrig","madMainOrig","magMainOrig",...
        "nUsedMainOrig","tMainOrig","templateNumberOrig","vlpIOrig","minLowAmp",...
        "vlpI","ratioThreshold","ccThresh","madThresh","badI","-v7.3");
end

%%
clearvars -except *Orig minLowAmp vlpI ratioThreshold ccThresh madThresh badI S
nDays = 7;
vlpI = vlpIOrig;
causalFlag = false;
rate = t2r(tMainOrig(vlpI|~vlpI),days(nDays),[],causalFlag);
figure('units','normalized','outerposition',[0 0 3/4 1]);
plot(tMainOrig(vlpI|~vlpI),rate/nDays,'.'); zoom on; grid on; hold on;

otherIndex = true;
ratioThreshold2 = 2e-3; %0.1;
ratioThreshold3 = 3e-3; %0.2;
if ~otherIndex
    legend(sprintf('$FI \\ge %g$',ratioThreshold));
else
    vlpI = vlpIOrig & fIndexOrig >= ratioThreshold2;
    rate = t2r(tMainOrig(vlpI),days(nDays));
    plot(tMainOrig(vlpI),rate/nDays,'.'); zoom on; grid on;

    vlpI = vlpIOrig & fIndexOrig >= ratioThreshold3;
    rate = t2r(tMainOrig(vlpI),days(nDays));
    plot(tMainOrig(vlpI),rate/nDays,'.'); zoom on; grid on;
    legend('all',...
        sprintf('$FI \\ge %g$',ratioThreshold2),...
        sprintf('$FI \\ge %g$',ratioThreshold3));
end
title(sprintf("Daily Rate (%d-Day Smoothing)",nDays));

%
vlpI = vlpIOrig & fIndexOrig >= ratioThreshold2;
rate = t2r(tMainOrig(vlpI),days(nDays),[],causalFlag);

nDayVec = unique(1./(7:7:180)');
rateRatio = [];
tMain2 = tMainOrig(vlpI);
for i = 1:length(nDayVec)
    nDays2 = 1./nDayVec(i);
    rate2 = t2r(tMain2,days(nDays2),[],false)/nDays2;
    rate1 = t2r(tMain2,days(nDays2),[],true)/nDays2;
    rateRatio = [rateRatio rate2]; %./rate1];
end
Nfilt = 3;
w = nDayVec/sum(nDayVec);
rateRatio = convFilter(rateRatio,3,true);
stackRate1 = median(rateRatio,2,"omitnan");
figure();
semilogy(tMain2,stackRate1,'.'); zoom on; grid on; hold on;


dailyRate1 = [];
nWinVec = (2:10)';
maxNwin = max(nWinVec);
t3 = seconds(tMain2 - tMain2(1));
for i = 1:length(nWinVec)
    nWin = nWinVec(i);
    G = [1; zeros(nWin-2,1); -1];
    t4 = fftfilt(G,t3);

    %dRate = nWin.*86400./t4(2*nWin:end); %causal
    dRate = nWin.*86400./t4(nWin:end-nWin); %acausal
    dRate_ = dRate(maxNwin-nWin+1:end-maxNwin+nWin);
    dailyRate1 = [dailyRate1 dRate_];
end
w = nWinVec./sum(nWinVec);
tdum = tMain2(maxNwin:end-maxNwin);

dailyRate1 = convFilter(dailyRate1,3,true);
rateWA = sum((dailyRate1).*w',2,"omitnan"); %weighted average
rateMedian = median(dailyRate1,2,"omitnan"); %average rate ratio (snr) using multiple windows
figure();
ax(1) = subplot(211);
semilogy(tdum,rateWA,'.','linewidth',2); zoom on; grid on; hold on;
ax(2) = subplot(212);
semilogy(tdum,rateMedian,'.','linewidth',2); zoom on; grid on; %legend('weighted stack','mean');
linkaxes(ax,'x');

%
figure('units','normalized','outerposition',[0 0 1 1]);
axS(1) = subplot(211);
SS = scatter(tMainOrig(vlpI),rate/nDays,8*exp(1.5*(log10(lowAmpsOrig(vlpI))-1)),lowAmpsOrig(vlpI),'o','filled');
c = colorbar('eastoutside');
cPosition = c.Position;
c.Position = [1.03*cPosition(1) cPosition(2:4)];
ax = gca; set(ax,'ColorScale','log'); clim([minLowAmp 1e2*minLowAmp]);
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; ax.Box = 'on';
grid on; zoom on; hold on; colormap turbo; SS.MarkerFaceAlpha = 0.7;
ylabel(sprintf("Tasa Diaria (Media Movil de %d dias)",nDays));
title('VLPs de Cotopaxi');
%yyaxis right;
%plot(tMainOrig(vlpI),1:sum(vlpI),'.');

axS(2) = subplot(212);
semilogy(tMainOrig(vlpI),lowAmpsOrig(vlpI),'.');
zoom on; grid on;
ylabel('Amplitud');
%c2 = colorbar('eastoutside');
%c2.Visible = 'off';
yyaxis right;
plot(tMainOrig(vlpI),1:sum(vlpI),'.');
yyaxis(axS(2),'left')
linkaxes(axS,'x');

%
writeCatalogFlag = true;
if writeCatalogFlag
    %outFile = "~/research/now/cotopaxi/CotopaxiVLPCatalog_v2.txt";
    outFile = "~/masa/template_search/vlps_cotopaxi/CotopaxiVLPCatalog_v3.txt";
    %utFile = "~/research/now/cotopaxi/CotopaxiVLPCatalog_v3.txt";
    formatSpec = '%04d %02d %02d %02d %02d %f %f %f %f %f %f %02d';
    sum(vlpI)
    [yyyy,mm,dd,HH,MM,SS] = datevec(tMainOrig(vlpI));
    str = compose(formatSpec,yyyy,mm,dd,HH,MM,SS,fIndexOrig(vlpI),lowAmpsOrig(vlpI),highAmpsOrig(vlpI),ccMainOrig(vlpI),madMainOrig(vlpI),nUsedMainOrig(vlpI));
    str = string(str);
    fileID = fopen(outFile,'w');
    fprintf(fileID,'%s\n',str);
    fclose(fileID);
end

tDay = dateshift(tMainOrig(vlpI),'start','day');
table(unique(tDay),groupcounts(tDay))