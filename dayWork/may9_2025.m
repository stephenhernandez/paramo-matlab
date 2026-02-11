%may9_2025
clear; 
custom_prefix = "ggp_no_pino";
dayStart = datetime(2020,01,01);
[evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
    timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
    locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1(dayStart);
minLat = -0.25;
maxLat = -0.09;
minLon = -78.70;
maxLon = -78.5;
minMag = -0.5;

%%
variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";...
    "templateNumber";"madMain";"nUsedMain";"polMain";"evidMain"};
[tMain,ccMain,ampMain,dMag,magMain,templateNumber,madMain,nUsedMain,polMain,evidMain] = ...
    loadRepeaterCatalog(custom_prefix,variableNames,dayStart);
[tMain,madMain,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain] = ...
    filterCatalog(tMain,madMain,60,ampMain,dMag,magMain,templateNumber,ccMain,nUsedMain,polMain,evidMain);

%%
ccI = (ccMain >= 0.17 & madMain >= 17) & magMain >= -1 & nUsedMain >= 3 & tMain >= dayStart;

close all;

% Figure 1
axL = linkedPlot(tMain(ccI),[(0:sum(ccI)-1)' ccMain(ccI) ampMain(ccI) madMain(ccI) ...
    magMain(ccI) templateNumber(ccI) nUsedMain(ccI)]); zoom on;
axL(3).YScale = 'log';
axL(4).YScale = 'log';

% Figure 2
figure();
plot(tMain(ccI),nUsedMain(ccI),'.','linewidth',3); zoom on; grid on;
ylabel('n_used')

% Figure 3
figure();
histogram(tMain(ccI),...
    (dateshift(min(tMain(ccI)),'start','hour'):hours(1):dateshift(max(tMain(ccI)),'end','hour'))');
zoom on; title('numero de detecciones cada hora'); grid on;

%%
ggpI = eqlat >= minLat & eqlat <= maxLat & eqlon >= minLon & eqlon <= maxLon & scbullmag >= minMag;
outerBoxAdder = 0.3; 
ggpI2 = eqlat >= minLat-outerBoxAdder & eqlat <= maxLat+outerBoxAdder & eqlon >= minLon-outerBoxAdder & eqlon <= maxLon+outerBoxAdder & scbullmag >= -1;

%%
tSH = tMain(ccI);
magSH = magMain(ccI);

% SC inner square
tSC5_1 = t(ggpI);
magSC5_1 = scbullmag(ggpI);
eqlatSC5_1 = eqlat(ggpI);
eqlonSC5_1 = eqlon(ggpI);
eqdepthSC5_1 = eqdepth(ggpI);

% SC inner+outer square
tSC5_2 = t(ggpI2);
magSC5_2 = scbullmag(ggpI2);
eqlatSC5_2 = eqlat(ggpI2);
eqlonSC5_2 = eqlon(ggpI2);
eqdepthSC5_2 = eqdepth(ggpI2);

% SC just the outer "ring"
[~,~,~,just_t2I] = ...
    synchronizeCatalog(tSC5_1,tSC5_2,2,true);
tSC5_outside = tSC5_2(just_t2I);
magSC5_outside = magSC5_2(just_t2I);
eqlatSC5_outside = eqlatSC5_2(just_t2I);
eqlonSC5_outside = eqlonSC5_2(just_t2I);
eqdepthSC5_outside = eqdepthSC5_2(just_t2I);

%%
[t1_commonI,t2_commonI] = ...
    synchronizeCatalog(tSH,tSC5_outside,20,false);
delI = false(size(tSH));
delI(t1_commonI) = true;
tSHOrig = tSH;
magSHOrig = magSH;
tSH(delI) = [];
magSH(delI) = [];

% remove events that are likely in the outer ring...
tSC5_removed = tSC5_outside(t2_commonI);
magSC5_removed = magSC5_outside(t2_commonI);
eqlatSC5_removed = eqlatSC5_outside(t2_commonI);
eqlonSC5_removed = eqlonSC5_outside(t2_commonI);
eqdepthSC5_removed = eqdepthSC5_outside(t2_commonI);

%%
[~,t2_commonI,just_t1I,just_t2I] = ...
    synchronizeCatalog(tSH,tSC5_1,20,false);

tSC5_unified = [tSC5_1(t2_commonI); tSC5_1(just_t2I)];
magSC5_unified = [magSC5_1(t2_commonI); magSC5_1(just_t2I)];
[tSC5_unified,sortI] = sort(tSC5_unified);
magSC5_unified = magSC5_unified(sortI);

% merge SC% and stephens recent custom cat...
tFinal = [tSC5_1(t2_commonI); tSH(just_t1I); tSC5_1(just_t2I)];
magFinal = [magSC5_1(t2_commonI); magSH(just_t1I); magSC5_1(just_t2I)];
[tFinal,sortI] = sort(tFinal);
magFinal = magFinal(sortI);

%%
variableNames = {"tMain";"ccMain";"ampMain";"dMag";"magMain";"templateNumber";"madMain";"nUsedMain"};
[tGGP,ccGGP,ampGGP,dMagGGP,magGGP,templateNumberGGP,madGGP,nUsedGGP] = ...
    loadRepeaterCatalog("ggp",variableNames,dayStart);
[tGGP,ccGGP,ampGGP,dMagGGP,magGGP,templateNumberGGP,madGGP,nUsedGGP] = ...
    filterCatalog(tGGP,ccGGP,10,ampGGP,dMagGGP,magGGP,templateNumberGGP,madGGP,nUsedGGP);

ampThresh = 1e3;
%ggpI = templateNumberGGP == 1 & madGGP >= 11.5 & ccGGP >= 0.165 & nUsedGGP > 3 & tGGP >= dayStart & ampGGP >= ampThresh;
ggpI1 = (templateNumberGGP == 1 & madGGP >= 11.5 & ccGGP >= 0.165 & nUsedGGP > 3 & tGGP >= dayStart & ampGGP >= ampThresh);
dMagGGP(ggpI1) = dMagGGP(ggpI1) - 3.8;
ggpI2 = (templateNumberGGP == 2 & madGGP >= 18 & ccGGP >= 0.18 & nUsedGGP > 3 & tGGP >= dayStart & ampGGP >= ampThresh);
dMagGGP(ggpI2) = dMagGGP(ggpI2) - 6.35;
ggpI = ggpI1 | ggpI2;

[tGGP,ccGGP,ampGGP,dMagGGP,magGGP,templateNumberGGP,madGGP,nUsedGGP] = ...
    deal(tGGP(ggpI),ccGGP(ggpI),ampGGP(ggpI),dMagGGP(ggpI),magGGP(ggpI),templateNumberGGP(ggpI),madGGP(ggpI),nUsedGGP(ggpI));

%%
[t1_commonI,t2_commonI,~,just_t2I] = synchronizeCatalog(tFinal,tGGP,5,true);
tFinal(t1_commonI) = tGGP(t2_commonI);
magFinal(t1_commonI) = dMagGGP(t2_commonI);
tFinal = [tFinal; tGGP(just_t2I)];
magFinal = [magFinal; dMagGGP(just_t2I)];
[tFinal,sortI] = sort(tFinal);
magFinal = magFinal(sortI);
fI = magFinal >= minMag;
tFinal = tFinal(fI);
magFinal = magFinal(fI);

%%
% Figure 4
figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(tFinal,...
    (dateshift(min(tFinal),'start','day'):days(1):dateshift(max(tFinal),'end','day'))');
ax(1).XTick = [];
yyaxis right;
plot(tFinal,(0:length(tFinal)-1)','.');
ax(1).YAxis(1).Label.String = "Numero Diario";
ax(1).YAxis(2).Label.String = "Numero Acumulativo";

ax(2,1) = nexttile();
plot(tFinal,magFinal,'o','linewidth',1); zoom on;
ax(2).YLabel.String = "magnitud";
linkaxes(ax,"x");
grid(ax,"on");
xlim tight;
hold(ax(2),"on");

%%
% figure(4);
% plot(ax(2),tSC5_removed,magSC5_removed,'.'); zoom on; grid on;
figure(4);
plot(ax(2),tSC5_unified,magSC5_unified,'.'); zoom on; grid on;

% %%
% % figure(4);
% % plot(ax(2),tGGP(t2_commonI),dMagGGP(t2_commonI),'.'); zoom on; grid on;
% figure(4);
% plot(ax(2),tGGP(just_t2I),dMagGGP(just_t2I),'.'); zoom on; grid on;

%%
nDays = 30;
[rate,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin] = t2r(tFinal,days(nDays),magFinal);
axL = linkedPlot(tFinal,[rate/nDays,meanMagsFixedTimeWin,medianMagsFixedTimeWin,sumEnergyFixedTimeWin/nDays],'.');
axL(4).YScale = "log";

%%
[N,edges] = histcounts(tFinal,(dateshift(min(tFinal),'start','day'):days(1):dateshift(max(tFinal),'end','day'))');
figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(ax,'BinEdges',edges','BinCounts',N); zoom on; grid on;
ax(2,1) = nexttile();
plot(ax(2),tFinal,rate/nDays,'.'); zoom on; grid on;
linkaxes(ax,'x');

%% figure 7
[N,edges] = histcounts(tFinal,(dateshift(min(tFinal),'start','day'):days(1):dateshift(max(tFinal),'end','day'))');
figure();
tiledlayout(2,1,'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile();
histogram(ax,'BinEdges',edges','BinCounts',N); zoom on; grid on; 
ylabel("Numero Diario"); 
xlim([datetime(2013,01,01) datetime(2025,05,30)]); 
yyaxis right;
plot(tFinal(magFinal>= 3),magFinal(magFinal>=3)','p');  
ax(1).YAxis(2).Label.String = "Magnitud";

ax(2,1) = nexttile();
histogram(ax(2),'BinEdges',edges','BinCounts',N); zoom on; grid on; 
xlim([datetime(2022,01,01) dateshift(datetime("now"),"end","day")]); 
ylabel("Numero Diario"); 
yyaxis right;
plot(tFinal(magFinal>= 3),magFinal(magFinal>=3)','p');  ax(2).YAxis(2).Label.String = "Magnitud";
ax(2).YLim = [3 4.5];

length(magFinal)
length(magFinal)/(days(datetime(2025,05,30)-datetime(2013,01,01))+1)