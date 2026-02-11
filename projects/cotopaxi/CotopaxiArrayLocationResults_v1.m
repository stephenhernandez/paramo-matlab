clear; close all; clc;

[minLon,maxLon,minLat,maxLat] = deal(-78.48,-78.385,-0.715,-0.64);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
demData = double(demData);

%%
colorbarRightShift = 0.93;
demContrast = 1/3;
ilumAz = 90;
wMark = 1.8; %"brightness"... higher number means more washed out tones
elevAngle = 65;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',demContrast,'Azimuth',ilumAz,'noplot',...
    'Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1],'LandColor',demcmap([500 4000])); %flipud(summer(512)));
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

%----------necessary to make sure axes are synched later------
cbar1 = colorbar(ax(laxes));
cbar1.Visible = 'off';
cbar1.Position = [colorbarRightShift,0.1,0.0083333333333333,0.35];
%--------------------------------------------------------------

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

refEllipse = referenceEllipsoid('wgs84');
refShifter = (maxLat-minLat)*(0.0095/0.3);
scaleRefLat0 = minLat + refShifter;
scaleRefLon0 = minLon + refShifter;

[latout1,lonout1] = reckon(scaleRefLat0,scaleRefLon0,1*1e3,90,refEllipse);
[latout5,lonout5] = reckon(scaleRefLat0,scaleRefLon0,4*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
scaleRefLon2 = scaleRefLon0;
plot(ax(laxes),[scaleRefLon0 lonout5],[scaleRefLat0 scaleRefLat0],'k-','linewidth',3);

textShifter = -refShifter/2;
textFontSize = 24;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',3);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout1 lonout1],[scaleRefLat0 scaleRefLat2],'k-','linewidth',3);
text(textShifter+lonout1,-textShifter+scaleRefLat2,"1","FontSize",textFontSize);

plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',3);
text(textShifter+lonout5,-textShifter+scaleRefLat2,"4 km","FontSize",textFontSize);

labelFontSize = 28;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);

%% permanent stations
StationMarkerSize = 22;
contour(ax(laxes),lon,lat,demData',(2000:100:6000)','w-','linewidth',0.2);
[C,h] = contour(ax(laxes),lon,lat,demData',(2000:200:6000)','k-','linewidth',1);
clabel(C,h); %,'FontSize',15,'Color','red')

ax(laxes).XAxis.FontSize = labelFontSize;
ax(laxes).YAxis.FontSize = labelFontSize;
ax(laxes).YTickLabelRotation = 30;
ax(laxes).XTickLabelRotation = 30;
ax(laxes).LineWidth = 1;

%%
cd ~/igdata/

CTX = load('CTXAll.mat');
tCTX = CTX.tMain - minutes(0);
doaCTX = CTX.DOA;
avCTX = CTX.AV;
CTXlats = CTX.elementLats;
CTXlons = CTX.elementLons;
fOrig = CTX.fOrig;
SW = CTX.SW;
lfc2 = 0.4;
hfc2 = lfc2*4;
fI = fOrig >= lfc2 & fOrig <= hfc2;
meanLinearity = sum(SW(fI,:))'/sum(fI);

clear CTX
TAM = load('TAMAll.mat');
tTAM = TAM.tMain - minutes(0);
doaTAM = TAM.DOA;
avTAM = TAM.AV;
TAMlats = TAM.elementLats;
TAMlons = TAM.elementLons;
[lia,locb] = ismember(tTAM,tCTX);
tCommon = tTAM(lia);
diaCTX2 = doaCTX(locb(lia));
doaCTX2 = doaCTX(locb(lia));
avCTX2 = avCTX(locb(lia));
meanLinearity2 = meanLinearity(locb(lia));

plot(ax(laxes),CTXlons,CTXlats,'v','MarkerFaceColor','w','MarkerEdgeColor','k');
plot(ax(laxes),TAMlons,TAMlats,'v','MarkerFaceColor','w','MarkerEdgeColor','k');

%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

[newlats,newlons] = gcxgc(repmat(mean(CTXlats),size(doaCTX2)),...
    repmat(mean(CTXlons),size(doaCTX2)),doaCTX2,repmat(mean(TAMlats),...
    size(doaTAM)),repmat(mean(TAMlons),size(doaTAM)),doaTAM);

tremLon = newlons(:,2);
tremLat = newlats(:,2);
avDiff = abs(avTAM - avCTX2);
tremI = avCTX2 >= 950 & avTAM >= 950 & avCTX2 <= 1250 & avTAM <= 1250 & ...
    avDiff <= 200 & tremLon < 0 & tremLat < 0;

SS = bubblechart(ax(laxes),tremLon(tremI),tremLat(tremI),...
    (meanLinearity2(tremI)),avCTX2(tremI));
bubblesize([2 12]);
%bubblesize(15*[min(meanLinearity2(tremI)) max(meanLinearity2(tremI))]);
SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.75;
SS.MarkerEdgeAlpha = 0.75;
SS.LineWidth = 1;

blgd = bubblelegend('linearity','Location','northwest');
%blgd.LimitLabels = {num2str(round(min(eqmag)*100)/100),num2str(round(max(eqmag)*100)/100)};
blgd.Style = 'horizontal';
blgd.NumBubbles = 3;

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

colormap(ax(laxes),'parula');
cbar2 = colorbar(ax(laxes));
cbar2.Position = [0.83 0.115 0.0065 0.81];
%clim(ax(laxes),[1050 1150]);
clim(ax(laxes),[900 1200]);


ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

%cd ~/FigCatchAll/
%print -djpeg CotopaxiArray_initialLocations

%%
figure('units','normalized','outerposition',[0 0 1 1]);
ax__ = subplot(211);
plot(tCommon(tremI),tremLon(tremI),'.'); zoom on; grid on; hold on;
ylabel('lon.');
ax__(2,1) = subplot(212);
plot(tCommon(tremI),tremLat(tremI),'.'); zoom on; grid on; hold on;
ylabel('lat.');
linkaxes(ax__,'x');

figure('units','normalized','outerposition',[0 0 1 1]);
ax2__ = subplot(211);
plot(tCommon(tremI),[doaTAM(tremI) doaCTX2(tremI)],'.'); zoom on; grid on; hold on;
ylabel('doa');
ax2__(2,1) = subplot(212);
plot(tCommon(tremI),[avTAM(tremI) avCTX2(tremI)],'.'); zoom on; grid on; hold on;
ylabel('app. vel.');
linkaxes(ax2__,'x');

%%
if exist('tBell','var')
    figure(); 
    clear axCmp; 
    axCmp(1) = subplot(221); 
    plot(tBell,belDOA(:,2),'.'); zoom on; grid on; hold on; plot(tCommon(tremI),doaCTX2(tremI),'.'); 
    ylabel('baz'); 
    legend('bell','hrndz'); 
    
    axCmp(2) = subplot(222); 
    plot(tBell,bellAV,'.'); zoom on; grid on; hold on; 
    plot(tCommon(tremI),avCTX2(tremI),'.'); 
    ylabel('app. vel.');
    
    axCmp(3) = subplot(223); 
    plot(tBell,bellLon,'.'); zoom on; grid on; hold on; 
    plot(tCommon(tremI),tremLon(tremI),'.'); 
    ylabel('lon'); legend('bell','hrndz'); 
    
    axCmp(4) = subplot(224); 
    plot(tBell,bellLat,'.'); zoom on; grid on; hold on; 
    plot(tCommon(tremI),tremLat(tremI),'.'); 
    ylabel('lat'); 
    linkaxes(axCmp,'x');
end

