% esmeraldas2022 from different velocity models
clear; close all; clc;

cd ~/research/now/esmeraldas/hypodd/
models = ["esmeraldas2022.m3d_2.reloc";...
    "esmeraldas2022.mleon.reloc";...
    "esmeraldas2022.leon.reloc";...
    "esmeraldas2022.central.reloc";...
    "esmeraldas2022.macquet.reloc";...
    "esmeraldas2022.soto.reloc";...
    "esmeraldas2022.mfont.reloc"];

titleStrs = ["M3D_2";"MLEON";"LEON";"CENTRAL";"MACQUET";"SOTO";"MFONT"];
for i = 1:length(models)+1
    if i == length(models)+1
        load orig_cat.txt;
        t = datetime(orig_cat(:,2:7));
        eqlat = orig_cat(:,8);
        eqlon = orig_cat(:,9);
        eqdepth = orig_cat(:,10);
        eqmag = orig_cat(:,11);
        plot_fig_(eqlat,eqlon,eqdepth,eqmag,sprintf("%s, N=%d\n","IASPEI",length(t)));
        continue;
    end

    [t,eqlat,eqlon,eqdepth,eqmag,id] = readHypoDD(models(i),true);
    [eqmag,sI] = sort(eqmag,'descend');
    t = t(sI);
    eqlat = eqlat(sI);
    eqlon = eqlon(sI);
    eqdepth = eqdepth(sI);
    id = id(sI);

    refVMdepth = 0;
    eqdepth = eqdepth-refVMdepth;
    plot_fig_(eqlat,eqlon,eqdepth,eqmag,sprintf("%s, N=%d",titleStrs(i),length(t)));
end

function plot_fig_(eqlat,eqlon,eqdepth,eqmag,titleStr)
[minLon,maxLon,minLat,maxLat] = deal(-79.9,-79.65,0.7,1.05);
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);

contrast = 1/4;
ilumAz = 90;
wMark = 3;
elevAngle = 50;

figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 0.5 1]);
ax(laxes) = axes; %subplot(8,1,1:7);
ax(laxes).Parent = fig(figNumber);

I = dem(lon,lat,demData','Contrast',contrast,'Azimuth',ilumAz,'noplot','Watermark',wMark,'Elevation',elevAngle,'NaNColor',[1 1 1]);
I = I.rgb;
imagesc(ax(laxes),lon,lat,I);

cbar1 = colorbar(ax(laxes));
colormap(cbar1);
clim([min(min(demData)) max(max(demData))]);
cbar1.Visible = 'off';

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);

stationOutlineColor = 'k';
plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
plot(ax(laxes),lon_noec,lat_noec,'-','linewidth',2,"Color",[0.5 0.5 0.5]);
grid on;
zoom on;

%% scaleruler
refShifter = (maxLat-minLat)*(0.007/0.3);
scaleRefLat0 = minLat + refShifter;
scaleRefLon0 = minLon + refShifter;
refEllipse = referenceEllipsoid('wgs84');

[~,lonout1] = reckon(scaleRefLat0,scaleRefLon0,1*1e3,90,refEllipse);
[~,lonout5] = reckon(scaleRefLat0,scaleRefLon0,5*1e3,90,refEllipse);

scaleRefLat2 = scaleRefLat0 + refShifter/2;
plot(ax(laxes),[scaleRefLon0 lonout5],[scaleRefLat0 scaleRefLat0],'k-','linewidth',2);

textShifter = -refShifter/2;
textFontSize = 18;
plot(ax(laxes),[scaleRefLon0 scaleRefLon0],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+scaleRefLon0,-textShifter+scaleRefLat2,"0","FontSize",textFontSize);

plot(ax(laxes),[lonout1 lonout1],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout1,-textShifter+scaleRefLat2,"1","FontSize",textFontSize);

plot(ax(laxes),[lonout5 lonout5],[scaleRefLat0 scaleRefLat2],'k-','linewidth',2);
text(textShifter+lonout5,-textShifter+scaleRefLat2,"5 km","FontSize",textFontSize);

%%
labelFontSize = 24;
ylabel(ax(laxes),'Latitude','FontSize',labelFontSize);
xlabel(ax(laxes),'Longitude','FontSize',labelFontSize);
title(ax(laxes),titleStr,'FontSize',labelFontSize)
magFact = 14;

%%
laxes = laxes + 1;
ax(laxes) = axes;

SS = bubblechart(ax(laxes),eqlon,eqlat,magFact*exp(eqmag),eqdepth); %,'filled');
SS.MarkerEdgeColor = 'k';
SS.MarkerFaceAlpha = 0.8;
SS.MarkerEdgeAlpha = 0.8;
SS.LineWidth = 1;
blgd = bubblelegend('Magnitude','Location','northwest');
blgd.LimitLabels = {num2str(round(10*min(eqmag))/10),num2str(round(10*max(eqmag))/10)};
blgd.Style = 'vertical';
blgd.NumBubbles = 2;

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

colormap turbo;
cbar2 = colorbar(ax(laxes));

clim([14 20]);
cbar2.Label.String = 'Depth [km]';

StationMarkerSize = 20;
[stlaH,stloH] = metaDataFromStationList(["HB10";"HB12";"HB15";"HB17";"HB18";"HB19";"HB56"]);
h1 = plot(ax(laxes),stloH,stlaH,'v','markeredgecolor','k','markerfacecolor','w','MarkerSize',StationMarkerSize);

allStationKstnms = ["AAT1","AATC","AES1","AES2","ALIT","ALOR","ANGU",...
    "ANTG","ANTI","ANTM","ANTS","AOTA","APED","ARA2","ASAM","ASDO","ATON",...
    "BBAC","BMAS","BMOR","BNAS","BONI","BPAT","BREF","BTAM","BTER","BVC2",...
    "CAB1","CAYR","CERN","CHL1","CHL2","CHMA","CHSH","CIVI","COTA","CUIC",...
    "CUSE","CUSW","ECEN","ESM1","FLF1","GGPC","GGPT","ILLI","IMBA",...
    "JAMA","JIPI","JSCH","JUA2","LITA","LNGL","MONB","NAS2","OTAV","PAC1",...
    "PAS1","PECV","PINO","PIS1","PITA","PORT","PTGL","PUEM","PULU","REVS",...
    "SLOR","SNLR","SRAM","SUCR","TAMB","TAMH","TING","TOMA","URCU",...
    "VC1","VCES","YAHU","YANA"];

[stla,stlo] = metaDataFromStationList(allStationKstnms);
h2 = plot(ax(laxes),stlo,stla,'^','markeredgecolor','w','markerfacecolor','k','MarkerSize',StationMarkerSize);

GPS_STLA(1) = 0.78153;
GPS_STLO(1) = -80.03036;
GPS_STNM(1) = "PTGL";

GPS_STLA(2,1) = 0.86165;
GPS_STLO(2,1) = -79.85933;
GPS_STNM(2,1) = "TRAT";

GPS_STLA(3,1) = 0.97878;
GPS_STLO(3,1) = -79.62541;
GPS_STNM(3,1) = "TRTH";

h3 = scatter(ax(laxes),GPS_STLO,GPS_STLA,StationMarkerSize.^2,"magenta",'s','filled'); %%,'markeredgecolor','k','markerfacecolor','m','MarkerSize',StationMarkerSize);
h3.MarkerEdgeColor = 'k';
h3.MarkerFaceAlpha = 0.85;
[~,legIcons] = legend([h1;h2;h3],'Seismic (Temp.)','Seismic (Perm.)','GPS','Location','NorthWest');
legIcons = findobj(legIcons, '-not', 'Marker', 'none');
legIcons = findobj(legIcons, '-property', 'Marker', '-and', '-not', 'Marker', 'none');
legIcons(4).MarkerSize = StationMarkerSize;

linkaxes(ax);
ax(laxes).Visible = 'off';
axis(ax(laxes),[minLon maxLon minLat maxLat]);

%%
laxes = laxes + 1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

plot(ax(laxes),NaN,NaN,'.');

cbar4 = colorbar(ax(laxes));
cbar4.Visible = 'off';
cbar4.Position = [0.9,0.01,0.0083333333333333,0.8];

hold(ax(laxes),'on');
axis(ax(laxes),'equal');

ax(laxes).Visible = 'off';
linkaxes(ax);
axis(ax(laxes),[minLon maxLon minLat maxLat]);
titleSplit = split(titleStr,",");
print('-djpeg',strcat("Figure_",titleSplit(1),"_2022"));
end
