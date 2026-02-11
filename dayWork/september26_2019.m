clear;
close all;
clc;

%%
[stla,stlo,stel] = metaDataFromStationList(["FER1" "FER2"]);
lo = -91.4500605 - 0.05;
la = -0.3797075;

%%
refEllipse = referenceEllipsoid('wgs84');
[t,eqlat,eqlon,eqdepth,eqmag,ids] = readCat1();

%%
d_ = distance(eqlat,eqlon,la,lo,refEllipse)*1e-3;
dI = d_ < 40;

%%
figure('units','normalized','outerposition',[0 0 1 0.75]);
subplot(211)
plot(stlo,stla,'^');
axis equal
load ~/igdata/ec_boundaries.mat
hold on;
plot(lonEC,latEC,'k','linewidth',2); axis equal; zoom on;
axis([-92 -91 -0.8 0])
S = scatter(eqlon(dI),eqlat(dI),2*exp(eqmag(dI)),eqdepth(dI),'filled'); 
c = colorbar;
S.MarkerEdgeColor = 'k';
S.MarkerFaceAlpha = 0.2;
c.Label.String = 'Depth [km]';
c.Label.Interpreter = 'Latex';
caxis([0 20]);

subplot(212);
S2 = scatter(t(dI),1:length(t(dI)),8*exp(eqmag(dI)),d_(dI),'o','filled'); 
c2 = colorbar;
S2.MarkerEdgeColor = 'k';
S2.MarkerFaceAlpha = 0.2;
c2.Label.String = 'Distance to Reference Point [km]';
c2.Label.Interpreter = 'Latex';