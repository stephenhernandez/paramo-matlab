cd ~/research/now/imantag/
clear; close all; clc;
magFact = 6;
refVMDepth = 0;
tStart = datetime(2016,01,01);
tEnd = datetime(2016,06,01);
minLon = -90; maxLon = -70; minLat = -10; maxLat = 10;

oldData = load('orig_cat.txt');
origID = oldData(:,1);
told = dn2dt(datenum(oldData(:,2:7)));
origLat = oldData(:,8);
origLon = oldData(:,9);
origDepth = oldData(:,10);
origMag = oldData(:,11);

%%
figure(); 
ax(1) = subplot(121);
SOld = scatter(origLon,origLat,magFact*exp(origMag),refVMDepth-origDepth,'filled','linewidth',0.02,'markeredgecolor','k');
c = colorbar;
caxis([-5 refVMDepth]);
c.Label.String = 'Depth [km]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.Label.FontSize = 18;


%%
oldData = load('reloc_cat.txt');
newID = oldData(:,1);
newt = dn2dt(datenum(oldData(:,2:7)));
newLat = oldData(:,8);
newLon = oldData(:,9);
newDepth = oldData(:,10);
newMag = oldData(:,11);

ax(2) = subplot(122);
scatter(newLon,newLat,magFact*exp(newMag),refVMDepth-newDepth,'filled','linewidth',0.02,'markeredgecolor','k'); zoom on;
c = colorbar;
caxis([-4 refVMDepth]);
c.Label.String = 'Depth [km]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.Label.FontSize = 18;
%axis 'equal';
linkaxes(ax,'xy');