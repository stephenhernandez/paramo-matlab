clear; close all; clc;

regionName = 'quito';
[~,~,minLat,maxLat,minLon,maxLon] = get_region_dimensions(regionName);
tStart = datetime(2019,05,01);
tEnd = datetime(2019,06,07);
minMag = 1;
maxDepth = 15;
E = extractSC3Catalog(tStart,tEnd,minLon,maxLon,minLat,maxLat,minMag,maxDepth);

workingDir = ['~/research/now/',regionName];
outfile = [workingDir,'/relocation/phase_',regionName,'.dat'];
phase_input_file = ['ph2dt_',regionName,'.inp'];
depCorr = 0;
scthree2hypodd(E,outfile,0);

close all;
cd([workingDir,'/relocation/']);
unix(['~/soft/hypodd/src/ph2dt/ph2dt ',phase_input_file]);
unix('~/soft/hypodd/src/hypoDD/hypoDD hypoDD.inp');
unix(['bash ~/scripts/separateCatalogs.sh ',outfile]);
unix('wc reloc_cat.txt');

%%
plotOldAndNewLocations();

%%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(newLon,-newDepth,'o','linewidth',2); zoom on;
ylim([-15 0]);
xlim([minLon maxLon]);
hold on;
plot(origLon,-origDepth,'s','linewidth',2);
ylabel('Profundidad'); xlabel('Longitud');
legend('New Relocations','Original Locations');

%%
figure('units','normalized','outerposition',[0 0 1 1]); 
plot(newLat,-newDepth,'o','linewidth',2); zoom on;
ylim([-15 0])
xlim([minLat maxLat])
hold on;
plot(origLat,-origDepth,'s','linewidth',2);
ylabel('Profundidad'); xlabel('Latitud');
legend('New Relocations','Original Locations');