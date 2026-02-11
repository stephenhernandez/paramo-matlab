clear; close all; clc;

%%
cd ~/research/now/sn_eruption/IGUANA_PICKS_20180624_20180627/BULLETINS/
files = dir('*.bulletin');
for i = 1:length(files)
    E(i,1) = readSCBulletin(files(i).name);
end
clear files i
regionName = 'sierra_negra';
outfile = ['~/regions_old/',regionName,'/hypodd/phase_',regionName,'.dat'];
depCorr = 2;

%%
scthree2hypodd();
cd ~/regions_old/sierra_negra/hypodd/
phase_input_file = ['ph2dt_',regionName,'.inp'];
phase_data_file = ['phase_',regionName,'.dat'];

%%
cd(['~/regions_old/',regionName,'/hypodd/']);
unix(['~/soft/hypodd/src/ph2dt/ph2dt ',phase_input_file]);
unix('~/soft/hypodd/src/hypoDD/hypoDD hypoDD.inp');
unix(['bash ~/scripts/separateCatalogs.sh ',phase_data_file]);
unix('wc reloc_cat.txt');

%%
close all; 
[t,eqlat,eqlon,eqdepth,eqmag,id,D,volcanoName] = generateSeismicityAnimation(datetime(2018,06,23),datetime(2018,06,29),regionName,0);
cd ~/regions_old/sierra_negra/animation/
!bash ~/scripts/mkAnimation.sh diasRelocation

%%
close all;
plotOldAndNewLocations();
[t,eqlat,eqlon,eqdepth,eqmag,id,D,volcanoName] = generateSeismicityAnimation(datetime(2018,06,23),datetime(2018,06,29),regionName,0);
t = pull(E,'t');
eqlat = pull(E,'lat');
eqlon = pull(E,'lon');
eqmag = pull(E,'mag');
eqdepth = pull(E,'depth');
figure(); S = scatter(eqlon,eqlat,2*exp(eqmag),eqdepth,'filled'); colorbar; axis equal; caxis([-2 5]);
S.MarkerFaceAlpha = 0.5;
S.MarkerEdgeColor = 'k';
