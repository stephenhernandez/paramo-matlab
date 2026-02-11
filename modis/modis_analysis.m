%modis analysis
%this code takes 15 minutes to run, its too slow, find a way to speedup
clear; close all; clc;

cd ~/masa/modis_burnt_area/
lats = []; lons = []; t = []; qa = []; first_days = []; last_days = []; uncert = [];
files = dir('*.hdf');
lFiles = length(files);
for i = 1:lFiles
    %[lats_,lons_,t_,qa_,first_days_,last_days_,uncert_] = read_modis_hdf_file(modis_hdf_file);
    modis_hdf_file = files(i).name;
    disp(modis_hdf_file);
    [lats_,lons_,t_,qa_,first_days_,last_days_,uncert_] = read_modis_hdf_file(modis_hdf_file);
    lats = [lats; lats_];
    lons = [lons; lons_];
    t = [t; t_];
    qa = [qa; qa_];
    first_days = [first_days; first_days_];
    last_days = [last_days; last_days_];
    uncert = [uncert; uncert_];
end

%%
[t,sI] = sort(t);
lats = lats(sI);
lons = lons(sI);
qa = qa(sI);
first_days = first_days(sI);
last_days = last_days(sI);
uncert = uncert(sI);
load ~/igdata/ec_boundaries.mat;
in = inpolygon(lons,lats,lonEC,latEC);
xy = [lons(in) lats(in)];
[C,ia,ic] = unique(xy,'rows');
close all
figure(); plot(C(:,1),C(:,2),'.'); zoom on; grid on; axis equal;
axis equal
load ~/igdata/ec_boundaries.mat
hold on; plot(lonEC,latEC,'k-','linewidth',2); axis equal;
tEC = t(in);
fdEC = first_days(in);
ldEC = last_days(in);
figure(); plot(tEC,days(ldEC-fdEC),'.'); zoom on;
figure(); plot(tEC,(0:length(tEC)-1)','.'); zoom on; grid on;

nDays = 365;
rate = t2r(tEC,days(nDays)); figure(); plot(tEC,rate/nDays,'.'); zoom on; grid on;