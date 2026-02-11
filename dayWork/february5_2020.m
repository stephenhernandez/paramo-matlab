% february5_2020
% this code loads tiffs that pedro espin gave me and overlays topography 
% contours so that i can see the relation of the deformation to elevation
% patterns (and location relative to caldera)

cd ~/research/now/fernandina/

clear; close all; clc;
[fig,boundaryBox] = loadBasemap('fernandina');
fig.Visible = 'on';
tif_file = 'cm_28_03FEBDESCE.tif';
% data = geotiffread(tif_file);
% proj = geotiffinfo(tif_file);
% r = proj.SpatialRef.RasterSize(1);
% c = proj.SpatialRef.RasterSize(2);
% [ROW,COL] = meshgrid(1:r,1:c);
% %[X,Y] = pix2map(proj.RefMatrix,ROW,COL);
% [X,Y] = intrinsicToWorld(proj.RefMatrix,ROW,COL);
% lon = X(:,1);
% lat = Y(1,:)';

[data,R] = readgeoraster(tif_file);
[lat,lon] = geographicGrid(R);
lat = lat(:,1);
lon = lon(1,:)';
figure(); h = imagesc(lon,lat,data); axis xy; axis equal;
c = colorbar;
colormap parula
axis equal;
zoom on;
figure(1); ax1 = gca;
ax1Chil = ax1.Children;
figure(2);
ax = gca;
ax(2) = axes;
newHandle = copyobj(ax1Chil,ax(2));
c2 = colorbar(ax(2));
c2.Visible = 'off';
ax(2).Visible = 'off';
axis(ax(1),'equal');
axis(ax(2),'equal');
axis(ax(1),boundaryBox);
axis(ax(2),boundaryBox);
linkaxes(ax,'xy');
colormap(ax(2),'bone');
caxis(ax(1),5*[-1 1]);
