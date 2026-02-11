clear
close all; 

data = importdata('~/igdata/wsm_020.0_regularization_Corrm_exponential_030_km_sm0_000_bayesian_sol_coupling.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
figure(); contour(xq,yq,vq,50:10:90,'linewidth',2);

cd ~/research/now/pedernales/
data = readtable('isc_nas_model.gmt','FileType','text');

isclon = data.x_lon;
isclat = data.lat;
isc = data.ISC_;

[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);

hold on;
contour(xq,yq,vq,50:10:90,'linewidth',2);
load ~/igdata/ec_boundaries.mat
hold on; plot(lonEC,latEC,'k-','linewidth',2); axis equal; zoom on; grid on;