%clear;
%cd ~/regions_old/pedernales/coulomb/
cd ~/research/now/pedernales_coulomb/
fNames = dir('ss_dc3d_xy_*');

%
load('~/igdata/ec_boundaries','latEC','lonEC');
load('~/igdata/soam_noec');
load('~/igdata/ecTrench','latTrench','lonTrench');
S = shaperead('~/igdata/volcanes/volcanes_2.shp');
[center_lat,center_lon] = get_region_dimensions('pichincha');
minLon = -82;
maxLon = -78;
minLat = -2;
maxLat = 2;

%
maxDepthI = 60;
fNames = fNames(1:maxDepthI);
printFlag = false;
startLonI1 = 190;
startLatI = 102;
depthSliceI = 3;

%
hx = 2000;
hy = hx;
hz = 1000;

%
load(fNames(1).name,'DC3D','XGRID','YGRID');
ZGRID = -(0:maxDepthI-1)';
YGRID = YGRID';
XGRID = XGRID';

%
lX = length(XGRID);
lY = length(YGRID);
lZ = length(ZGRID);
lL = lX*lY;

%
UXbig = NaN(lY,lX,maxDepthI);
UYbig = UXbig;
UZbig = UXbig;

%
UXX = UXbig;
UYY = UXbig;
UZZ = UXbig;
UYZ = UXbig;
UXZ = UXbig;
UXY = UXbig;

%
shears = UXbig;
normals = UXbig;
coulombs = UXbig;

%
strike_m = 27;
dip_m = 23;
rake_m = 118;
friction_m = 0.4;

satLev = 1e-6;
for i = 1:length(fNames)
    load(fNames(i).name,'DC3D');
    
    % dependent variables
    UXbig(:,:,i) = reshape(DC3D(:,6),[lY lX]);
    UYbig(:,:,i) = reshape(DC3D(:,7),[lY lX]);
    UZbig(:,:,i) = reshape(DC3D(:,8),[lY lX]);
    
    % strains
    UXX(:,:,i) = satLev*reshape(DC3D(:,9),[lY lX]);
    UYY(:,:,i) = satLev*reshape(DC3D(:,10),[lY lX]);
    UZZ(:,:,i) = satLev*reshape(DC3D(:,11),[lY lX]);
    
    % more strains
    UYZ(:,:,i) = satLev*reshape(DC3D(:,12),[lY lX]);
    UXZ(:,:,i) = satLev*reshape(DC3D(:,13),[lY lX]);
    UXY(:,:,i) = satLev*reshape(DC3D(:,14),[lY lX]);
    
    n = length(DC3D);
    c1 = zeros(n,1) + strike_m;
    c2 = zeros(n,1) + dip_m;
    c3 = zeros(n,1) + rake_m;
    c4 = zeros(n,1) + friction_m;
    
    tic;
    [shear_,normal_,coulomb_] = calc_coulomb(c1,c2,c3,c4,DC3D(:,9:14)');
    toc;
    
    shears(:,:,i) = reshape(1e5*shear_,[lY lX]);
    normals(:,:,i) = reshape(1e5*normal_,[lY lX]);
    coulombs(:,:,i) = reshape(1e5*coulomb_,[lY lX]);
end
clear DC3D

% UXX = gradient(UXbig,hx);
% [~,UYY] = gradient(UYbig,hy);
% [~,~,UZZ] = gradient(UZbig,hz);

%%
coulombs = 1e-3*coulombs; %convert to kPa
shears = 1e-3*shears; %convert to kPa
normals = 1e-3*normals; %convert to kPa
strainTrace = (UXX+UYY+UZZ)/satLev;

%
satLev = 1;
startLonI = 1+lY*(startLonI1-1);

%
[XX,YY] = meshgrid(XGRID,YGRID,ZGRID);
XX_ = squeeze(XX(:,:,1)); XX_ = XX_(:);
YY_ = squeeze(YY(:,:,1)); YY_ = YY_(:);

%
hI = (startLatI:lY:lL)';
refEllipse = referenceEllipsoid('wgs84');

%%
[lonRecovered,latRecovered] = xy2latlon(XX_,YY_,-80,0,refEllipse);
latVec = reshape(latRecovered,[lY lX]);
lonVec = reshape(lonRecovered,[lY lX]);
latVec = latVec(:,1);     % this is vertical (Y) component
lonVec = lonVec(1,:)';    % this is horizontal (X) component

%
lonRecovered = reshape(lonRecovered,[lY lX]);
latRecovered = reshape(latRecovered,[lY lX]);
lonRecovered = repmat(lonRecovered,[1 1 lZ]);
latRecovered = repmat(latRecovered,[1 1 lZ]);

% close all;
cmp1 = squeeze(UXbig(startLatI,:,:))';
cmp2 = squeeze(UZbig(startLatI,:,:))';

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
imagesc(XGRID,ZGRID,cmp1);
hold on;
axis xy;
colorbar;
caxis([-2 2]);
ylabel('Depth [km]');
title('Horizontal Displacements');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

subplot(2,2,2);
imagesc(XGRID,ZGRID,cmp2);
hold on;
axis xy;
colorbar;
caxis([-2 2]);
ylabel('Depth [km]');
title('Vertical Displacements');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

fig2 = figure('units','normalized','outerposition',[0 0 1 1]); % figure 2
strain_ = squeeze(strainTrace(startLatI,:,:))';
subplot(2,1,1);
imagesc(XGRID,ZGRID,strain_);
hold on;
contour(XGRID',ZGRID,strain_,(0.1:0.1:0.9)','linewidth',1,'color',[0.5 0.5 0.5]);
contour(XGRID',ZGRID,strain_,[1 1],'linewidth',2,'color',[0.5 0.5 0.5]);
contour(XGRID',ZGRID,strain_,[0 0],'-.','linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
axis equal;
%colorbar;
caxis(satLev*[-1 1]);
ylabel('Depth [km]');
xlabel('Distance from -80$^\circ$ Longitude [km]');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

%
cmp1 = squeeze(UYbig(:,startLonI1,:))';
cmp2 = squeeze(UZbig(:,startLonI1,:))';

figure(fig1); %figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,3);
imagesc(YGRID,ZGRID,cmp1);
hold on;
contour(YGRID',ZGRID,cmp1,[0 0],'linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
colorbar;
caxis([-0.05 0.05]);
ylabel('Depth [km]');
title('Horizontal Displacements');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

subplot(2,2,4);
imagesc(YGRID,ZGRID,cmp2);
hold on;
contour(YGRID',ZGRID,cmp2,[0 0],'linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
colorbar;
caxis([-0.05 0.05]);
ylabel('Depth [km]');
title('Vertical Displacements');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

figure(fig2);
strain_ = squeeze(strainTrace(:,startLonI1,:))';
subplot(2,1,2);
imagesc(YGRID,ZGRID,strain_);
hold on;
contour(YGRID',ZGRID,strain_,(0.1:0.1:0.9)','linewidth',1,'color',[0.5 0.5 0.5]);
contour(YGRID',ZGRID,strain_,[1 1],'linewidth',2,'color',[0.5 0.5 0.5]);
contour(YGRID',ZGRID,strain_,[0 0],'-.','linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
axis equal;
%colorbar;
caxis(satLev*[-1 1]);
ylabel('Depth [km]');
xlabel('Distance from 0$^\circ$ Latitude [km]');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

%
cmp1 = squeeze(UXbig(:,:,depthSliceI));
cmp2 = squeeze(UYbig(:,:,depthSliceI));
cmp3 = squeeze(UZbig(:,:,depthSliceI));

%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1);
imagesc(lonVec,latVec,cmp1);
hold on;
contour(lonVec',latVec,cmp1,[0 0],'linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
axis equal;
colorbar;
caxis(0.5*[-1 1]);
plot(center_lon,center_lat,'^');
plot(lonEC,latEC,'k','linewidth',1);
plot([-82 -78.5],[center_lat center_lat],'--','color',[0.5 0.5 0.5],'linewidth',1);
plot([center_lon center_lon],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',1);
plot([S.X],[S.Y],'k--','linewidth',1);
axis([minLon maxLon minLat maxLat]);
title('East-West Displacement');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

subplot(1,3,2);
imagesc(lonVec',latVec,cmp2);
hold on;
contour(lonVec',latVec,cmp2,[0 0],'linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
axis equal;
colorbar;
caxis(0.5*[-1 1]);
plot(center_lon,center_lat,'^');
plot(lonEC,latEC,'k','linewidth',1);
plot([-82 -78.5],[center_lat center_lat],'--','color',[0.5 0.5 0.5],'linewidth',1);
plot([center_lon center_lon],[-2 2],'--','color',[0.5 0.5 0.5],'linewidth',1);
plot([S.X],[S.Y],'k--','linewidth',1);
axis([minLon maxLon minLat maxLat]);
title('North-South Displacement');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

subplot(1,3,3);
imagesc(lonVec',latVec,cmp3);
hold on;
contour(lonVec',latVec,cmp3,[0 0],'linewidth',1,'color',[0.5 0.5 0.5]);
axis xy;
axis equal;
colorbar;
caxis(0.5*[-1 1]);
plot(center_lon,center_lat,'^');
plot(lonEC,latEC,'k','linewidth',1);
plot([minLon maxLon],[center_lat center_lat],'--','color',[0.5 0.5 0.5],'linewidth',1);
plot([center_lon center_lon],[minLat maxLat],'--','color',[0.5 0.5 0.5],'linewidth',1);
axis([minLon maxLon minLat maxLat]);
title('Vertical Displacement');
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;

%
figure('units','normalized','outerposition',[0 0 1 1]);
strain_ = squeeze(strainTrace(:,:,depthSliceI));
imagesc(lonVec,latVec,strain_);
hold on;
axis xy;
axis equal;
c = colorbar;
caxis(satLev*[-1 1]);
hh = plot([center_lon center_lon],[minLat maxLat],'-.','linewidth',1);
plot([minLon maxLon],[center_lat center_lat],'-.','linewidth',1,'Color',hh.Color);
plot(center_lon,center_lat,'^','markeredgecolor','k','markerfacecolor','w');
contour(lonVec',latVec,strain_,(0.1:0.1:0.9)','linewidth',1,'color',[0.5 0.5 0.5]); %,'ShowText','on');
contour(lonVec',latVec,strain_,[1 1],'linewidth',2,'color',[0.5 0.5 0.5]);
contour(lonVec',latVec,strain_,[0 0],'-.','linewidth',1,'color',[0.5 0.5 0.5]);
plot(lonTrench,latTrench,'--','Color','k','linewidth',1);
plot([S.X],[S.Y],'k--','linewidth',1);
c.Label.String = 'microstrain [$\mu\epsilon$]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = 18;
plot(lonEC,latEC,'k','linewidth',2);
plot(lon_noec,lat_noec,'k','linewidth',1);
ax = gca;
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
axis([minLon maxLon minLat maxLat]);
ylabel('Latitude');
xlabel('Longitude');

%%
load('~/igdata/ecSlabModelUSGS');
figNumber = 5;
[XXLon,YYLat,ZZ] = meshgrid(lonVec,latVec,-(0:59)');
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]); 
ax = gca; 
s = slice(ax,XXLon,YYLat,ZZ,coulombs,xq,yq,newDepth);
shading interp
c = colorbar;
c.Label.String = '$\delta$CFF [kPa]';
c.Label.Interpreter = 'latex';

hold on;
load('~/igdata/ec_boundaries','latEC','lonEC');
plot3(ax,lonEC,latEC,zeros(size(lonEC)),'k','linewidth',4);
caxis([-100 100]);

view(0,90);
axis equal;
xlim([-82 -78]);
ylim([-2 2]);
hold on;

newX = s.XData;
newY = s.YData;
newZ = s.ZData;
newC = s.CData;

contour(ax,newX,newY,newC,(-100:10:100)','linewidth',2,'color',[0.5 0.5 0.5],'showtext','on');
contour(ax,newX,newY,newZ,-(0:10:80)','k--','linewidth',1);

axis equal;
xlim([-82 -78]);
ylim([-2 2]);
ylabel('Latitude');
xlabel('Longitude');

disp('contouring coseismic slip');
data = importdata('~/regions_old//pedernales/slip_53_s_m.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
disp('done contouring coseismic slip');

laxes = length(ax);
laxes = laxes +1;
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);
contour(ax(laxes),xq,yq,vq,1:7,'linewidth',2);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
colormap(ax(laxes),'hot');
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';
linkaxes(ax(1:laxes),'xy');
axis(ax(1:laxes),'equal');
axis(ax(1:laxes),[minLon maxLon minLat maxLat]);

%%
if printFlag
    for i = 1:4
        figure(i);
        fname = strcat('strainModelingFigure_',num2str(i));
        print('-depsc',fname);
    end
end
