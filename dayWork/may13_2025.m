%may13_2025
clear; close all; clc;
cd ~/Desktop/
Tdeep = readtable("displacement_esmeraldas_deep_2.cou","FileType","text");
xgrid = Tdeep.Var1;
ygrid = Tdeep.Var2;

yinc = median(diff(ygrid(1:10)));
xinc = yinc;

loninc = 0.01;
latinc = loninc;

minLon = -80.3;
minLat = 0.6;

LON_PER_X = loninc / xinc;
LAT_PER_Y = latinc / yinc;

grid = [min(xgrid);min(ygrid);max(xgrid);max(ygrid)];

dx = xgrid - grid(1);
dy = ygrid - grid(2);

lonNew = minLon + dx.*LON_PER_X;
latNew = minLat + dy.*LAT_PER_Y;

load ~/igdata/ec_boundaries.mat;

close all;
figure();
plot(lonEC,latEC,'k','linewidth',4);
hold on; ax = gca;
ax.ColorOrderIndex = 1;
Qdeep = quiver(lonNew,latNew,Tdeep.Var4,Tdeep.Var5,0,"LineWidth",2);
zoom on; axis equal;

%% fault patch
%strike = 36;
dip = 25;
rake = 125;

xstart = -8.397163911607818;
ystart = -4.031629747979384;
xfinish = 1.242514225988740;
yfinish = 9.236248959769753;


% %%MW6.2
% xstart = -7.311359261370606;
% ystart = -3.366916303958663;
% xfinish = 0.947023533338642;
% yfinish = 7.999772467009350;

xwidth = xfinish-xstart;
yheight = yfinish-ystart;
strike = 90-atan2d(yheight,xwidth);
fault_length = sqrt(xwidth^2 + yheight^2);
fault_width = 9.77; %8.68;
projected_fault_width_meters = 1e3*fault_width*cosd(dip);

refEllipse = referenceEllipsoid('wgs84');
faultLon = minLon + (xstart-grid(1)) .* LON_PER_X;
faultLon(2,1) = minLon + (xfinish-grid(1)) .* LON_PER_X;
faultLat = minLat + (ystart-grid(2)) .* LAT_PER_Y;
faultLat(2,1) = minLat + (yfinish-grid(2)) .* LAT_PER_Y;

[lat_,lon_] = reckon(faultLat(2),faultLon(2),projected_fault_width_meters,90+strike,refEllipse);
faultLat(3,1) = lat_;
faultLon(3,1) = lon_;


[lat_,lon_] = reckon(faultLat(1),faultLon(1),projected_fault_width_meters,90+strike,refEllipse);
faultLat(4,1) = lat_;
faultLon(4,1) = lon_;

%%
plot([faultLon;faultLon(1)],[faultLat;faultLat(1)],'linewidth',2);
epicenterLon = -79.69034;
epicenterLat = 1.09784;
ESMR_lon = -79.724375;
ESMR_lat = 0.934665;
plot(epicenterLon,epicenterLat,'kp',"LineWidth",2);
plot(ESMR_lon,ESMR_lat,'kv','MarkerFaceColor','w');
text(ESMR_lon,ESMR_lat,"ESMR",'FontSize',14);

%%
ew_deep_orig = Tdeep.Var4;
ns_deep_orig = Tdeep.Var5;
F_ew_deep = scatteredInterpolant(lonNew,latNew,ew_deep_orig,"linear","nearest");
esmr_ew_deep_predicted = F_ew_deep(ESMR_lon,ESMR_lat);
F_ns_deep = scatteredInterpolant(lonNew,latNew,ns_deep_orig,"linear","nearest");
esmr_ns_deep_predicted = F_ns_deep(ESMR_lon,ESMR_lat);

%%
Tshallow = readtable("displacement_esmeraldas_shallow_2.cou","FileType","text");
xgrid = Tshallow.Var1;
ygrid = Tshallow.Var2;

yinc = median(diff(ygrid(1:10)));
xinc = yinc;
LON_PER_X = loninc / xinc;
LAT_PER_Y = latinc / yinc;
grid = [min(xgrid);min(ygrid);max(xgrid);max(ygrid)];

dx = xgrid - grid(1);
dy = ygrid - grid(2);
lonNew = minLon + dx.*LON_PER_X;
latNew = minLat + dy.*LAT_PER_Y;

ax.ColorOrderIndex = 2;
Qshallow = quiver(lonNew,latNew,Tshallow.Var4,Tshallow.Var5,0,"LineWidth",1);

scale=1;
hU1 = get(Qdeep,'UData');
hV1 = get(Qdeep,'VData');
set(Qdeep,'UData',scale*hU1,'VData',scale*hV1)
hU2 = get(Qshallow,'UData');
hV2 = get(Qshallow,'VData');
set(Qshallow,'UData',scale*hU2,'VData',scale*hV2);

%% interpolate to ESMR site
ew_shallow_orig = Tshallow.Var4;
ns_shallow_orig = Tshallow.Var5;
F = scatteredInterpolant(lonNew,latNew,ew_shallow_orig,"linear","nearest");
esmr_ew_shallow_predicted = F(ESMR_lon,ESMR_lat);
F = scatteredInterpolant(lonNew,latNew,ns_shallow_orig,"linear","nearest");
esmr_ns_shallow_predicted = F(ESMR_lon,ESMR_lat);