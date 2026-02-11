%clear; close all; clc;
refEllipse = referenceEllipsoid('wgs84');
maplegend2(1) = 1081;
maplegend2(2) = 1;
maplegend2(3) = -79.9;

minLon = -79.9; maxLon = -79.6; minLat = 0.7; maxLat = 1;
boundaryBox = [minLon; maxLon; minLat; maxLat];

plat1 = 0.91;
plon1 = -79.9;
plat2 = 0.73;
plon2 = -79.6;
lslices = length(plat1);

[DEMlon,DEMlat,demData] = cutDEM(boundaryBox,true);
[r,c] = meshgrid(DEMlon,DEMlat);
r = r'; c = c';
F = griddedInterpolant(r,c,demData,"spline","nearest");

geolat = [];
geolon = geolat;
[hw,az_] = distance(plat1(1),plon1(1),plat2(1),plon2(1),refEllipse);
latout = plat1(1);
lonout = plon1(1);
nRange = 2e3;
for j = 1:nRange
    [latout,lonout] = reckon(latout,lonout,hw/nRange,az_,refEllipse);
    geolat = [geolat; latout];
    geolon = [geolon; lonout];
    [~,az_] = distance(latout,lonout,plat2(1),plon2(1),refEllipse);
end

z = F((geolon),(geolat));
magFact = 8;
rang = hw*(0:nRange-1)*1e-3/nRange;

%%
for i = 1:lslices
    [~,~,lat_,lon_] = mapprofile(demData,maplegend2,[plat1(i) plat2(i)],[plon1(i) plon2(i)],refEllipse);
    [~,az_] = distance(plat1(i),plon1(i),plat2(i),plon2(i),refEllipse);

    az1 = az_ + 90;
    if az1 >= 360
        az1 = az1 - 360;
    end

    az2 = az_ - 90;
    if az2 >= 360
        az2 = az2 - 360;
    end

    [yv(1),xv(1)] = reckon(plat2(i),plon2(i),hw,az1,refEllipse);
    [yv(2),xv(2)] = reckon(plat2(i),plon2(i),hw,az2,refEllipse);
    [yv(3),xv(3)] = reckon(plat1(i),plon1(i),hw,az2,refEllipse);
    [yv(4),xv(4)] = reckon(plat1(i),plon1(i),hw,az1,refEllipse);

    figure(2); hold on;
    plot(lon_,lat_,'r');
    plot([xv(:);xv(1)],[yv(:);yv(1)],'-','linewidth',2);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(rang,z*1e-3,'-','linewidth',2);
    xlim([0 max(rang)]);
    max(rang)
    [~,sI] = sort(eqmag,'descend');
    in = inpolygon(eqlon,eqlat,xv,yv);
    in = sI(in);

    [erot,nrot] = rotate2d(eqlon-plon1(1),eqlat-plat1(1),az_);
    hold on;
    SS = scatter(111.19*nrot(in),-eqdepth(in),magFact*exp(eqmag(in)),...
        -eqdepth((in)),'filled','linewidth',1,'markeredgecolor','k');
    SS.MarkerFaceAlpha = 0.8;
    zoom on;
    ylim([-35 0.5]);
    colorbar;
    clim([-21 -16]);
    grid on;
end
