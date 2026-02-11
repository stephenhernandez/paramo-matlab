%plotNewCoastAttenuation;

close all;
%newmeanmag = newMag;
id = ids;
%cd ~/regions/pedernales/;
%load('~/igdata/ecuador_coast_dem.mat','lon','lat','dem');
boundaryBox(1) = -84;
boundaryBox(2) = -74;
boundaryBox(3) = -6; 
boundaryBox(4) = 4;

[lon,lat,demdata] = cutDEM(boundaryBox,true);
load('~/igdata/ec_boundaries','latEC','lonEC');
load('~/igdata/ecTrench','latTrench','lonTrench');

markerFaceAlpha = 0.5;
refEllipse = referenceEllipsoid('wgs84');
maplegend2(1) = 743;
maplegend2(2) = 2;
maplegend2(3) = -83;
hw = 15;

magFact = 2;
minTime = dn2dt(floor(min(datenum(t))));
maxTime = dn2dt(ceil(max(datenum(t))));

ilumAz = 180-25;
demdataBlur = imgaussfilt(demdata,2);
I = dem(lon,lat,demdataBlur','Contrast',1,'Azimuth',ilumAz,'Interp','noplot','Watermark',3,'Elevation',55);

figNumber = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(I.x,I.y,I.rgb); axis xy; axis equal;
ax = gca;
hold on;
plot(lonEC,latEC,'k-','linewidth',2);
plot(lonTrench,latTrench,'--','linewidth',2,'color',[0.5 0.5 0.5]);
contour(lon,lat,demdata,-(1000:1000:5000),'k','linewidth',1);

newmeanmag = newMag(tI);
newmeanmag = round(newmeanmag*10)/10;
[magSort,sI] = sort(newmeanmag,'descend');
S = scatter(origlon(sI),origlat(sI),magFact*exp(magSort),origdepth(sI),'filled','linewidth',0.1,'markeredgecolor','k');
colorbar;
clim([0 20]);
S.MarkerFaceAlpha = markerFaceAlpha;

plat1 = [1.4; 1.4; 0.85; 0.5; 0.05; -0.5; -0.4];
plon1 = [-79.7; -80.15; -80.65; -81; -81.2; -80.95; -81.2];

plat2 = [0.75; 0.65; 0.42; 0.14; -0.57; 0.9; 0.05];
plon2 = [-79.45; -79.7; -79.7; -79.75; -80.05; -80.25; -80.05];
lslices = length(plat1);

disp('contouring coseismic slip');
data = importdata('~/research/now/pedernales/slip_53_s_m.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
laxes = length(ax);
disp('done contouring coseismic slip');

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

figure(1);
linkaxes(ax(1:laxes),'xy');
axis(ax(1:laxes),'equal');
%axis(ax(1:laxes),[minLon maxLon minLat maxLat]);
%axis(ax(1:laxes),[-81.3 -79.2 -1 1]);
% axis equal;
axis([-81.3 -79.2 -0.85 1.5]);

%%
for i = 1:lslices
    [z,rang,lat_,lon_] = mapprofile(demdata,maplegend2,[plat1(i) plat2(i)],[plon1(i) plon2(i)],refEllipse);
    rang = rang*1e-3;
    [~,az_] = distance(plat1(i),plon1(i),plat2(i),plon2(i),refEllipse);
    az1 = az_ + 90;
    if az1 >= 360
        az1 = az1 - 360;
    end
    
    az2 = az_ - 90;
    if az2 >= 360
        az2 = az2 - 360;
    end
    
    [yv(1),xv(1)] = reckon(plat2(i),plon2(i),hw*1e3,az1,refEllipse);
    [yv(2),xv(2)] = reckon(plat2(i),plon2(i),hw*1e3,az2,refEllipse);
    [yv(3),xv(3)] = reckon(plat1(i),plon1(i),hw*1e3,az2,refEllipse);
    [yv(4),xv(4)] = reckon(plat1(i),plon1(i),hw*1e3,az1,refEllipse);
    figure(1); hold on;
    plot(lon_,lat_,'r');
    plot([xv(:);xv(1)],[yv(:);yv(1)],'r-','linewidth',2);
    
    %figure i+1
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(rang,z*1e-3,'r','linewidth',2);
    xlim([0 max(rang)]);
    max(rang)
    in = inpolygon(origlon,origlat,xv,yv);
    
%     tStart = datetime(2015,10,01);
%     tEnd = datetime(2018,09,01);
%     regionName = 'pedernales';
%     depthCorrection = 0;
%     refVMDepth = depthCorrection;
%     scBulletinMatName = ['~/regions/',regionName,'/hypodd/',regionName,'.mat'];
%     exportSCBulletin(id(in),scBulletinMatName);
%     outfile = ['~/regions/',regionName,'/hypodd/phase_',regionName,'.dat'];
%     cd(['~/regions/',regionName,'/hypodd/']);
%     scthree2hypodd(regionName,outfile,depthCorrection);
%     phase_input_file = ['ph2dt_',regionName,'.inp'];
%     phase_data_file = ['phase_',regionName,'.dat'];
%     cd(['~/regions/',regionName,'/hypodd/']);
%     unix(['rm ',scBulletinMatName]);
%     unix(['~/soft/hypodd/src/ph2dt/ph2dt ',phase_input_file]);
%     unix('~/soft/hypodd/src/hypoDD/hypoDD hypoDD.inp');
%     unix(['bash ~/scripts/separateCatalogs.sh ',phase_data_file]);
%     unix('wc reloc_cat.txt');
%     
%     oldData = load('reloc_cat.txt');
%     newID = oldData(:,1);
%     newt = dn2dt(datenum(oldData(:,2:7)));
%     newLat = oldData(:,8);
%     newLon = oldData(:,9);
%     newDepth = refVMDepth-oldData(:,10);
%     newMag = oldData(:,11);
%     
%     [~, nrot] = rotate2d(newLon,newLat,az_);
%     hold on;
%     S = scatter(111.19*(nrot-min(nrot)),newDepth,magFact*exp(newMag),datenum(newt),'filled','linewidth',1,'markeredgecolor',[0.5 0.5 0.5]);
    
    [~, nrot] = rotate2d(origlon(in),origlat(in),az_);
    hold on;
    S = scatter(111.19*(nrot-min(nrot)),-origdepth(in),magFact*exp(newmeanmag(in)),datenum(t(in)),'filled','linewidth',1,'markeredgecolor',[0.5 0.5 0.5]);
    
    caxis(datenum([minTime maxTime]));
    S.MarkerFaceAlpha = markerFaceAlpha;
    
    c = colorbar;
    c.TickLabels = datestr(c.Ticks);
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
    
    axis equal;
    
    title(num2str(az_));
    xlim([0 max(rang)]);
    ylim([-40 2]);
end

%%
figure(1);
axis equal;
axis([-81.3 -79.2 -0.85 1.5]);
% 
% for i = 1:8
%     figure(i);
%     fname = '/home/shernandez/regions/pedernales/crossSections/basic_layout';
%     if i > 1
%         fname = strcat('/home/shernandez/regions/pedernales/crossSections/CS_',num2str(i-1));
%     end
%     print('-depsc',fname);
% end
