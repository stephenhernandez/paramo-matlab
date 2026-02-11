function plotFinalDiagnostics2(t,origmag,nMLv,origlon,origlat,...
    newmeanmag,newstdmag,usedPhases,origdepth,minLon,minLat,maxLon,maxLat,printFlag,formatFlag)
newmeanmag = round(newmeanmag*10)/10;
load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
trenchParallelFlag = true;
iscLevels = 0:10:100;
iscLevels = iscLevels(iscLevels>=20);
iscCMap = 'cool';
magThresh = 5.5;
magFact = 3;
pauseTime = 2;
minTime = dn2dt(floor(min(datenum(t)))); %datenum(2015,10,01);
maxTime = dn2dt(ceil(max(datenum(t)))); %datenum(2018,09,15);
winterMap = parula(256);
ilumAz = 25;
wMark = 4;
rotAngle = 20;
torig = datetime(2016,04,16,00,00,00);
tdur = t-torig;
tdur = seconds(tdur);
alpha1 = 0.5;
alpha2 = 0.25;
alpha3 = 0.2;
alpha4 = 0.4;
alpha5 = 0.9;
alpha6 = 0.6;
alpha7 = 0.2;
alpha8 = 0.3;

%%
slope = 0.9136;
mw = slope*(newmeanmag + 0.43) + 0.0085; %newmeanmag;
m0 = m02mw(mw,0); %1.5*mw+9.1;
m0I = t > datetime(2016,04,17) & newmeanmag >= 3;
m0 = m0(m0I);
m0 = cumsum(m0);
m0 = m0*1e-16;

figNumber = 1;
clear ax c;
magI = newmeanmag < 5.5;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 0.8 1]);
ax(1) = subplot(211);
S = scatter(t(magI),newmeanmag(magI),1.5*exp(newmeanmag(magI)),origlat(magI),'filled','linewidth',0.2,'markeredgecolor','k');
hold on;
S = scatter(t(~magI),newmeanmag(~magI),1.5*exp(newmeanmag(~magI)),origlat(~magI),'p','filled','linewidth',0.2,'markeredgecolor','k');
colormap('parula');
c = colorbar;
c.Label.String = 'Latitude';
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
caxis([-2 2]);
ylabel('Magnitude');
ax(1).XTickLabelRotation = rotAngle;
ax(1).YLim = [1 8.5];
yyaxis(ax(1),'right');
%hh = plot([dn2dt(minTime);t],[0;m0],'.-');
hh = plot([datetime(2016,04,17);t(m0I)],[0;m0],'k.-');
ylabel(ax(1),'Cumulative Moment [x$10^{16} Nm$]');
ax(1).XLim = [minTime maxTime];
ax(1).YColor = 'k';
%hh.Visible = 'off';
%ax(1).YLim = [0 4e26];

ax(2) = subplot(212);
edges = dn2dt(floor(datenum(min(t))):ceil(datenum(max(t))));
Ncount = histcounts(t,edges);
stairs(ax(2),edges(1:end-1),Ncount,'.-');
ylabel(ax(2),'Daily Count');
yyaxis(ax(2),'right');
plot(t,1:length(t),'.-');

ax(2).XLim = ([minTime maxTime]);
ax(2).XTickLabelRotation = rotAngle;
ylabel(ax(2),'Cumulative Number');
c(1) = colorbar(ax(1));
c(2) = colorbar(ax(2));
c(2).Visible = 'off';
linkaxes(ax,'x');

S.MarkerFaceAlpha = alpha1; %0.5;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha2; %0.25;

%% figure 2a
disp('contouring coseismic slip');
data = importdata('~/research/now/pedernales/slip_53_s_m.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
disp('done contouring coseismic slip');

cd ~/research/now/pedernales/;
% load('~/igdata/ecuador_coast_dem.mat','lon','lat','dem');
% demdata = dem;
% clear dem;
[lon,lat,demdata] = cutDEM([minLon maxLon minLat maxLat],true);
demdata = demdata';
clear ax c;

figNumber = figNumber+1;
clear ax c;
magI = newmeanmag < magThresh & t < datetime(2016,04,16,23,58,00);

fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
laxes = 0;
laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
cmap = gray(512);

I = dem(lon,lat,demdata,'Contrast',2,'Azimuth',ilumAz,'Interp','noplot','Watermark',wMark,'Elevation',55);
imagesc(ax(laxes),I.x,I.y,I.rgb);
axis xy; pause(pauseTime);
hold(ax(laxes),'on');

if trenchParallelFlag
    strike = 28;
    refEllipse = referenceEllipsoid('wgs84');
    reflat = -2;
    reflon = -82;
    AprimeL1 = -80.15;
    AprimeL2 = 1.9;
    trenchExtension = 500;
    [lat_tmp,lon_tmp] = reckon(reflat,reflon,trenchExtension*1000,strike,refEllipse);
    plot(ax(laxes),[reflon lon_tmp],[reflat lat_tmp],'k','linewidth',2);
    text(ax(laxes),reflon,reflat,'A','FontSize',16)
    text(ax(laxes),AprimeL1,AprimeL2,'$A\prime$','FontSize',16);
end

plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
ylabel(ax(laxes),'Latitude');
xlabel(ax(laxes),'Longitude');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
plot(ax(laxes),lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes + 1;
winterMap = parula(256);
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),datenum(t(magI)),'filled','marker','o','linewidth',1);
%S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),log10(abs(origdepth(magI))),'filled','marker','o','linewidth',1);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
colormap(ax(laxes),winterMap);
c(laxes) = colorbar(ax(laxes));
caxis(ax(laxes),datenum([minTime maxTime]));
%caxis(ax(laxes),[0 2]);
c(laxes).TickLabels = datestr(c(laxes).Ticks);
c(laxes).TickLabelInterpreter = 'latex';
c(laxes).Label.Interpreter = 'latex';

%-----
S.MarkerFaceAlpha = alpha3; %0.2;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha4; %0.4;
%-----

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
contour(ax(laxes),xq,yq,vq,1:7,'linewidth',2);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
colormap(ax(laxes),'hot');
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
magI = newmeanmag > magThresh & newmeanmag < 7.5 & t < datetime(2016,04,16,23,58,00);
S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),datenum(t(magI)),'filled','marker','p','linewidth',1.5,'markeredgecolor','k');
%S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),log10(abs(origdepth(magI))),'filled','marker','p','linewidth',1.5,'markeredgecolor','k');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
colormap(ax(laxes),winterMap);
c(laxes) = colorbar(ax(laxes));
caxis(ax(laxes),datenum([minTime maxTime]));
%caxis(ax(laxes),[0 2]);
c(laxes).Visible = 'off';

%-----
S.MarkerFaceAlpha = alpha5; %0.8;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha6; %0.5;
%-----

magI = newmeanmag > 7.5;
laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,1,ax(laxes));
ax(laxes).Parent = fig(figNumber);
plot(ax(laxes),origlon(magI),origlat(magI),'p','markerfacecolor','w','markeredgecolor','k','markersize',20);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

linkaxes(ax(1:laxes),'xy');
axis(ax(1:laxes),'equal');
axis(ax(1:laxes),[minLon maxLon minLat maxLat]);

%% figure 2b
disp('contouring coseismic slip');
data = importdata('~/research/now/pedernales/slip_53_s_m.dat');
isclon = data(:,1);
isclat = data(:,2);
isc = data(:,3);
[xq,yq] = meshgrid(min(isclon):0.01:max(isclon),min(isclat):0.01:max(isclat));
F = scatteredInterpolant(isclon,isclat,isc,'natural','none');
vq = F(xq,yq);
disp('done contouring coseismic slip');

cd ~/research/now/pedernales/;
% load('~/igdata/ecuador_coast_dem.mat','lon','lat','dem');
% demdata = dem;
% clear dem;
clear ax c;

magI = newmeanmag < magThresh & t >= datetime(2016,04,16,23,58,00);
laxes = 0;
laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
cmap = gray(512);
I = dem(lon,lat,demdata,'Contrast',2,'Azimuth',ilumAz,'Interp','noplot','Watermark',wMark,'Elevation',55);
imagesc(ax(laxes),I.x,I.y,I.rgb);
axis xy; pause(pauseTime);
hold(ax(laxes),'on');

if trenchParallelFlag
    strike = 28;
    refEllipse = referenceEllipsoid('wgs84');
    reflat = -2;
    reflon = -82;
    AprimeL1 = -80.15;
    AprimeL2 = 1.9;
    trenchExtension = 500;
    [lat_tmp,lon_tmp] = reckon(reflat,reflon,trenchExtension*1000,strike,refEllipse);
    plot(ax(laxes),[reflon lon_tmp],[reflat lat_tmp],'k','linewidth',2);
    text(ax(laxes),reflon,reflat,'A','FontSize',16)
    text(ax(laxes),AprimeL1,AprimeL2,'$A\prime$','FontSize',16);
end

plot(ax(laxes),lonEC,latEC,'k','linewidth',2);
ylabel(ax(laxes),'Latitude');
xlabel(ax(laxes),'Longitude');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
plot(ax(laxes),lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes + 1;
winterMap = parula(256);
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),datenum(t(magI)),'filled','marker','o','linewidth',1);
%S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),log10(abs(origdepth(magI))),'filled','marker','o','linewidth',1);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
colormap(ax(laxes),winterMap);
c(laxes) = colorbar(ax(laxes));
caxis(ax(laxes),datenum([minTime maxTime]));
%caxis(ax(laxes),[0 2]);
c(laxes).TickLabels = datestr(c(laxes).Ticks);
c(laxes).TickLabelInterpreter = 'latex';
c(laxes).Label.Interpreter = 'latex';

%-----
S.MarkerFaceAlpha = alpha3; %0.2;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha4; %0.4;
%-----

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
contour(ax(laxes),xq,yq,vq,1:7,'linewidth',2);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
colormap(ax(laxes),'hot');
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
magI = newmeanmag > magThresh & newmeanmag < 7.5 & t >= datetime(2016,04,16,23,58,00);
S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),datenum(t(magI)),'filled','marker','p','linewidth',1.5,'markeredgecolor','k');
%S = scatter(ax(laxes),origlon(magI),origlat(magI),magFact*exp(newmeanmag(magI)),log10(abs(origdepth(magI))),'filled','marker','p','linewidth',1.5,'markeredgecolor','k');
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
colormap(ax(laxes),winterMap);
c(laxes) = colorbar(ax(laxes));
caxis(ax(laxes),datenum([minTime maxTime]));
%caxis(ax(laxes),[0 2]);
c(laxes).Visible = 'off';

%-----
S.MarkerFaceAlpha = alpha5; %0.8;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha6; %0.5;
%-----

magI = newmeanmag > 7.5 & t >= datetime(2016,04,16,23,58,00);
laxes = laxes +1;
ax(laxes) = axes;
subplot(1,2,2,ax(laxes));
ax(laxes).Parent = fig(figNumber);
plot(ax(laxes),origlon(magI),origlat(magI),'p','markerfacecolor','w','markeredgecolor','k','markersize',20);
axis(ax(laxes),'equal');
axis(ax(laxes),[minLon maxLon minLat maxLat]);
ax(laxes).Visible = 'off';
c(laxes) = colorbar(ax(laxes));
c(laxes).Visible = 'off';

linkaxes(ax(1:laxes),'xy');
axis(ax(1:laxes),'equal');
axis(ax(1:laxes),[minLon maxLon minLat maxLat]);

%%
if trenchParallelFlag
    [refDist,refAz] = distance(reflat,reflon,origlat,origlon,refEllipse);
    refDist = refDist/1000;
    origlat = refDist.*cosd(refAz-strike);
    tdiff = abs(t - datetime(2016,04,16,23,58,34.624));
    nlI = tdiff == min(tdiff);
    newRefLat = origlat(nlI); %0.30945;
    origlat = origlat-newRefLat;
    origlon = refDist.*sind(refAz-strike);
end

tI = tdur >= 0 & newmeanmag < magThresh;
figNumber = figNumber+1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'renderer','painters')
ax(1) = axes;
scatter(ax(1),tdur(tI)/86400,origlat(tI),magFact*exp(newmeanmag(tI)),origlon(tI),'.','linewidth',0.1);
c(1) = colorbar(ax(1));
c(1).Visible = 'off';
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(2) = axes;
S = scatter(ax(2),tdur(tI)/86400,origlat(tI),magFact*exp(newmeanmag(tI)),origlon(tI),'filled','linewidth',1);
ax(1).XScale = 'log';
ax(1).XLim = [1 1000];
ax(2).XScale = 'log';
ax(2).XScale = 'log';
ax(2).XLim = [1 1000];
ax(2).XTick = [1 10 100 600];
ax(2).XAxisLocation = 'top';
newLabels = {datestr(1+datenum(2016,04,16),1);datestr(10+datenum(2016,04,16),1);datestr(100+datenum(2016,04,16),1);datestr(600+datenum(2016,04,16),1)};
ax(2).XTickLabel = newLabels;
ax(2).XTickLabelRotation = 15;
grid on
ax(1).YTickLabel = '';
hold(ax(2),'on');

S.MarkerFaceAlpha = alpha7; %0.2;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = alpha8; %0.3;

tI = tdur >= 0 & newmeanmag >= magThresh;
S = scatter(ax(2),tdur(tI)/86400,origlat(tI),magFact*exp(newmeanmag(tI)),origlon(tI),'filled','linewidth',1.5,'marker','p');
S.MarkerFaceAlpha = alpha7; %0.2;
S.MarkerEdgeColor = 'k';
S.MarkerEdgeAlpha = 0.9;

c(2) = colorbar(ax(2));
if trenchParallelFlag
    caxis(ax(2),[0 100]);
    c(2).Label.String = 'Distance from Trench [km.]';
else
    caxis(ax(2),[-82 -79]);
    c(2).Label.String = 'Longitude';
end
c(2).TickLabelInterpreter = 'latex';
c(2).Label.Interpreter = 'latex';

if trenchParallelFlag
    ax(1).YLim = [-300 150];
    ax(2).YLim = [-300 150];
    ax(2).YLabel.String = 'Distance along Trench [km.]';
else
    ax(1).YLim = [minLat maxLat];
    ax(2).YLim = [minLat maxLat];
    ax(2).YLabel.String = 'Latitude';
end
ax(1).XLim = [1 1000];
ax(2).XLim = [1 1000];
ax(1).XLabel.String = ['Time Since ',datestr(torig),' [days]'];
if trenchParallelFlag
    [refDist,refAz] = distance(reflat,reflon,0.47,-80.12,refEllipse);
    refDist = refDist*1e-3;
    maxTr = refDist.*cosd(refAz-strike);
    [refDist,refAz] = distance(reflat,reflon,-0.58,-80.22,refEllipse);
    refDist = refDist*1e-3;
    minTr = refDist.*cosd(refAz-strike);
    plot(ax(2),[1e0 1e0],[minTr maxTr]-newRefLat,'color',[0.5 0.5 0.5],'linewidth',5);
    plot(ax(2),1e0,0,'p','markerfacecolor','w','markeredgecolor','k','markersize',20);
else
    plot(ax(2),[1e0 1e0],[-0.58 0.47],'color',[0.5 0.5 0.5],'linewidth',5);
end
%colormap jet;

%%
if formatFlag
    formatStr = '-djpeg';
else
    formatStr = '-dsvg';
end

if printFlag
    for i = figNumber:-1:1
        fname = ['~/regions/pedernales/paper/figure_',num2str(i)];
        figure(fig(i));
        pause(pauseTime);
        print(formatStr,fname,'-painters');
        %export_fig(fname,'-eps','-transparent');
    end
end
end