%%
clear;
close all;
%cd ~/research/now/sn_eruption/geodesy/gala_daily_tseries/
cd ~/research/now/sierra_negra/gps/

files = dir('*.tseries');
n = 0;
for i = 1:length(files)
    disp(i);
    n = n+1;
    fname = files(i).name;
    stnm = fname(1:4);
    data = load(fname);
    t = datetime(data(:,12),data(:,13),data(:,14),data(:,15),data(:,16),data(:,17));
    tseries = 1e3*data(:,2:4);
    figure('units','normalized','outerposition',[0 0 1 0.75]);
    h = plot(t,tseries,'o');
    title(stnm);
    lgd = legend(h,'E','N','Z','location','NorthWest');
    lgd.Orientation = 'horizontal';
end

%%
clear;
close all;

printFlag = false;
plotFlag = true;

regionName = 'sierra_negra';
demName = 'sierraNegra';
minLat= -0.88;
maxLat= -0.76;
minLon= -91.19;
maxLon= -91.08;
animFlag = false;
fNamePrefix = '~/regions/sierra_negra/animation/GPS_frame_';

%cd ~/regions/sierra_negra/Sierra_Negra_March2018/
%cd ~/research/now/sn_eruption/geodesy/gala_daily_tseries/
cd ~/research/now/sierra_negra/gps/
files = dir('*.tseries');
I = [2 4:11];
minT = datetime(2017,01,01);
maxT = datetime(2018,03,10);
tMaster = hours(12) + minT:days(1):maxT;
newTseries = NaN(length(tMaster),3);
E = NaN(length(tMaster),length(I));
N = E;
Z = E;

Ndegree = 3;
n = 0;
for i = I
    disp(i);
    n = n+1;
    fname = files(i).name;
    data = load(fname);
    stnm = fname(1:4);
    t = datetime(data(:,12),data(:,13),data(:,14),data(:,15),data(:,16),data(:,17));
    tI = t >= minT;
    tseries = 1e3*data(tI,2:4);
    tnew= t(tI);
    clear residuals;
    residuals = tseries;
    
    for j = 1:3
        p = polyfit(datenum(tnew)-min(datenum(tnew)),tseries(:,j),Ndegree);
        newTseries(:,j) = polyval(p,datenum(tMaster)-min(datenum(tMaster)));
        residuals(:,j) = tseries(:,j) - polyval(p,datenum(tnew)-min(datenum(tnew)));
        tseries(:,j) = tseries(:,j) - newTseries(1,j);
        newTseries(:,j) = newTseries(:,j) - newTseries(1,j);
    end
    if plotFlag
        figure('units','normalized','outerposition',[0 0 1 0.75]);
        subplot(3,3,[1 4 7]);
        h = plot(tnew,tseries,'o');
        hold on;
        plot(tMaster,newTseries,'k-','linewidth',3);
        ylabel('[mm.]');
        title(['Degree ',num2str(Ndegree),' polynomial fit']);
        lgd = legend(h,'E','N','Z','location','NorthWest');
        lgd.Orientation = 'horizontal';
        
        subplot(3,3,2);
        plot(tnew,residuals(:,1),'o');
        hold on;
        plot([min(tnew) max(tnew)],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2);
        title('E residuals');
        ylabel('[mm.]');
        subplot(3,3,5);
        plot(tnew,residuals(:,2),'o');
        hold on;
        plot([min(tnew) max(tnew)],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2);
        title('N residuals');
        ylabel('[mm.]');
        subplot(3,3,8);
        plot(tnew,residuals(:,3),'o');
        hold on;
        plot([min(tnew) max(tnew)],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2);
        ylabel('[mm.]');
        title('Z residuals');
        
        subplot(3,3,3);
        histogram(residuals(:,1),20);
        title(['E res. dist., rms = ',num2str(rms(residuals(:,1)))]);
        xlabel('[mm.]');
        subplot(3,3,6);
        histogram(residuals(:,2),20);
        title(['N res. dist., rms = ',num2str(rms(residuals(:,2)))]);
        xlabel('[mm.]');
        subplot(3,3,9);
        histogram(residuals(:,3),20);
        title(['Z res. dist., rms = ',num2str(rms(residuals(:,3)))]);
        xlabel('[mm.]');
        suptitle(stnm);
        
        if printFlag
            fname = strcat(stnm,'_polyfit');
            print('-djpeg',fname);
        end
    end
    
    [lia,locb] = ismember(tnew,tMaster);
    newTseries(locb,:) = tseries;
    E(:,n) = newTseries(:,1);
    N(:,n) = newTseries(:,2);
    Z(:,n) = newTseries(:,3);
end

maxT = dateshift(maxT,'start','day');
stnms = ["GV01","GV03","GV04","GV05","GV06","GV07","GV08","GV09","GV10"];
[stla,stlo,stel] = metaDataFromStationList(stnms);

if plotFlag
    figure('units','normalized','outerposition',[0 0 1 0.8]);
    quiver(stlo,stla,E(end,:)',N(end,:)');
    axis equal;
    load('~/igdata/ec_boundaries','latEC','lonEC');
    hold on;
    plot(lonEC,latEC,'k','linewidth',2);
end

%%
close all;
figNumber = 100;
laxes = 1;
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
ax(laxes) = axes;
ax(laxes).Parent = fig(figNumber);

load('~/igdata/soam_noec','lat_noec','lon_noec');
load('~/igdata/ec_boundaries','latEC','lonEC');
%load(strcat('~/regions/',regionName,'/dem/',demName),'lon','lat','demData');

refEllipse = referenceEllipsoid('wgs84');

% cut the dem data to custom bounds
[lon,lat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);
%[lon,lat,demData] = cut_regional_dem(lon,lat,demData,minLat,maxLat,minLon,maxLon,true);

minElev = min(min(demData));
minElev = round(minElev/100)*100;
minElev = max([minElev 0]);

maxElev = max(max(demData));
maxElev = round(maxElev/100)*100;

levInc = maxElev - minElev;
levInc = levInc/1000;
levInc = 10*floor(levInc);
levInc = max([levInc 10]);

zlevs = minElev:levInc:maxElev;
contour(ax(laxes),lon,lat,demData',zlevs);

disp(strcat('Min. Elev.: ',num2str(minElev)));
disp(strcat('Max. Elev.: ',num2str(maxElev)));
disp(strcat('Contour Increment: ',num2str(levInc)));

hold(ax(laxes),'on');
axis(ax(laxes),'xy');
stationOutlineColor = 'k';

ylabel(ax(laxes),'Latitud');
xlabel(ax(laxes),'Longitud');
axis(ax(laxes),'equal');
%geoshow('~/igdata/fallas2008completas.shp','Color','r','linewidth',2);
%S = shaperead('~/igdata/volcanes_2.shp');
%plot(ax(laxes),[S.X],[S.Y],'k--','linewidth',1);
text(ax(laxes),stlo+0.001,stla+0.001,stnms,'fontsize',15,'FontName','Castellar');
colorbar;
colormap copper

plot(ax(laxes),stlo,stla,'v','markeredgecolor',stationOutlineColor,'markerfacecolor','w','markersize',15,'linewidth',1);
plot(ax(laxes),lonEC,latEC,'k','linewidth',3);
axis(ax(laxes),[minLon maxLon minLat maxLat]);

ampFact1 = 5.e-5;
ampFact2 = 0.5*ampFact1;

quiver(ax(laxes),-91.188,-0.762,ampFact1*500,0,0,'k','linewidth',3); text(ax(laxes),-91.16,-0.762,'500 mm. horiz.','fontsize',15,'FontName','Castellar');
quiver(ax(laxes),-91.188,-0.764,ampFact2*500,0,0,'r','linewidth',3); text(ax(laxes),-91.174,-0.764,'500 mm. vert.','fontsize',15,'FontName','Castellar');

for i = 1:size(N,1)
    figure(1);
    hq = quiver(ax(laxes),stlo,stla,ampFact1*E(i,:)',ampFact1*N(i,:)',0,'k','linewidth',3);
    hold on;
    hv = quiver(ax(laxes),stlo,stla,zeros(length(I),1),ampFact2*Z(i,:)',0,'r','linewidth',3);
    ht = title(datestr(dateshift(tMaster(i),'start','day')));
    pause(0.1);
    if animFlag
        if i < 10
            frameName = strcat(fNamePrefix,'000',num2str(i));
        elseif i < 100
            frameName = strcat(fNamePrefix,'00',num2str(i));
        elseif i < 1000
            frameName = strcat(fNamePrefix,'0',num2str(i));
        else
            frameName = strcat(fNamePrefix,num2str(i));
        end
        print('-djpeg',frameName);
    end
    delete(hq);
    delete(hv);
    delete(ht);
end
