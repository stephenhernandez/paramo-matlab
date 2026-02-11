clear; close all;
cd ~/research/now/sierra_negra/NLL_TEPP_20190904/BULLETINS;

files = dir('dias*');

E = populateSeisCompStructure(length(files));
for i = 1:length(files)
    E(i) = readSCBulletin(files(i).name,true);
end

t = pull(E,'t');
eqmag = pull(E,'mag');
eqlat = pull(E,'lat');
eqlon = pull(E,'lon');
eqdepth = pull(E,'depth');

[t,sI] = sort(t);
eqmag = eqmag(sI);
eqlat = eqlat(sI);
eqlon = eqlon(sI);
eqdepth = eqdepth(sI);

refEllipse = referenceEllipsoid('wgs84');
referenceLat = -0.815;
referenceLon = -91.13;

%%
clearvars -except t eq* ref* ss* %reference* min* max* dem*
defFontSize = 16;
set(groot,'defaultAxesFontSize',defFontSize);

GV04_lat = -0.81;
GV04_lon = -91.14;

dataGPS = load('~/research/now/sierra_negra/gps/geodesy/gala_daily_tseries/GV04.tseries');
tGPS = datetime(dataGPS(:,12),dataGPS(:,13),dataGPS(:,14),dataGPS(:,15),dataGPS(:,16),dataGPS(:,17));
tseries = 1e3*dataGPS(:,2:4);
zGPS = tseries(:,3);

minOpacity = 0.25;
minLat = -0.87;
maxLat = -0.77;
minLon = -91.185;
maxLon = -91.085;
maxDepth = 6;
minT = datetime(2018,04,20);
maxT = datetime(2018,08,20);
minMag = 1;

[demLon,demLat,demData] = cutDEM([minLon,maxLon,minLat,maxLat],true);

fI = eqlon >= minLon & eqlon <= maxLon & eqlat >= minLat ...
    & eqlat <= maxLat & eqdepth <= maxDepth ...
    & eqmag >= minMag & t >= minT & t <= maxT;

t2 = t(fI);
zI = tGPS >= dateshift(min(t2),'start','day') & tGPS <= dateshift(max(t2),'end','day');
zStair = interp1(datenum(tGPS(zI)),zGPS(zI),datenum(t2),'pchip');
zStair = zStair - zStair(1);

[dist,azs] = distance(referenceLat,referenceLon,eqlat(fI),eqlon(fI),refEllipse);
azs(azs > 60) = azs(azs > 60) - 360;

close all;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig.Visible = 'on';

ax(1) = subplot(4,3,[1 4]);
contour(ax(1),demLon,demLat,demData,(500:50:1150)','Color',[0.15 0.15 0.15]);
hold(ax(1),'on');
h1 = plot(ax(1),referenceLon,referenceLat,'r.','MarkerSize',20);
h2 = plot(ax(1),GV04_lon,GV04_lat,'p','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','g');
ax(1).XAxisLocation = 'top';

axis equal; grid on;
xlim([minLon maxLon]);
ylim([minLat maxLat]);
xlab1 = xlabel('Lon.','FontSize',defFontSize);
ylab1 = ylabel('Lat.','FontSize',defFontSize);
ax(1).Color = [0.85 0.85 0.85];

ax(2) = subplot(4,3,[2 3 5 6]);
ax(2).Color = [0.85 0.85 0.85];
c2 = colorbar(ax(2)); c2.Visible = 'off';
grid(ax(2),'on');
hold(ax(2),'on');
ax(2).XAxisLocation = 'top';
ylab2 = ylabel('Cum. Uplift [mm.] Since 21-Apr-2018','FontSize',defFontSize);
%ylab2 = ylabel('Magnitude [Ml]','FontSize',defFontSize);

ax(3) = subplot(4,3,[7 8 9 10 11 12]);
ax(3).Color = [0.85 0.85 0.85];
grid(ax(3),'on');
hold(ax(3),'on');

xlab3 = xlabel('Azimuth','FontSize',defFontSize);
ylab3 = ylabel('Depth [km.]','FontSize',defFontSize);

text(-225,1.65,'Southern TDF','FontSize',defFontSize);
text(-105,1.65,'Western TDF','FontSize',defFontSize);
text(-10,1.65,'Northern TDF','FontSize',defFontSize);

pause(2);
opacityDecayRate = -0.125;
fI = find(fI);
lg = length(fI);

%if ~exist('ss1','var')
ss1 = gobjects(lg,1);
ss2 = ss1;
ss3 = ss1;
%end
set(gcf, 'InvertHardCopy', 'off');

for i = 1:lg
    disp(i)
    ss1(i) = scatter(ax(1),eqlon(fI(i)),eqlat(fI(i)),5*exp(eqmag(fI(i))),...
        datenum(t(fI(i))),'p','filled');
    ss1(i).MarkerFaceAlpha = 1;
    ss1(i).MarkerEdgeColor = 'k';
    ss1(i).MarkerEdgeAlpha = 1;
    caxis(ax(1),datenum([minT maxT]));
    
    %     ss2(i) = scatter(ax(2),t(fI(i)),eqmag(fI(i)),5*exp(eqmag(fI(i))),...
    %         datenum(t(fI(i))),'p','filled');
    %     ss2(i).MarkerFaceAlpha = 1;
    %     ss2(i).MarkerEdgeColor = 'k';
    %     ss2(i).MarkerEdgeAlpha = 1;
    %     ylim(ax(2),[minMag 5]);
    
    ss2(i) = scatter(ax(2),t(fI(i)),zStair(i),50,datenum(t(fI(i))),'o','filled');
    ss2(i).MarkerFaceAlpha = 1;
    ss2(i).MarkerEdgeColor = 'k';
    ss2(i).MarkerEdgeAlpha = 1;
    caxis(ax(2),datenum([minT maxT]));
    xlim(ax(2),[minT maxT]);
    
    ss3(i) = scatter(ax(3),azs(i),-eqdepth(fI(i)),5*exp(eqmag(fI(i))),...
        datenum(t(fI(i))),'p','filled');
    ss3(i).MarkerFaceAlpha = 1;
    ss3(i).MarkerEdgeColor = 'k';
    ss3(i).MarkerEdgeAlpha = 1;
    
    xlim(ax(3),[-300 60]);
    ylim(ax(3),[-maxDepth 2]);
    caxis(ax(3),datenum([minT maxT]));
    c3 = colorbar;
    c3.TickLabels = datestr(c3.Ticks);
    
    %%
    if i > 1
        ss1(i-1).Marker = 'o';
        %ss2(i-1).Marker = 'o';
        ss3(i-1).Marker = 'o';
        
        for j = 1:i-1 %just change the opacity
            age_ = days(t(fI(i)) - t(fI(j)));
            opacity_ = exp(opacityDecayRate*age_);
            ss1(j).MarkerFaceAlpha = max([opacity_ minOpacity]);
            ss2(j).MarkerFaceAlpha = max([opacity_ minOpacity]);
            ss3(j).MarkerFaceAlpha = max([opacity_ minOpacity]);
            
            ss1(j).MarkerEdgeAlpha = max([opacity_ minOpacity]);
            ss2(j).MarkerEdgeAlpha = max([opacity_ minOpacity]);
            ss3(j).MarkerEdgeAlpha = max([opacity_ minOpacity]);
            
            ss1(j).SizeData = 4*exp(eqmag(fI(j)));
            %ss2(j).SizeData = 4*exp(eqmag(fI(j)));
            ss3(j).SizeData = 4*exp(eqmag(fI(j)));
        end
    end
    
    suptitle(['N: ',num2str(i),', t: ',datestr(t(fI(i)),'yyyy-mm-dd HH:MM:SS'),', Ml: ',num2str(eqmag(fI(i)))]);
    pause(2);
    if i < 10
        counterStr = ['000',num2str(i)];
    elseif i < 100
        counterStr = ['00',num2str(i)];
    elseif i < 1000
        counterStr = ['0',num2str(i)];
    else
        counterStr = num2str(i);
    end
    
    %%
    delete(h1);
    delete(h2);
    h1 = plot(ax(1),referenceLon,referenceLat,'r.','MarkerSize',20);
    h2 = plot(ax(1),GV04_lon,GV04_lat,'p','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','g');
    
    %%
    fname = ['~/Desktop/frames/frame_',counterStr];
    print(fname,'-djpeg');    
end
