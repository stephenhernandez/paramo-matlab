%function
%thermal_animation
clear; close all; clc;

%cd ~/research/now/reventador/ash_shape_analysis/
%fname = "thermal_data.xlsx";
%fname = "suomi_virs.xls";

cd ~/igdata
fname = "frp_reventador_silvia.xlsx";

if strcmp(fname,"thermal_data.xlsx")
    [num,txt] = xlsread(fname);
    timestr = char(num2str(num(:,4)));
    min_ = str2double(string(timestr(:,end-1:end)));
    hour_ = str2double(string(timestr(:,1:end-2)));

    acqtime = dn2dt(num(:,3) + 693960);
    acqtime = acqtime+hours(hour_) + minutes(min_);
    t = acqtime;
    tlat = num(:,1);
    tlon = num(:,2);
    temp = num(:,7);
    frp = num(:,8);
    daynight = string(txt(:,9));
    lT = length(t);
    tStart = min(t);
    tEnd = max(t);
else
    [num,txt] = xlsread(fname);
    timestr = char(num2str(num(:,7)));
    min_ = str2double(string(timestr(:,end-1:end)));
    hour_ = str2double(string(timestr(:,1:end-2)));

    %latitude	longitude	brightness	scan	track	acq_date	acq_time	Daily Alerts	satellite	instrument	confidence	version	bright_t31	frp	daynight	type
    acqtime = dn2dt(num(:,6) + 693960);
    acqtime = acqtime+hours(hour_) + minutes(min_);
    t = acqtime;
    tlat = num(:,1);
    tlon = num(:,2);
    temp = num(:,13);
    frp = num(:,14);
    daynight = string(txt(:,14));
    lT = length(t);
    tStart = min(t);
    tEnd = max(t);
end

%% some parameteres for hillshade option
fstring = '-djpeg';
boundaryBox = getRegionSpatialDimensions('reventador');
ilumAz = -70;
wMark = 2;
elevAngle = 60;
fontSize = 18;
magFact = 10;

figNumber = 100;
[lon,lat,demData] = cutDEM(boundaryBox,true);
fig(figNumber) = figure('units','normalized','outerposition',[0 0 1 1]);
fig(figNumber).Visible = 'off';

ax = subplot(5,1,[1 2 3 4]);
I = dem(lon,lat,demData','Contrast',1,'Azimuth',-70,'Interp','noplot','Watermark',2,'Elevation',60);

laxes = 1;
imagesc(ax(laxes),I.x,I.y,I.rgb);
hold(ax(laxes),'on');
axis(ax(laxes),'xy');
axis(ax(laxes),'equal')
axis(ax(laxes),'tight')
ylabel(ax(laxes),'Latitud','FontSize',fontSize);
xlabel(ax(laxes),'Longitud','FontSize',fontSize);

%ax(lT+2) = subplot(5,1,5);

%%
pointColormap = 'turbo';
allFrames = true;
markerFaceAlpha = 0.75;
LT = false;
timeColorCode = true;
dI = tlon >= boundaryBox(1) &tlon <= boundaryBox(2) & tlat >= boundaryBox(3) & tlat <= boundaryBox(4);
dI = find(dI);
lDI = length(dI);

if allFrames
    lDI = 1:lDI;
end
for i = lDI
    disp(i);
    daynight_ = daynight(i);
    frp_ = frp(i);
    lat_ = tlat(i);
    lon_ = tlon(i);
    t_ = t(i);
    t31_ = temp(i);

    %%
    ax(laxes+1) = axes;
    subplot(5,1,[1 2 3 4],ax(laxes+1));


    S = scatter(ax(laxes+1),tlon(dI(1:i)),tlat(dI(1:i)),magFact*(frp(dI(1:i))),datenum(t(dI(1:i))),'filled');
    clim(ax(laxes+1),datenum([tStart tEnd]));
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.TickLabels = datestr(c.Ticks);
    colormap(ax(laxes+1),pointColormap);

    hold(ax(laxes+1),'on');
    axis(ax(laxes+1),'equal');

    %%
    if LT
        titleStr = [datestr(t(dI(i))),' (Tiempo Local), N=',num2str(i)];
    else
        titleStr = [datestr(t(dI(i))),' (UTC), N=',num2str(i),', FRP=',num2str(frp_),', $T_{31}$=',num2str(t31_)];
    end
    h = title(ax(laxes),titleStr,'FontSize',fontSize);
    S.MarkerFaceAlpha = markerFaceAlpha;
    S.MarkerEdgeColor = 'k';
    S.MarkerEdgeAlpha = markerFaceAlpha;

    %%
    linkaxes(ax(laxes:laxes+1));
    ax(laxes+1).Visible = 'off';
    axis(ax(laxes+1),boundaryBox); %[minLon maxLon minLat maxLat]);

    if i < 10
        fname = strcat('~/research/now/reventador/ash_shape_analysis/animations/frame_000',num2str(i));
    elseif i < 100
        fname = strcat('~/research/now/reventador/ash_shape_analysis/animations/frame_00',num2str(i));
    elseif i < 1000
        fname = strcat('~/research/now/reventador/ash_shape_analysis/animations/frame_0',num2str(i));
    else
        fname = strcat('~/research/now/reventador/ash_shape_analysis/animations/frame_',num2str(i));
    end


    %% bottom panel
    ax(laxes+2) = axes;
    ax(laxes+2).Parent = fig(figNumber);
    subplot(5,1,5,ax(laxes+2));

    S = scatter(ax(laxes+2),t(dI(1:i)),temp(dI(1:i)),magFact*(frp(dI(1:i))),datenum(t(dI(1:i))),'filled','linewidth',0.2,'markeredgecolor','k');
    caxis(ax(laxes+2),datenum([tStart tEnd]));
    colormap(ax(laxes+2),pointColormap);

    c = colorbar;
    xlim(ax(laxes+2),[tStart tEnd]);
    ylabel(ax(laxes+2),'Temperature [K]','FontSize',fontSize);
    c.Visible = 'off';
    S.MarkerFaceAlpha = markerFaceAlpha;

    %%
    if allFrames
        disp(fname);
        print(fstring,fname);
        if i < max(lDI)
            delete(h);
            delete(ax(laxes+1));
            delete(ax(laxes+2));

        end
    else
        disp(fname);
    end
end
