%function plotEcuadorSeismicity(tRef,data,boundaryBox,urbanFlag,printFlag)
function plotEcuadorSeismicity(t,data,boundaryBox,urbanFlag,printFlag)
fig = loadBasemap(boundaryBox,'pink',urbanFlag,true);
cbar1 = findobj(fig.Children,'Type','Colorbar');
ax = findobj(fig.Children,'Type','Axes');

%%
%Zlimits = [-5476; 5957];
%fig.Colormap = demcmap(Zlimits,256,winter(64),flipud(pink(64)));
if urbanFlag
    geoshow(ax,'~/igdata/ZonaUrbana/ZonaUrbana.shp');
end

%%
vFlag = 0;
magFact = 8;
MarkerFaceAlpha = 0.4;
cbarLocation = cbar1.Location;
MarkerEdgeColor = 'k';
ScatterLineWidth = 0.4;
maxDepth = 200;

if vFlag
    fig.Visible = 'on';
end

%%
[t,sI] = sort(t);
eqlat = data(sI,1);
eqlon = data(sI,2);
eqdepth = data(sI,3);
eqmag = data(sI,4);

%%
nDays = 7;
stride = days(nDays);
tStart1 = dateshift(min(t),'start','day');
tEnd1 = dateshift(max(t),'end','day');
tStarts = (tStart1:tEnd1-stride+1)';

%%
for i = 1:length(tStarts)
    disp(i);
%     if mod(i,2)
%         fig.Visible = 'on';
%     else
%         fig.Visible = 'off';
%     end
    tStart_ = tStarts(i);
    tEnd_ = dateshift(tStart_ + stride - 1,'end','day');
    fName = strcat('~/products/summaries/EcuadorWeeklySeismicity_',datestr(tStart_,'yyyymmdd'));

    laxes = 2;
    tI = t >= tStart_ & t <= tEnd_;
    ntotal = sum(tI);
    if ~ntotal
        continue;
    end

    ax(laxes) = axes;
    maxmag = max(eqmag(tI));
    minmag = min(eqmag(tI));
    %medmag = median(eqmag(tI));

    %%
    S = scatter(ax(laxes),eqlon(tI),eqlat(tI),magFact*exp(eqmag(tI)),eqdepth(tI),'o','filled');
    zoom on;
    S.MarkerFaceAlpha = MarkerFaceAlpha;
    S.MarkerEdgeColor = MarkerEdgeColor;
    S.LineWidth = ScatterLineWidth;
    c = colorbar(ax(laxes),cbarLocation);
    axis(ax(1:laxes),'equal');

    cmap2 = (turbo(256)); %'parula';
    colormap(ax(laxes),cmap2);
    c.Label.String = 'depth [km.]';
    c.Label.Interpreter = 'latex';
    caxis([1 maxDepth]);
    set(gca,'ColorScale','log');
    h = title(ax(1),{strcat(datestr(tStart_,'yyyy/mm/dd'),' - ',datestr(tEnd_,'yyyy/mm/dd')); ...
        ['$N_{total}$=',num2str(ntotal),', Max. Mag.=',num2str(maxmag),', Min. Mag.=',num2str(minmag)]});
    xlabel(ax(1),'Longitude');
    ylabel(ax(1),'Latitude');

    %%
    ax(laxes).Visible = 'off';
    %minLon = boundaryBox;
    axis(ax(1:laxes),boundaryBox); %[minLon maxLon minLat maxLat]);
    linkaxes(ax(1:end),'xy');

    %%
    if printFlag
        print('-dpng',fName);
        delete(ax(2));
        delete(h);
    end
end
