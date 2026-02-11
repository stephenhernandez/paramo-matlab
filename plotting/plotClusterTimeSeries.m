function [multipletI,lf,singletonI,lSingletons,ax] = plotClusterTimeSeries(families,...
    tabs,maxAmpRMS,MAGFLAG,LOGFLAG)
if nargin < 4
    MAGFLAG = false;
end

if nargin < 5
    LOGFLAG = false;
end

%%
LABELTHRESHOLD = 10;
lGroups = length(families);
minT = NaN(lGroups,1);
t1 = tabs(1);
if isduration(t1) || isdatetime(t1)
    minT = NaT(lGroups,1);
    DTFLAG = true;
end
maxT = minT;
lf = NaN(lGroups,1);

%%
for i = 1:lGroups
    fam_ = families{i};
    lf(i) = length(fam_);
    t_ = tabs(fam_);
    minT(i) = min(t_);
    maxT(i) = max(t_);
end
nTot = sum(lf);

%%
multipletI = lf > 1;
lMultiplets = sum(multipletI);
singletonI = ~multipletI;
lSingletons = sum(singletonI);
minT = minT(multipletI);
maxT = maxT(multipletI);
lf = lf(multipletI);
multipletI = families(multipletI);

%%
[minT,sI] = sort(minT);
multipletI = multipletI(sI);
maxT = maxT(sI);
lf = lf(sI);
n_multiplet_members = sum(lf);

%%
figure("units","normalized","outerposition",[0.05 0.05 0.6 0.9]);
nSubplots = 5;
tile_spacing = "compact";
tiledlayout(nSubplots,1,"Padding","compact","TileSpacing",tile_spacing);
ax(1) = nexttile(1,[nSubplots-1,1]);
hold on; grid on; zoom on;

magFact = 50;
if MAGFLAG
    symbolSize = magFact*(maxAmpRMS-min(maxAmpRMS)+0.1);
else
    symbolSize = magFact*exp(log10(maxAmpRMS));
end
lw_thin = 0.1;
LW2 = 1;
MarkerFaceAlpha = 0.8;
MarkerEdgeColor = "k";
MarkerEdgeAlpha = 0.3;

%%
GLOBALMINT = min(tabs);
if DTFLAG
    GLOBALMINT = dateshift(GLOBALMINT,"start","day");
end

mI = cat(1,multipletI{:});
for i = 1:lMultiplets
    cI = multipletI{i};
    lMembers = lf(i);
    t_ = tabs(cI);
    amp_ = maxAmpRMS(cI);
    symbol_size = symbolSize(cI);
    if LOGFLAG
        if ~DTFLAG
            fprintf("cannot apply log time scale to non-datetime type...\n");
            return;
        end
        pp = plot(ax(1),seconds([minT(i) maxT(i)]-GLOBALMINT),[i i],'-','color',[0.5 0.5 0.5],"linewidth",lw_thin);
    else
        pp = plot(ax(1),[minT(i) maxT(i)],[i i],'-','color',[0.5 0.5 0.5],"linewidth",lw_thin);
    end
    pp.Color(4) = MarkerEdgeAlpha;

    y_ = i*ones(lf(i),1);
    if LOGFLAG
        ss = scatter(ax(1),seconds(t_-GLOBALMINT),y_,symbol_size,amp_,'filled');
    else
        ss = scatter(ax(1),t_,y_,symbol_size,amp_,'filled');
    end

    if lMembers >= LABELTHRESHOLD
        if LOGFLAG
            text(ax(1),seconds(maxT(i)-GLOBALMINT),i,sprintf("%d",lMembers),...
                "FontSize",14,"FontWeight","light");
        else
            text(ax(1),maxT(i),i,sprintf("%d",lMembers),...
                "FontSize",14,"FontWeight","light");
        end
    end

    ss.MarkerFaceAlpha = MarkerFaceAlpha;
    ss.MarkerEdgeColor = MarkerEdgeColor;
    ss.MarkerEdgeAlpha = MarkerEdgeAlpha;
    ss.LineWidth = LW2;
end

colorbar;
ax(1).Box = "on";
ax(1).XTickLabel = [];
titleStr = sprintf("total: %d; multiplets: %d; " + ...
    "multiplet members: %d; number singletons: %d",...
    nTot,lMultiplets,n_multiplet_members,lSingletons);
title(titleStr);
ylabel("sequential family number");

%%
ylim_ = ax(1).YLim;
ymax = ylim_(2);
if (ymax-lMultiplets)/ymax >= 0.1
    ymax = 5*ceil(lMultiplets/5);
    ax(1).YLim = [ylim_(1) ymax];
end

%%
ax(2) = nexttile();
hold on; grid on; zoom on;
if LOGFLAG
    plot(ax(2),seconds(tabs(mI)-GLOBALMINT),maxAmpRMS(mI),".");
else
    plot(ax(2),tabs(mI),maxAmpRMS(mI),".");
end
ax(2).Box = "on";
linkaxes(ax,"x");

%%
if ~lSingletons
    fprintf("no singletons... are you sure?\n");
    return;
end

singletonI = cat(1,families{singletonI});
t_ = tabs(singletonI);
amp_ = maxAmpRMS(singletonI);
symbol_size = symbolSize(singletonI);
if LOGFLAG
    pp = plot(ax(2),seconds(t_-GLOBALMINT),amp_,".","Color",[0.5 0.5 0.5],"MarkerSize",12);
    pp.Color(4) = MarkerEdgeAlpha;
    ax(2).XScale = "log";
    pp = plot(ax(1),seconds([min(t_) max(t_)]-GLOBALMINT),[0 0],"-","color",[0.5 0.5 0.5],"linewidth",lw_thin);
    pp.Color(4) = MarkerEdgeAlpha;
else
    pp = plot(ax(2),t_,amp_,".","Color",[0.5 0.5 0.5],"MarkerSize",12);
    pp.Color(4) = MarkerEdgeAlpha;
    pp = plot(ax(1),[min(t_) max(t_)],[0 0],"-","color",[0.5 0.5 0.5],"linewidth",lw_thin);
    pp.Color(4) = MarkerEdgeAlpha;
end

y_ = zeros(lSingletons,1);
if LOGFLAG
    ss = scatter(ax(1),seconds(t_ - GLOBALMINT),y_,symbol_size,amp_,"filled");
    ax(1).XScale = "log";
else
    ss = scatter(ax(1),t_,y_,symbol_size,amp_,"filled");
end

ss.MarkerFaceAlpha = MarkerFaceAlpha;
ss.MarkerEdgeColor = MarkerEdgeColor;
ss.MarkerEdgeAlpha = MarkerEdgeAlpha;
ss.LineWidth = LW2;
legend(ax(2),"multiplet","singletons",...
    "Orientation","horizontal","Location","best");