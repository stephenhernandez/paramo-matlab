function [nonsingletonsI,lf,singletonI,lf0] = plotTemplateTimeSeries(templateNumber,tMain,ampMain)
uniqTemplateNumber = unique(templateNumber);
lFamilies = length(uniqTemplateNumber);
lf = groupcounts(templateNumber);

minT = NaT(lFamilies,1);
maxT = minT;
families = cell(lFamilies,1);
for i = 1:lFamilies
    fam_ = find(templateNumber == uniqTemplateNumber(i));
    families{i} = fam_;
    t_ = tMain(fam_);
    minT(i) = min(t_);
    maxT(i) = max(t_);
end

%%
singletonI = find(lf == 1);
lSingletons = length(singletonI);
if lSingletons
    singletonI = cat(1,families{singletonI});
    lf0 = length(singletonI);
    minT0 = min(tMain(singletonI));
    maxT0 = max(tMain(singletonI));
else
    singletonI = [];
    lf0 = [];
end

nonsingletonsI = find(lf ~= 1);
minT = minT(nonsingletonsI);
maxT = maxT(nonsingletonsI);
lf = lf(nonsingletonsI);
uniqTemplateNumber = uniqTemplateNumber(nonsingletonsI);
nonsingletonsI = families(nonsingletonsI);

%%
[minT,sI] = sort(minT);
nonsingletonsI = nonsingletonsI(sI);
maxT = maxT(sI);
lf = lf(sI);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
lFamilies = length(nonsingletonsI);
if lSingletons
    i = 0;
    t_ = tMain(singletonI);
    amp_ = ampMain(singletonI);

    pp = plot([minT0 maxT0],[i i],'-','color',[0.5 0.5 0.5],'linewidth',2);
    pp.Color(4) = 0.9;

    y_ = i*ones(lf0,1);
    ss = scatter(t_,y_,5*exp(log10(amp_)),log10(amp_),'filled');
    ss.MarkerFaceAlpha = 0.4;
    ss.MarkerEdgeColor = 'k';
    ss.MarkerEdgeAlpha = 0.4;
    ss.LineWidth = 0.1;
end

for i = 1:lFamilies
    cI = nonsingletonsI{i};
    t_ = tMain(cI);
    amp_ = ampMain(cI);

    pp = plot([minT(i) maxT(i)],[uniqTemplateNumber(i) uniqTemplateNumber(i)],...
        '-','color',[0.5 0.5 0.5],'linewidth',2);
    pp.Color(4) = 0.9;

    y_ = uniqTemplateNumber(i)*ones(lf(i),1);
    ss = scatter(t_,y_,5*exp(log10(amp_)),log10(amp_),'filled');
    ss.MarkerFaceAlpha = 0.8;
    ss.MarkerEdgeColor = 'k';
    ss.MarkerEdgeAlpha = 0.8;
    ss.LineWidth = 0.2;
end
colorbar;
zoom on;