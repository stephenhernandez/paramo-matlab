function ax = linkedPlot(abscissa,ordinate,symbol,tile_spacing,line_width)
if nargin < 3
    symbol = '.';
end

if nargin < 4
    tile_spacing = "none";
end

if nargin < 5
    line_width = 2;
end

%%
[ny,nSubplots] = size(ordinate);
if isempty(abscissa)
    abscissa = (0:ny-1)';
end
nx = size(abscissa,1);

if nx ~= ny
    fprintf("Error, x and y do not have same size\n");
    return;
end

figure('units','normalized','outerposition',[0 0.1 0.5 0.9]);
tiledlayout(nSubplots,1,"Padding","compact","TileSpacing",tile_spacing);

ax = gobjects(nSubplots,1);
n = 0;
for i = 1:nSubplots
    n = n+1;
    ax(n,1) = nexttile;
    thisY = ordinate(:,i);
    plot(abscissa,thisY,symbol,"LineWidth",line_width); grid on;
end
linkaxes(ax,'x'); zoom on;

%%
if nSubplots < 2
    return;
end

%%
for i = 1:nSubplots
    if ~mod(i,2)
        ax(i).YAxisLocation = "right";
    end
    if i < nSubplots
        ax(i).XTickLabel = [];
    end
end