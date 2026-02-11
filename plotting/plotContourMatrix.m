function [xend,yend,xstart,ystart] = plotContourMatrix(ax,C,plotFlag)
if nargin < 3; plotFlag = true; end

c = size(C,2);
maxPoints = c;
n = 1;
cPos = ones(maxPoints,1);
level = cPos;
count = cPos;

while cPos(n) <= c
    level(n) = C(1,cPos(n));
    count(n) = C(2,cPos(n));
    n = n+1;
    cPos(n) = cPos(n-1)+count(n-1)+1;
end

n = n-1;
xend = zeros(n,1);
yend = xend;
xstart = xend;
ystart = xend;

%%
cPos = cPos(1:n);
level = level(1:n);
count = count(1:n);

%%
uniqueLevels = unique(level);
lul = length(uniqueLevels);
cmap = bone(lul);
cOld = uniqueLevels(1);
cTmp = 1;

for i = 1:n
    cLevel = level(i);
    if cOld ~= cLevel
        cOld = cLevel;
        cTmp = cTmp + 1;
    end
    xtmp(1:count(i),1) = C(1,cPos(i)+1:cPos(i)+count(i));
    ytmp(1:count(i),1) = C(2,cPos(i)+1:cPos(i)+count(i));
    if plotFlag
        hold(ax,'on');
        plot(ax,xtmp,ytmp,'--','color',cmap(cTmp,:),'linewidth',0.1);
    end

    xstart(i) = xtmp(1);
    ystart(i) = ytmp(1);
    xend(i) = xtmp(end);
    yend(i) = ytmp(end);
    clear xtmp ytmp
end