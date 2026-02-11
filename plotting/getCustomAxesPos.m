function ax = getCustomAxesPos(rows,cols,space)
if nargin < 1; rows = 1; end
if nargin < 2; cols = 1; end
if nargin < 3; space = 0.01; end %0.02 = 2% (percent)


%%

% Finds the dimensions of the figure as a whole by creating one big axes
% (which we will delete later) when we no longer need it.
BigAx = newplot;
set(BigAx,'Visible','off','color','none')
pos = get(BigAx,'Position');

% Calcuates the width and height of each subplot.
% Also where you can set the percent spacing between axes.
pos(1) = 0.10;
pos(2) = 0.10;
pos(3) = 0.97;
pos(4) = 0.89;
width = pos(3)/cols;             % width of each subplot
height = pos(4)/rows;            % height of each subplot

% Deleting the big axes
BigAxParent = get(BigAx,'Parent');
delete(BigAx);

% Makes handles to each of the subplot axes. Sets the position of each
% subplot to eliminate all the gray spaces between them. The handles to
% each subplot can be found in matrix 'ax'.
ax = gobjects(rows,cols);
for i=rows:-1:1
    for j=cols:-1:1
        axPos = [pos(1)+(j-1)*width,pos(2)+(rows-i)*height,width*(1-space),height*(1-space)];
        ax(i,j) = axes('Position',axPos,'parent',BigAxParent);
    end
end
