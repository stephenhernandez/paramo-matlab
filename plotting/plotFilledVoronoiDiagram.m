%function
%plotFilledVoronoiDiagram
clear; close all; clc;

tags = ["BRTU","BMAS","BPAT","BBIL","BRUN","BULB","POND"];
%tags = ["SN01","SN02","SN03","SN04","SN05","SN06","SN07","SN08","SN09","SN10","SN11","SN12","SN13","SN14","VCH1","PVIL","CEAZ","ALCE","FER1","FER2"];
minLat=-1.54;
maxLat=-1.35;
minLon=-78.54;
maxLon=-78.35;

[stla,stlo,stel] = metaDataFromStationList(tags);
lS = length(stlo);
stla = [stla; minLat-0.2; minLat-0.2; maxLat+0.2; maxLat+0.2];
stlo = [stlo; minLon-0.2; maxLon+0.2; minLon-0.2; maxLon+0.2];
[v,c] = voronoin([stlo stla]);
cmap = parula(lS);
load ~/igdata/ec_boundaries.mat

%%
close all; figure();
hold on;
plot(lonEC,latEC,'k','linewidth',2); axis equal; zoom on;
for i = 1:lS%:length(stlo)
disp(i);
patch('Xdata',v(c{i},1),'YData',v(c{i},2),'FaceColor',cmap(i,:)); 
plot(stlo(i),stla(i),'d','markerfacecolor','w','markeredgecolor','k');
h = text(stlo(i),stla(i),tags(i),'Color','w');
end

load ~/igdata/dem/tungurahua_dem.mat
contour(lon(:),lat(:),demData,(2000:100:5000)','k','linewidth',1);
colorbar;

axis equal;
axis([minLon maxLon minLat maxLat]);

% hold on;
% plot(stlo,stla,'d'); for i = 1:length(stlo); text(stlo(i),stla(i),tags(i)); end;
% axis equal