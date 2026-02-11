clear; close all;
cd ~/igdata/INSAR_SN/

T = readtable("20241020.txt");
amp = T.Var3;
lon = T.Var1;
lat = T.Var2;
goodI = isfinite(amp);

[C,IA,IC] = unique([lon(goodI) lat(goodI) amp(goodI)],'rows');
lon = C(:,1); 
lat = C(:,2); 
amp = C(:,3);

minLon = min(lon);
minLat = min(lat);
maxLon = max(lon);
maxLat = max(lat);

dxyDegrees = 0.05;
[Xq,Yq] = meshgrid((minLon:dxyDegrees:maxLon)',...
    (minLat:dxyDegrees:maxLat)');
dxyMeters = 111.317*1000*dxyDegrees;

dxyDegreesOrig = median([median(diff(lonG)),median(diff(latG))]);
dxyMetersOrig = 111.317*1000*dxyDegreesOrig;

minDat = min(min(SO2Smooth,[],1,"omitnan"),[],2,"omitnan");
SO2Smooth = 1 + SO2Smooth - minDat;
dI = ~isfinite(SO2Smooth);
SO2Smooth(dI) = 1;

[r,c] = meshgrid(lonG,latG);
r = r';
c = c';
SO2Smooth = SO2Smooth'; %SO2Orig';
F = griddedInterpolant(r,c,SO2Smooth,"makima","nearest");
Xq = Xq';
Yq = Yq';

SO2Interp = F(Xq,Yq);

close all; 
figure(); scatter(lon,lat,200,amp,'filled'); zoom on; colorbar;


% %
% T = readtable("20241113.txt");
% amp = T.Var3;
% lon = T.Var1;
% lat = T.Var2;
% goodI = isfinite(amp);
% [C,IA,IC] = unique([lon(goodI) lat(goodI) amp(goodI)],'rows');
% lon = C(:,1); lat = C(:,2); amp = C(:,3);
% %
% figure(); scatter(lon,lat,200,amp,'filled'); zoom on; colorbar;
% T = readtable("20241207.txt");
% amp = T.Var3;
% lon = T.Var1;
% lat = T.Var2;
% goodI = isfinite(amp);
% [C,IA,IC] = unique([lon(goodI) lat(goodI) amp(goodI)],'rows');
% lon = C(:,1); lat = C(:,2); amp = C(:,3);
% figure(); scatter(lon,lat,200,amp,'filled'); zoom on; colorbar;