clear; close all;

[center_lat,center_lon] = get_region_dimensions("cotopaxi");
horSwing = 8; %in km;
horInc = 0.4; %in km
zInc = 400; %in meters

% horSwing = 0.5; %in km;
% horInc = 0.05; %in km
% zInc = 50; %in meters
%horSwing = 10;
%horInc = 0.200;
%zInc = 100;

maxLat = center_lat + km2deg(horSwing);
minLat = center_lat - km2deg(horSwing);
maxLon = center_lon + km2deg(horSwing);
minLon = center_lon - km2deg(horSwing);

%%
minElev = -10000;
maxElev = 5900;

xvec = (minLon:km2deg(horInc):maxLon)';
yvec = (minLat:km2deg(horInc):maxLat)';
zvec = (maxElev:-zInc:minElev)';
[X,Y,Z] = meshgrid(xvec,yvec,zvec);

%%
[demLon,demLat,demData] = cutDEM([minLon maxLon minLat maxLat]); %local dem
X = X(:);
Y = Y(:);
Z = Z(:);

demData = demData';
uniqX = unique(X);
uniqY = unique(Y);

for i = 1:length(uniqX)
    x_ = uniqX(i);
    for j = 1:length(uniqY)
        y_ = uniqY(j);
        vq = interp2(demLon,demLat,demData,x_,y_);
        zI = X == x_ & Y == y_;
        zAll = Z(zI);
        vI = zAll < vq;
        zAll(~vI) = NaN;
        Z(zI) = zAll;
    end
end

%%
badI = ~isfinite(Z);
X(badI) = [];
Y(badI) = [];
Z(badI) = [];

%%
refEllipse = referenceEllipsoid('wgs84');
Model = load('~/research/now/cotopaxi/CotopaxiTremorAttenuationModel_v30012023.mat');
kstnmsOrig = Model.kstnmsOrig;
[stla,stlo,stel] = metaDataFromStationList(kstnmsOrig);

lX = length(X);
lStations = length(kstnmsOrig);
uniqX = unique(X);
uniqY = unique(Y);

distsCoarse = NaN(lX,lStations);
for i = 1:lStations
    stla_ = stla(i);
    stlo_ = stlo(i);
    stel_ = stel(i);
    for j = 1:length(uniqX)
        x_ = uniqX(j);
        for k = 1:length(uniqY)
            y_ = uniqY(k);
            horI = X == x_ & Y == y_;
            d_ = distance(y_,x_,stla_,stlo_,refEllipse)*1e-3;
            distsCoarse(horI,i) = sqrt(d_.^2 + ((Z(horI) - stel_)/1000).^2);
        end
    end
end

xvec = (-horInc:0.1:horInc)';
yvec = (-horInc:0.1:horInc)';
zvec = (zInc:-100:-zInc)';
[Xf,Yf,Zf] = meshgrid(xvec,yvec,zvec);

% finer grid
Xf = km2deg(Xf(:)); % deviation in dec. deg.
Yf = km2deg(Yf(:)); % deviation in dec. deg.
Zf = Zf(:);         % deviation in elevation in meters

clearvars -except X Y Z distsCoarse Zf Xf Yf stla stlo stel;
save("Cotopaxi3DCoarseTremorModel_2Steps");

%%
% figure();
% SS1 = scatter3(deg2km(X(:)-center_lon),deg2km(Y(:)-center_lat),Z(:)/1000,...
%     [],Z(:)/1000,'filled');
% zoom on; grid on; hold on;
% 
% [demLon2,demLat2] = meshgrid(demLon,demLat);
% SS2 = scatter3(deg2km(demLon2(:)-center_lon),deg2km(demLat2(:)-center_lat),...
%     demData(:)/1000,[],demData(:)/1000,'filled');
% axis equal;
% SS2.MarkerFaceAlpha = 0.85;
% xlabel("Easting");
% ylabel("Northing");
% 
% figure(1);
% plot3(deg2km(stlo-center_lon),deg2km(stla-center_lat),stel/1000,'kv');
