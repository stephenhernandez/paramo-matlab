clear; close all; clc;

X0 = 797590;
Y0 = 9879600;

dh = 200;
dz = 100;

%
maxDepth = 6000;
maxWidth = 12000;
nRows = maxWidth/dh+1;
nColumns = maxDepth/dz + 1;

%
theta = 35;
dx = dh*sind(theta);
dy = dh* cosd(theta);

X = NaN(nRows,1);
Y = X;
X(1) = X0;
Y(1) = Y0;
for i = 2:nRows
    X(i) = X(i-1) - dx;
    Y(i) = Y(i-1) - dy;
end

%
figure(); plot(X,Y,'.'); zoom on; grid on; axis equal;

%
filasX = NaN(nRows,nColumns);
filasY = filasX;
filasX(:,1) = X;
filasY(:,1) = Y;

%%
newDX = dz*cosd(35)/tand(60);
newDY = dz*sind(35)/tand(60);
Z = NaN(nColumns,1);
Z(1) = 0;
for j = 2:nColumns
    filasX(:,j) = filasX(:,1) - (j-1)*newDX;
    filasY(:,j) = filasY(:,1) + (j-1)*newDY;
    Z(j) = Z(j-1)+dz;
end

for i = 2:nColumns
    figure(1); hold on;
    plot(filasX(:,i),filasY(:,i),'.');
end

%%
figure(); hold on;
Zfinal = filasX;
for i = 1:nColumns
    plot3(filasX(:,i),filasY(:,i),-repmat(Z(i),nRows,1),'.');
    Zfinal(:,i) = repmat(Z(i),nRows,1);
end
xlabel('Lon');
ylabel('Lat');
zlabel('Prof.');

%%
Xfinal = filasX(:);
Yfinal = filasY(:);
Zfinal = Zfinal(:);
