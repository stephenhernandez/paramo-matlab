clear; close all; 

tic;
ampFactor = 75;
[stlat2,stlon2] = metaDataFromStationList("SAGA");
refEllipse = referenceEllipsoid('wgs84');
boundaryBox = getRegionSpatialDimensions('sangay');

npoles = 1;
newFs = 1;
%Sorig = loadWaveforms(datetime(2021,06,03),2,"SAGA",["HHZ";"HHN";"HHE"],"EC",["";"01"]);
Sorig = loadWaveforms(datetime(2021,04,27),2,["SAGA"],["HHZ";"HHN";"HHE";"BDF"],"EC",["";"01"]);
if isnat(Sorig(1).ref)
    disp('no data');
    return;
end
Sorig = interpolateWaveforms(detrendWaveforms(differentiateWaveforms(Sorig,1)));



%close all;
% S = taperWaveforms(detrendWaveforms(resampleWaveforms(detrendWaveforms(Sorig),newFs)),0.001);
% S = scaleWaveforms(transferWaveforms(resampleWaveforms(S([1 2 3]),newFs),1/160,1/20,npoles,newFs,'vel'),1e9/9.81);
% S = rotateWaveforms(S,35,3,2,true);

%plotWaveforms(S);
toc;

%%
S = detrendWaveforms(cutWaveforms(Sorig,datetime(2021,04,27,22,15,00),0,minutes(120)));
S = syncWaveforms(S);
Srot = rotateWaveforms(S,35,3,2,true);

%%
boundaryBox = getRegionSpatialDimensions('sangay');
[fig,boundaryBox] = loadBasemap(boundaryBox,'turbo');
fig.Visible = 'on';
ax = gca;

%%
% S2 = loadWaveforms(datetime(2021,04,27),02,"SAGA",["HHZ";"HHN";"HHE"]);
% Scut = detrendWaveforms(cutWaveforms(S2,datetime(2021,04,27,20,00,00),0,hours(12)));
%Scut3 = detrendWaveforms(cutWaveforms(resampleWaveforms(resampleWaveforms(Scut,0.1),1),datetime(2021,04,27,22,15,00),0,minutes(30)));

% S2 = loadWaveforms(datetime(2021,05,07),01,"SAGA",["HHZ";"HHN";"HHE"]);
% Scut3 = detrendWaveforms(cutWaveforms(resampleWaveforms(resampleWaveforms(S2,0.25),8),...
%     datetime(2021,05,07,04,15,00),0,minutes(80)));

Z = S(1);
N = S(2);
E = S(3);
amps = sqrt(N.d.^2 + E.d.^2);
azi = 90 - atan2d(N.d,E.d);
[LatOut,LonOut] = reckon(stlat2,stlon2,ampFactor*amps,azi,refEllipse);

figure(1);
hold on;

tt = getTimeVec(E);
tt = seconds(tt - min(tt));

ax(2) = axes; %subplot(211);
ax(2).Parent = fig;

plot(ax(2),LonOut,LatOut,'.');
hold(ax(2),'on');
plot(ax(2),stlon2,stlat2,'v','markersize',20,'linewidth',4);
axis(ax,'equal');

%%
[LatOut2,LonOut2] = reckon(stlat2,stlon2,ampFactor*amps,35*ones(size(azi)),refEllipse);
plot(ax(2),LonOut2,LatOut2,'.');
zoom on;

ax(2).Visible = 'off';
linkaxes(ax,'xy');
axis(ax,boundaryBox);

grid on;
axis(ax,'equal');
axis(ax,boundaryBox);
plotWaveforms(Srot);
