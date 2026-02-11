clear; close all; clc;

regionName = 'sierra_negra';
tStart = datetime(2018,04,20);
tEnd = datetime(2018,08,31);
regionName = char(regionName);
boundaryBox = [-91.3 -91.0 -0.95 -0.7];
minMag = 1;
scflag = false; %true is orig cat., false is relocated cat.
diasFlag = true;
magFact = 2;
depthPanel = false;

[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,depthPanel,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;

%%
tStart = datetime(2018,04,20);
tEnd = datetime(2018,06,26,09,00,00);
boundaryBox = [-91.2 -91.08 -0.88 -0.76];
[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,depthPanel,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;

%%
tStart = datetime(2018,06,26,09,00,00);
tEnd = datetime(2018,06,26,20,00,00);

[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,depthPanel,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;

%%
tStart = datetime(2018,06,26,20,00,00);
tEnd = datetime(2018,06,30);

[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,depthPanel,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;

%%
tStart = datetime(2018,06,30);
tEnd = datetime(2018,08,31);

[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,depthPanel,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;

%%
tStart = datetime(2018,06,24);
tEnd = datetime(2018,06,30);
boundaryBox = [-91.155 -91.12 -0.855 -0.83];
[~,~,~,~,~,id] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,scflag,...
    minMag,true,true,1000,false,false,false,magFact,18,'-djpeg',[],diasFlag);
zoom on;