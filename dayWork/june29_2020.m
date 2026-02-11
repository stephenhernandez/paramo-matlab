%clear; close all; clc;
cd ~/matlab/
PWD = pwd;
PWD = strcat(PWD,'/working/');

setenv('TZ','America/Guayaquil');
datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
javaaddpath([PWD,'/IRIS-WS-2.0.15.jar']);
javaaddpath([PWD,'/matTaup/lib/matTaup.jar']);
format long g;
addpath(genpath(PWD),'-end');
rng('shuffle');

set(groot,'defaultAxesFontSize',18);
set(groot,'defaultColorbarFontSize',12);
set(groot,'defaultLineLineWidth',0.1);
set(groot,'defaultScatterLineWidth',0.1);
set(groot,'defaultStemLineWidth',0.1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
set(groot,'defaultTextarrowshapeInterpreter','latex');
set(groot,'defaultTextboxshapeInterpreter','latex');
set(groot,'defaultLineMarkerSize',8);
set(groot,'defaultFigurePaperPositionMode','auto');
set(groot,'defaultFigureRenderer','Painters');
set(groot,'defaultFigureColor',[1 1 1]);

%%
cd ~/products/rsam/;

%%
S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"JUI6","SHZ",0,1,"EC","",0,0);
save('EC.JUI6..SHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"CAMI","SHZ",0,1,"EC","",0,0);
save('EC.CAMI..SHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"SAGA","HHZ",0,1,"EC","",0,0);
save('EC.SAGA..HHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BREF","BHZ",0,1,"EC","",0,0);
save('EC.BREF..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BTAM","BHZ",0,1,"EC","",0,0);
save('EC.BTAM..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BBIL","BHZ",0,1,"EC","",0,0);
save('EC.BBIL..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BRUN","BHZ",0,1,"EC","",0,0);
save('EC.BRUN..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BMAS","BHZ",0,1,"EC","",0,0);
save('EC.BMAS..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BPAT","BHZ",0,1,"EC","",0,0);
save('EC.BPAT..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BULB","BHZ",0,1,"EC","",0,0);
save('EC.BULB..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BVC2","BHZ",0,1,"EC","",0,0);
save('EC.BVC2..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BVC2","HHZ",0,1,"EC","",0,0);
save('EC.BVC2..HHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BNAS","BHZ",0,1,"EC","",0,0);
save('EC.BNAS..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BNAS","HHZ",0,1,"EC","",0,0);
save('EC.BNAS..HHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"BMOR","BHZ",0,1,"EC","",0,0);
save('EC.BMOR..BHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');

S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,10,-inf,"POND","HHZ",0,1,"EC","",0,0);
save('EC.POND..HHZ_10Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat');
