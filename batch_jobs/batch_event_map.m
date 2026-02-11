%function batchJob3(ids)
%if nargin < 1
%%
clear; close all; clc;
cd ~/matlab/
PWD = pwd;
PWD = strcat(PWD,'/working/');
setenv('TZ','America/Guayaquil');
datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
javaaddpath([PWD,'irisJava/IRIS-WS-2.20.1.jar']);
javaaddpath([PWD,'matTaup/lib/matTaup.jar']);
format long g;
addpath(genpath(PWD),'-end');
disp('running startup script, setting up custom default environment.');
rng('shuffle');
set(groot,'defaultAxesFontSize',22);
set(groot,'defaultColorbarFontSize',18);
set(groot,'defaultLineLineWidth',1);
set(groot,'defaultScatterLineWidth',1);
set(groot,'defaultStemLineWidth',1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
set(groot,'defaultTextarrowshapeInterpreter','latex');
set(groot,'defaultTextboxshapeInterpreter','latex');
set(groot,'defaultTextFontName','Helvetica');
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultLineMarkerSize',15);
set(groot,'defaultFigurePaperPositionMode','auto');
%set(groot,'defaultFigureRenderer','Painters');
set(groot,'defaultFigureColor',[1 1 1]);

%%
minMag = 6;
maxMag = 8;

% [evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
% timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
% locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1();
try
    [~,~,t,~,~,~,eqmag,ids,stderr,azgap,nPhases,nMLv,...
        ~,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,~,~,...
        ~,~,~,~,~,scbullmag] = readCat1();
catch ME
    warning(ME.message);
    try
        [~,~,t,~,~,~,eqmag,ids,stderr,azgap,nPhases,~,...
            ~,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,~,~,...
            ~,~,~,~,~,scbullmag] = readCat1();
    catch ME
        rethrow(ME);
    end
end

eqmag = round(eqmag*10)/10;
scbullmag = round(scbullmag*10)/10;

sError = sqrt(eqlaterr.^2 + eqlonerr.^2 + eqdeptherr.^2);
eI = eqmag <= maxMag & eqmag >= minMag & t > datetime(2025,01,31) & stderr > 0 & stderr < 2 & ...
    eqmagerr >= 0 & eqmagerr < 1 & sError <= 100 & azgap <= 360 & nMLv >= 8 & nPhases >= 8;

%%
eqmag = eqmag(eI);
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

[eqmag,eI] = sort(eqmag,'descend');
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

scbullmag(scbullmag < -100) = eqmag(scbullmag < -100);
eqmag = scbullmag;

eI = eqmag >= minMag & t > datetime(2010,01,01);

eqmag = eqmag(eI);
ids = ids(eI);

[~,eI] = sort(eqmag,'descend');
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

%%
npoles = 4;
hfc = 40;
mainDir = '~/products/EcuadorStrongMotion/';
vFlag = false;
saveRoot = fullfile('~','products','EcuadorStrongMotion');

for i = 1:length(ids)
    eventID_ = ids(i);
    yyyy = char(eventID_);
    yyyy = yyyy(6:9);
    files = dir(strcat(mainDir,yyyy,'/*/',eventID_,...
        '*/',eventID_,'.jpg'));
    lfiles = length(files);

    if lfiles > 0
        fprintf('already processed %s, mag: %f, skipping...\n',...
            eventID_,eqmag(i));
        continue;
    end

    fprintf('will process: <strong>%s</strong>\n',eventID_);
    E = readSCBulletin(eventID_);
    origmag = E.mag;
    origt = E.t;
    [yyyy,mm,dd] = datevec(origt);
    yyyyStr = sprintf("%04d",yyyy);
    mmStr = sprintf("%02d",mm);
    origmag = round(origmag*100)/100;
    eventDirName = sprintf("%s_M%3.2f_%s",eventID_,origmag,datestr(origt,30));
    writeDir = fullfile(saveRoot,yyyyStr,mmStr,eventDirName);
    if ~exist(writeDir,'dir')
        mkdir(writeDir);
    end
    plotEvents(E,vFlag,writeDir);
end