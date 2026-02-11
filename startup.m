function startup(latexFlag)
if nargin < 1
    latexFlag = false;
end

%%
cd ~/matlab/
PWD = pwd;
PWD = strcat(PWD,'/working/');

setenv('TZ','America/Guayaquil');
datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
javaaddpath([PWD,'irisJava/IRIS-WS-2.0.15.jar']);
javaaddpath([PWD,'matTaup/lib/matTaup.jar']);

format long g;
addpath(genpath(PWD),'-end');
disp('running startup script, setting up custom default environment.');
rng('shuffle');

cd ~;
set(groot,'defaultAxesFontSize',24);
set(groot,'defaultColorbarFontSize',20);
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultScatterLineWidth',1);
set(groot,'defaultStemLineWidth',1);
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultColorbarTickLabelInterpreter','none');
set(groot,'defaultTextarrowshapeInterpreter','none');
set(groot,'defaultTextboxshapeInterpreter','none');
set(groot,'defaultTextFontName','Helvetica');
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultLineMarkerSize',14);
set(groot,'defaultFigurePaperPositionMode','auto');
set(groot,'defaultFigureRenderer','opengl');
set(groot,'defaultFigureColor',[1 1 1]);
set(groot,'defaultAxesTitleFontWeight','normal')

if latexFlag
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultColorbarTickLabelInterpreter','latex');
    set(groot,'defaultTextarrowshapeInterpreter','latex');
    set(groot,'defaultTextboxshapeInterpreter','latex');
end