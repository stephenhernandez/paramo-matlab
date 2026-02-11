%function continuous_volcanic_monitoring()
clear; close all;
cd ~/matlab/
PWD = pwd;
PWD = strcat(PWD,'/working/');

setenv('TZ','America/Guayaquil');
datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
javaaddpath([PWD,'irisJava/IRIS-WS-2.20.1.jar']);
javaaddpath([PWD,'/matTaup/lib/matTaup.jar']);
format long g;
addpath(genpath(PWD),'-end');
rng('shuffle');

set(groot,'defaultAxesFontSize',22);
set(groot,'defaultColorbarFontSize',18);
set(groot,'defaultLineLineWidth',0.1);
set(groot,'defaultScatterLineWidth',0.1);
set(groot,'defaultStemLineWidth',0.1);
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultColorbarTickLabelInterpreter','none');
set(groot,'defaultTextarrowshapeInterpreter','none');
set(groot,'defaultTextboxshapeInterpreter','none');
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaultTextInterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultColorbarTickLabelInterpreter','latex');
% set(groot,'defaultTextarrowshapeInterpreter','latex');
% set(groot,'defaultTextboxshapeInterpreter','latex');
set(groot,'defaultTextFontName','Helvetica');
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultLineMarkerSize',8);
set(groot,'defaultFigurePaperPositionMode','auto');
%set(groot,'defaultFigureRenderer','Painters');
set(groot,'defaultFigureColor',[1 1 1]);


%% define snclGroups
% Sangay SNCL Group
snclGroups(1).snclList = [...
    "SAGA" "HHZ" "EC" "";...
    "SAGA" "BDF" "EC" "01";...
    "PUYO" "HHZ" "EC" "";...
    "TAIS" "HHZ" "EC" "";...
    "PKYU" "HHZ" "EC" "";...
    "TAMH" "HHZ" "EC" "";...
    "BPAT" "BHZ" "EC" "";...
    "BMAS" "HHZ" "EC" "";...
    "BULB" "HHZ" "EC" "";...
    "BRUN" "BHZ" "EC" "";...
    "PORT" "HHZ" "EC" "";...
    "SAG1" "BDF" "EC" "01";...
    "SAG1" "BDF" "EC" "02";...
    "SAG1" "BDF" "EC" "03";...
    "SAG1" "BDF" "EC" "04";...
    "SAG1" "BDF" "EC" "05";...
    "PIKA" "HHZ" "EC" "";...
    "SAG1" "HHZ" "EC" ""];
snclGroups(1).saveDir = 'sangay';
snclGroups(1).refEllipse = referenceEllipsoid('wgs84');
snclGroups(1).pauseTime = 2;
snclGroups(1).lfc = 0.6;
snclGroups(1).hfc = 1.2;
snclGroups(1).newFs = 20;
snclGroups(1).faceAlpha = 0.25;
snclGroups(1).templateFileName = '~/research/now/sangay/sangay_svd_basis_functions_withSNR_v8';
snclGroups(1).temporaryCatalogFile = '~/research/now/sangay/sangayTenSensorsUpdate_v4';
snclGroups(1).longTermCatalogFile = '~/research/now/sangay/SangayRegionalAnalysis_v7.mat';
snclGroups(1).shortSnippet = 60;
snclGroups(1).threshold = 0.2;
snclGroups(1).maxTemplates = 20;
snclGroups(1).recordLength = 150;
snclGroups(1).maxN = 2e3;
snclGroups(1).mpd = 30;
snclGroups(1).linearccnorm = true; %<-- must be true for Sangay for subspace detector!!!!!
snclGroups(1).plotFlag = false;
snclGroups(1).verboseFlag = false;
snclGroups(1).diffFlag = 0;
snclGroups(1).debuggingMode = false;
snclGroups(1).nHours = 8;
snclGroups(1).nDays = 1;

%% Reventador SNCL Group
snclGroups(2) = snclGroups(1); %copy previous group
snclGroups(2).snclList = [...
    "REVS" "HHZ" "EC" "";...
    "REVS" "BDF" "EC" "01";...
    "REVN" "HHZ" "EC" "";...
    "REVN" "BDF" "EC" "";...
    "BOSC" "HHZ" "EC" "";...
    "ANTS" "HHZ" "EC" "";...
    "CASC" "HHZ" "EC" "";...
    "CAYR" "SHZ" "EC" ""];
snclGroups(2).saveDir = 'reventador';
snclGroups(2).lfc = 0.25;
snclGroups(2).hfc = 1;
snclGroups(2).templateFileName = '~/research/now/reventador/reventador_svd_basis_functions'; %casc_svd_basis_functions_2';
snclGroups(2).temporaryCatalogFile = '~/research/now/reventador/revTmpUpdate_v2';
snclGroups(2).longTermCatalogFile = '~/research/now/reventador/ReventadorSubspaceDetectorResults_v6';
snclGroups(2).threshold = 0.12;
snclGroups(2).maxTemplates = 5;
snclGroups(2).recordLength = 120;
snclGroups(2).mpd = 20;

%% fernandina sncl group
snclGroups(3) = snclGroups(1); %copy previous group
snclGroups(3).snclList = [...
    "FER1" "BHZ" "EC" "";...
    "FER1" "BHN" "EC" "";...
    "FER1" "BHE" "EC" "";...
    "FER2" "HHZ" "EC" "";...
    "FER2" "HHN" "EC" "";...
    "FER2","HHE","EC","";...
    "PAYG" "BHZ" "IU" "00"];
snclGroups(3).saveDir = 'fernandina';
snclGroups(3).lfc = 4;
snclGroups(3).hfc = 12;
snclGroups(3).templateFileName = '~/research/now/fernandina/fernandina_svd_basisFunctions_v2.mat';
snclGroups(3).temporaryCatalogFile = '~/research/now/fernandina/ferTmpUpdate';
snclGroups(3).longTermCatalogFile = '~/research/now/fernandina/fernandinaSubspaceDetector_v2';
snclGroups(3).threshold = 0.1;
snclGroups(3).maxTemplates = 5;
snclGroups(3).recordLength = 40;
snclGroups(3).mpd = 10;
snclGroups(3).shortSnippet = 30;
snclGroups(3).maxN = 1e3;

%% Cotopaxi SNCL Group
snclGroups(4) = snclGroups(1); %copy previous group
snclGroups(4).saveDir = 'cotopaxi';
snclGroups(4).lfc = 3/4; %5/8;
snclGroups(4).hfc = 3; %5;
saveHome = "~/research/now/cotopaxi/";
templateFileName = fullfile(saveHome,"cotopaxi_svd_basis_functions_v5.mat");
temporaryCatalogFile = fullfile(saveHome,"cotoTmpUpdate_v4")';
longTermCatalogFile = fullfile(saveHome,"cotopaxiSubspaceDetectorBREF_v8");
load(templateFileName,'kstnm','chan');
snclList = [kstnm chan repmat("EC",length(kstnm),1) repmat("",length(kstnm),1)];
snclList = [snclList; ...
    "CO1V" "HHZ" "EC" "";...
    "COSE" "HHZ" "EC" ""];
snclGroups(4).snclList = snclList;
snclGroups(4).templateFileName = templateFileName;
snclGroups(4).temporaryCatalogFile = temporaryCatalogFile;
snclGroups(4).longTermCatalogFile = longTermCatalogFile;
snclGroups(4).threshold = 0.17;
snclGroups(4).maxTemplates = 15; %15;
snclGroups(4).recordLength = 35; %60;
snclGroups(4).mpd = 10;

%%
htmlDir = '~/public_html/helis/';

%%
dayStart = datetime(2025,03,25);
dayEnd = datetime(2025,03,31);
dayInc = 1;
% dayVec = (dayStart:dayInc:dayEnd)';
% lDays = length(dayVec);
% n = 0;
R = load('~/research/now/sangay/reprocessDays','reprocessDays');
dayVec = R.reprocessDays; %(dayStart:dayInc:dayEnd)';
lDays = length(dayVec);
%n = 0;
%while n < lDays
for day_ = dayVec'
    cd(htmlDir);
    for i = 1%:length(snclGroups)
        snclGroup = snclGroups(i);
        if i == 1
            fprintf("----------------\n");
            fprintf("PROCESSING SANGAY\n");
            fprintf("----------------\n");
        elseif i == 2
            fprintf("----------------\n");
            fprintf("PROCESSING REVENTADOR\n");
            fprintf("----------------\n");
        elseif i == 3
            fprintf("----------------\n");
            fprintf("PROCESSING FERNANDINA\n");
            fprintf("----------------\n");
        else
            fprintf("----------------\n");
            fprintf("PROCESSING COTOPAXI\n");
            fprintf("----------------\n");
        end

        try
            %n = n+1;
            %day_ = dayVec(n);
            fprintf("%s\n",day_);
            processSNCLGroup(snclGroup,day_);
            fprintf("\n");
            fprintf("\n");
        catch ME
            warning(ME.message);
            disp('for some reason an error was thrown, cant quite figure out why');
        end
    end

    %
    % try
    %     processHelicorders(); %(sncl_list,true);
    %     pause(5);
    %     updateGlobalCatalog();
    % catch
    %     fprintf("\n");
    %     fprintf("this iteration didnt work\n");
    %     fprintf("\n");
    %     continue;
    % end
end
