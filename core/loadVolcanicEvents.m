function [tsheet,type,Tsp,amp,magnitude,codaDuration,period,refStation] = ...
    loadVolcanicEvents(volcanoName,versionNumber)
if nargin < 2; versionNumber = '6'; end
volcanoName = lower(volcanoName);
if strcmp(volcanoName,'guagua') || strcmp(volcanoName,'pichincha')
    volcanoName = 'pichincha';
    volcanoSheet = ['V.GuaguaP_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'coto') || strcmp(volcanoName,'cotopaxi')
    volcanoName = 'cotopaxi';
    volcanoSheet = ['V.Cotopaxi_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'cayambe') || strcmp(volcanoName,'caya')
    volcanoName = 'cayambe';
    volcanoSheet = ['V.Cayambe_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'tungu') || strcmp(volcanoName,'tungurahua')
    volcanoName = 'tungurahua';
    volcanoSheet = ['V.Tungurahua_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'cota') || strcmp(volcanoName,'cuic') || strcmp(volcanoName,'cotacachi') || strcmp(volcanoName,'cuicocha')
    volcanoName = 'cota';
    volcanoSheet = ['V.Cota.Cuic_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'saga') || strcmp(volcanoName,'sangay')
    volcanoName = 'sangay';
    volcanoSheet = ['V.Sangay_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'chiles') || strcmp(volcanoName,'cerronegro')
    volcanoName = 'chiles';
    volcanoSheet = ['V.CCN_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'antisana') || strcmp(volcanoName,'anti')
    volcanoName = 'antisana';
    volcanoSheet = ['V.Antisana_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'reve') || strcmp(volcanoName,'reventador')
    volcanoName = 'reventador';
    volcanoSheet = ['V.Reventador_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'chimbo') || strcmp(volcanoName,'chimborazo')
    volcanoName = 'chimborazo';
    volcanoSheet = ['V.Chimborazo_v',num2str(versionNumber),'.xlsx'];
elseif strcmp(volcanoName,'galapagos')
    volcanoName = 'galapagos';
    volcanoSheet = ['V.Galapagos_v',num2str(versionNumber),'.xlsx'];
else
    disp('volcano not found, sorry.');
    return;
end

%%
spreadsheet_dir = fullfile("~","masa","spreadsheets");
fileName = fullfile(spreadsheet_dir,volcanoSheet);
% [num,txt] = xlsread(fileName,'Ev.NOLoc');
% if ismac
%     tsheet = datetime(datestr(num(:,1)+693960)); %695422));
% else
%     tsheet = datetime(datestr(num(:,1)+693960));
% end
%FechaHora               Fecha        Tipo     T_S_P_    AmpM_x    AmpM_n    AmpUnidad    Coda      Periodo_seg_       Magnitud    Energia    Estacion    Evol_dr_cm2    Evol_dr_sm2       Volcan
%tsheet,type,Tsp,amp,magnitude,codaDuration,period,refStation

%%
% [tsheet,sI] = sort(tsheet);
% txt = txt(sI,:);
% num = num(sI,:);
if strcmp(volcanoName,'tungurahua') || strcmp(volcanoName,'cotopaxi') ...
        || strcmp(volcanoName,'reventador') || strcmp(volcanoName,'pichincha') ...
        || strcmp(volcanoName,'sangay') || strcmp(volcanoName,'guagua') ...
        || strcmp(volcanoName,'galapagos') || strcmp(volcanoName,'chiles') ...
        || strcmp(volcanoName,'antisana') || strcmp(volcanoName,'cota') ...
        || strcmp(volcanoName,'chimborazo')

    T = readtable(fileName,'Sheet','Ev.NOLoc');
    type = string(T.Tipo);
    tsheet = T.FechaHora;
    Tsp = T.T_S_P_;
    amp = T.AmpM_x;
    magnitude = T.Magnitud;
    codaDuration = T.Coda;
    period = T.Periodo_seg_;
    refStation = T.Estacion;
else
    Tsp = num(:,4);
    amp = num(:,5);
    codaDuration = num(:,6);
    period = num(:,7);
    magnitude = num(:,8);
    type = string(txt(:,3));
    refStation = string(txt(:,10));
end
