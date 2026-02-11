function [t,eqlat,eqlon,eqdepth,eqmag,id,D,volcanoName] = generateSeismicityAnimation(tStart,tEnd,volcanoName,boundaryBox,SC3Flag)
if strcmp(volcanoName,'chilesExpanded')
    regionName = 'chilesExpanded';
    demName = 'chilesExpanded_dem'; %Expanded';
    minMag = 0.;
    maxDepth = 20;
    addSimbologia = true;
    allFrames = false;
    hillShadeFlag = true;
    magFact = 20;
    depthPanel = 1;
    cumPanel = 0;
elseif strcmp(volcanoName,'sierraNegra') %southern isabela island (SN+CA)
    regionName = 'sierraNegra'; %'galapagos';
    demName = volcanoName; %'southernIsabelaDEM';
    minMag = 2;
    maxDepth = 100;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 3;
    depthPanel = 1;
    cumPanel = 0;
elseif strcmp(volcanoName,'cuicocha') % cuicocha
    regionName = 'cuicocha';
    demName = 'cuicocha_dem'; %volcanoName; %'southernIsabelaDEM';
    minMag = 0;
    maxDepth = 20;
    addSimbologia = true;
    allFrames = false;
    hillShadeFlag = true;
    magFact = 20;
    depthPanel = true;
    cumPanel = true;
elseif strcmp(volcanoName,'atuntaqui') % cuicocha
    regionName = 'atuntaqui';
    demName = 'atuntaqui_dem'; %volcanoName; %'southernIsabelaDEM';
    minMag = 0;
    maxDepth = 50;
    addSimbologia = true;
    allFrames = false;
    hillShadeFlag = true;
    magFact = 20;
    depthPanel = true;
    cumPanel = true;
elseif strcmp(volcanoName,'sierra_negra') %sierraNegraCaldera
    regionName = volcanoName;
    demName = 'sierraNegra';
    minMag = 0.5;
    maxDepth = 20;
    addSimbologia = 0;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 4;
    depthPanel = 0;
    cumPanel = 1;
elseif strcmp(volcanoName,'southernIsabela') %sierraNegraCaldera
    regionName = volcanoName;
    demName = 'sierraNegra';
    minMag = 2.;
    maxDepth = 20;
    addSimbologia = true;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 3;
    depthPanel = 1;
    cumPanel = 0;
elseif strcmp(volcanoName,'fernandina') %sierraNegraCaldera
    regionName = volcanoName;
    demName = 'fernandina_dem';
    minMag = 1;
    maxDepth = 20;
    addSimbologia = true;
    allFrames = 0;
    hillShadeFlag = false;
    magFact = 8;
    depthPanel = true;
    cumPanel = false;
elseif strcmp(volcanoName,'quito_valley') %sierraNegraCaldera
    regionName = volcanoName;
    demName = 'quito_dem';
    minMag = 2;
    maxDepth = 30;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = true;
    magFact = 4;
    depthPanel = false;
    cumPanel = true;
elseif strcmp(volcanoName,'balao')
    regionName = 'balao';
    demName = 'balao_dem';
    minMag = 1;
    maxDepth = 100;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 5;
    depthPanel = 1;
    cumPanel = 0;
elseif strcmp(volcanoName,'pichincha')
    regionName = 'pichincha';
    demName = 'pichincha_dem';
    minMag = 0.;
    maxDepth = 20;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 5;
    depthPanel = true;
    cumPanel = false;
elseif strcmp(volcanoName,'cotopaxi')
    regionName = 'cotopaxi';
    demName = 'cotopaxi_dem';
    minMag = 2.;
    maxDepth = 20;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = false;
    magFact = 5;
    depthPanel = 1;
    cumPanel = 0;
elseif strcmp(volcanoName,'esmeraldas')
    regionName = 'esmeraldas';
    demName = 'esmeraldas_dem';
    minMag = 2.;
    maxDepth = 55;
    addSimbologia = false;
    allFrames = true;
    hillShadeFlag = true;
    magFact = 10;
    depthPanel = 1;
    cumPanel = 0;    
elseif strcmp(volcanoName,'tungurahua')
    regionName = 'tungurahua';
    demName = 'tungurahua_dem';
    minMag = 2.;
    maxDepth = 20;
    addSimbologia = false;
    allFrames = false;
    hillShadeFlag = 0;
    magFact = 5;
    depthPanel = 1;
    cumPanel = 0;
% else
%     disp(strcat('volcano name: ',volcanoName,' not found, sorry'));
%     return;
end

% get the job done
[t,eqlat,eqlon,eqdepth,eqmag,id,D] = ...
    generateSeismicityAnimationFrames(tStart,tEnd,regionName,boundaryBox,SC3Flag,...
    minMag,depthPanel,cumPanel,maxDepth,addSimbologia,allFrames,...
    hillShadeFlag,magFact);
