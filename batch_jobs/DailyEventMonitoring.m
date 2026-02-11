% clear; close all;
% updateCotopaxiTremorCatalogBinaryClassifier(); toc;

%%
clear; close all;
dayStart = datetime(2022,10,01);
dayEnd = datetime(2026,02,10);
dayInc = 1;

%
tic;
days = (dayStart:dayInc:dayEnd)';
lDays = length(days);
for i = 1:lDays
    day_ = days(i);
    disp(day_);

    [M_,magErr_,t_,d_,Mcorr_] = loopCotopaxiTremorLocation(day_,day_,false,false); toc;
    load("~/masa/old/research/now/cotopaxi/CotopaxiTremorLocationResults_v5");
    tI = t >= day_+1;
    d(tI,:) = [];
    M(tI,:) = [];
    %McorrOrig(tI,:) = [];
    Mcorr(tI,:) = [];
    t(tI,:) = [];
    magErr(tI,:) = [];
    M = [M; M_];
    d = [d; d_];
    magErr = [magErr; magErr_];
    Mcorr = [Mcorr; Mcorr_];
    t = [t; t_];
    McorrOrig = Mcorr;
    clear Mcorr_ d_ M_ magErr_ t_;
    save("~/masa/old/research/now/cotopaxi/CotopaxiTremorLocationResults_v5","t","d","M","magErr","Mcorr","-v7.3");
    clearvars -except days dayStart dayEnd lDays dayInc i
end
clearvars -except dayStart dayEnd dayInc
toc;

% %%
% reventadorRepeaterSearch(dayStart,dayEnd,dayInc); toc;
% %cotopaxiRepeaterSearch(dayStart,dayEnd,dayInc); toc;
% batchChilesMarch2023(dayStart,dayEnd,dayInc); toc; %plotCVCCNBatchPlots();
% batchJob10(dayStart,dayEnd,dayInc); toc; %guagua
% PlataRepeaterSearch_v2(dayStart,dayEnd,dayInc); toc;
% CotoLPSearch(dayStart,dayEnd,dayInc); toc;
% SAG1_infrasound_array_processing_v2(dayStart,dayEnd,dayInc); toc;
% batchJob14(dayStart,dayEnd,dayInc); toc; %coto vlps
% 
% %% REVENTADOR
% thresh = 0.065;
% basisFunctionFileName = "~/igdata/ReventadorMultiplexedBasisFunctions_v1.mat";
% basisFunctions = load(basisFunctionFileName);
% reve_kstnms = ["REVN";"REVS";"CAYR";"CASC";"BONI";"YAHU";"IMBA";"ANTS";...
%     "ANTG";"CUSW";"URCU";"COTA";"CUIC";"TULM";"PULU";"GGPC";"GGPT"];
% reve_channels = ["HHZ";"HHN";"HHE";...
%     "BHZ";"BHN";"BHE";...
%     "SHZ";"SHN";"SHE"];
% writeFlag = true;
% maxEvents = 1e3;
% maxBasisFunctions = 20;
% linearccnorm = true;
% diffFlag = true;
% custom_prefix = "ReventadorMultiplexed_v1";
% 
% runSubspaceDetector(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
%     reve_kstnms,reve_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
%     diffFlag,custom_prefix); toc;
% 
% %% SANGAY
% thresh = 0.07;
% basisFunctionFileName = '~/igdata/SangayMultiplexedBasisFunctions_v2.mat';
% basisFunctions = load(basisFunctionFileName);
% sangay_kstnms = ["SAGA";"BPAT";"BMAS";"BRTU";"BULB";"BBIL";"BRUN";"PUYO";...
%     "TAMH";"PORT";"PKYU";"TAIS"];
% sangay_channels = ["HHZ";"HHN";"HHE";...
%     "BHZ";"BHN";"BHE"];
% writeFlag = true;
% maxEvents = 1e4;
% maxBasisFunctions = 20;
% linearccnorm = false;
% diffFlag = true;
% custom_prefix = "SangayMultiplexed_v1";
% 
% runSubspaceDetector(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
%     sangay_kstnms,sangay_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
%     diffFlag,custom_prefix); toc;
% 
% %% COTOPAXI
% thresh = 0.07;
% basisFunctionFileName = "~/igdata/CotopaxiMultiplexedBasisFunctions_v2.mat";
% basisFunctions = load(basisFunctionFileName);
% coto_kstnms = ["CO1V";"BREF";"BVC2";"BTAM"];
% coto_channels = ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"];
% 
% writeFlag = true;
% maxEvents = 1e4;
% maxBasisFunctions = 20;
% diffFlag = false;
% linearccnorm = true;
% custom_prefix = "CotopaxiMultiplexed_v1";
% 
% runSubspaceDetector(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
%     coto_kstnms,coto_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
%     diffFlag,custom_prefix); toc;
% 
% %%
% clear; close all; saveSangayCatalog(); SangayTimePredictability();
% GGP_swarm_parameters(); toc;
% plotCotopaxiVLP();
% 
% %%
% function [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
%     obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
%     runSubspaceDetector(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
%     kstnms,kcmpnms,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
%     diffFlag,customPrefix)
% 
% dayVec = (dayEnd:-dayInc:dayStart)';
% lDays = length(dayVec);
% for i = 1:lDays
%     tic;
%     day_ = dayVec(i);
%     S = loadWaveforms(day_,dayInc,...
%         kstnms,...
%         kcmpnms,...
%         "EC",""); %,true,true);
%     toc;
% 
%     badI = isnat(pull(S,'ref'));
%     S(badI) = [];
%     lS = length(S);
%     if ~lS
%         fprintf("not enough data for day: %s\n",day_);
%         continue;
%     end
% 
%     try
%         [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
%             obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
%             multiplexedSubspaceDetector(S,thresh,basisFunctions,...
%             maxBasisFunctions,maxEvents,writeFlag,diffFlag,linearccnorm,...
%             customPrefix);
%         toc;
%     catch
%         fprintf("something went wrong on day: %s\n",day_);
%         continue;
%     end
% end
% end