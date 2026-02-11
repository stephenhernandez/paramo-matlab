function runSubspaceDetector(dayStart,dayEnd,dayInc)

%% REVENTADOR
thresh = 0.065;
basisFunctionFileName = '~/igdata/ReventadorMultiplexedBasisFunctions_v1.mat';
basisFunctions = load(basisFunctionFileName);
reve_kstnms = ["REVN";"REVS";"CAYR";"CASC";"BONI";"YAHU";"IMBA";"ANTS";...
    "ANTG";"CUSW";"URCU";"COTA";"CUIC";"TULM";"PULU";"GGPC";"GGPT"];
reve_channels = ["HHZ";"HHN";"HHE";...
    "BHZ";"BHN";"BHE";...
    "SHZ";"SHN";"SHE"];
writeFlag = true;
maxEvents = 1e3;
maxBasisFunctions = 20;
linearccnorm = true;
diffFlag = true;
customPrefix = 'ReventadorMultiplexed_v1';

runSubspaceDetector_(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
    reve_kstnms,reve_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
    diffFlag,customPrefix); toc;

%% SANGAY
thresh = 0.07;
basisFunctionFileName = '~/igdata/SangayMultiplexedBasisFunctions_v2.mat';
basisFunctions = load(basisFunctionFileName);
sangay_kstnms = ["SAGA";"BPAT";"BMAS";"BRTU";"BULB";"BBIL";"BRUN";"PUYO";...
    "TAMH";"PORT";"PKYU";"TAIS"];
sangay_channels = ["HHZ";"HHN";"HHE";...
    "BHZ";"BHN";"BHE"];
writeFlag = true;
maxEvents = 1e4;
maxBasisFunctions = 20;
linearccnorm = false;
diffFlag = true;
customPrefix = 'SangayMultiplexed_v1';

runSubspaceDetector_(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
    sangay_kstnms,sangay_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
    diffFlag,customPrefix); toc;

%% COTOPAXI
thresh = 0.07;
basisFunctionFileName = '~/igdata/CotopaxiMultiplexedBasisFunctions_v2.mat';
basisFunctions = load(basisFunctionFileName);
coto_kstnms = ["CO1V";"BREF";"BVC2";"BTAM"];
coto_channels = ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"];

writeFlag = true;
maxEvents = 1e4;
maxBasisFunctions = 20;
diffFlag = false;
linearccnorm = true;
customPrefix = 'CotopaxiMultiplexed_v1';

runSubspaceDetector_(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
    coto_kstnms,coto_channels,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
    diffFlag,customPrefix); toc;
end

%%
function [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
    obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
    runSubspaceDetector_(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
    kstnms,kcmpnms,writeFlag,maxEvents,maxBasisFunctions,linearccnorm,...
    diffFlag,customPrefix)

dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);
for i = 1:lDays
    tic;
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        kstnms,...
        kcmpnms,...
        "EC",""); %,true,true);
    toc;

    badI = isnat(pull(S,'ref'));
    S(badI) = [];
    lS = length(S);
    if ~lS
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end

    try
        [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
            obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
            multiplexedSubspaceDetector(S,thresh,basisFunctions,...
            maxBasisFunctions,maxEvents,writeFlag,diffFlag,linearccnorm,...
            customPrefix);
        toc;
    catch
        fprintf('something went wrong on day: %s',day_);
        continue;
    end
end
end
