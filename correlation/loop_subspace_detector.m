function [ll,tabs,energyRatio,Neff,z2p,p2rms,kurt,skew,ReconstructCoefficients,...
    obsAmpRatio,synthAmpRatio,t,eratio,eratioOrig,theseU,dLong] = ...
    loop_subspace_detector(dayStart,dayEnd,dayInc,thresh,basisFunctions,...
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
        "EC","",false,false);
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