function batchJob10_2(dayStart,dayEnd,dayInc)
% clear; close all;
% dayStart = datetime(2024,08,26);
% dayEnd = datetime(2024,08,26);
% dayInc = 1;

cd ~/research/now/pichincha/pichincha_nxcorr/
load('GGPTemplateStructStack_7.mat','T');

%%
sizeT = size(T);
T = reshape(T,sizeT);
Torig = T;

lfc = 2;
hfc = 16;
newFs = 100;
thresh = 7;

dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

writeFlag = true;
plotFlag = ~true;
verboseFlag = true;

%%
minNpts = max(pull(Torig(:),'npts'));
for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        ["BTER";"PINO";"YANA";"JUA2";"GGPC";"GGPT"],...
        ["HHZ";"HHN";"HHE";"SHZ";"SHN";"SHE"],"EC");

    badI = false(size(S));
    kstnms = pull(S,'kstnm');
    kcmpnms = pull(S,'kcmpnm');
    pinoI = strcmp(kstnms,"PINO");
    lP = sum(pinoI);
    if lP
        sensorType = char(kcmpnms);
        sensorType = string(sensorType(:,1:2));
        if day_ >= datetime(2024,01,229)
            pinoHH = pinoI & strcmp(sensorType,"HH");
            lHH = sum(pinoHH);
            if lHH
                %
                pinoSH = pinoI & strcmp(sensorType,"SH");
                lSH = sum(pinoSH);
                if lSH
                    badI(pinoSH) = true; %clobber SH, but only if HH exists...
                end
                S_ = S(pinoHH);
                S_ = scaleWaveforms(S_,-1);
                lS_ = length(S_);
                for j = 1:lS_
                    kcmpnm_ = char(S_(j).kcmpnm);
                    S_(j).kcmpnm = string(['S',kcmpnm_(2:end)]);
                end
                S(pinoHH) = S_;
            end
        end
    end
    S(badI) = [];

    %%
    badI = isnat(pull(S,'ref')) | max(pull(S,'npts')) < minNpts;
    S(badI) = [];
    lS = length(S);
    if ~lS %min number of traces is 1
        fprintf("not enough data for day: %s\n",string(day_));
        continue;
    end

    %%
    Sf = differentiateWaveforms(S);
    Sf = padWaveforms(filterWaveforms(detrendWaveforms(Sf),lfc,hfc));

    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);

    customPrefix = "GGP";
    Tsncls = strcat(pull(T,'knetwk'),pull(T,'kstnm'),pull(T,'khole'),pull(T,'kcmpnm'));
    Ssncls = strcat(pull(Sf,'knetwk'),pull(Sf,'kstnm'),pull(Sf,'khole'),pull(Sf,'kcmpnm'));

    [lia,locb] = ismember(Ssncls,Tsncls);
    if sum(lia)
        [ccMain,tMain,ampMain,dMag,evidMain,magMain,templateNumber,madMain,nUsedMain,ccnorm,tLong] = ...
            templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,verboseFlag,thresh,customPrefix);
    end
end