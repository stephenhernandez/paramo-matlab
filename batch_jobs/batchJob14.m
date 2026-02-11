function batchJob14(dayStart,dayEnd,dayInc)
template_function_directory = fullfile("~","masa","template_search","template_functions");
T = load(fullfile(template_function_directory,"FortyCotopaxiVLPs.mat"),"T");
T = T.T';

%%
rename_cmps = ["BHZ";"BHN";"BHE"];
old_cmps = ["HHZ";"HHN";"HHE"];
for j = 1:length(old_cmps)
    old_ = old_cmps(j);
    new_ = rename_cmps(j);
    hh_kcmpnmI = strcmp(pull(T,"kcmpnm"),old_);
    if ~sum(hh_kcmpnmI)
        fprintf("no %s\n",old_);
        continue;
    end
    T_ = push(T(hh_kcmpnmI),"kcmpnm",repmat(new_,sum(hh_kcmpnmI,"all"),1),true);
    T(hh_kcmpnmI) = T_;
end

%%
lfc = 1/5;
hfc = 1/2;
newFs = 10;
madThresh = 8;

%%
if dayEnd < dayStart
    fprintf("dayEnd should be dayStart, fix!\n");
    return;
end
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

writeFlag = true;
plotFlag = false;
verboseFlag = true;
custom_prefix = "vlps_cotopaxi";
CANUSEGPU = canUseGPU();
fileVersionNumber = 1;

%%
MINNPTS = max(pull(T,"npts"),[],"all");
Tsncls = strcat(pull(T(:,1),"knetwk"),pull(T(:,1),"kstnm"),...
    pull(T(:,1),"khole"),pull(T(:,1),"kcmpnm"));
chanList = ["BHZ";"BHN";"BHE";"HHZ";"HHN";"HHE"];
for i = 1:lDays
    day_ = dayVec(i);
    if day_ >= datetime(2025,01,23)
        trial_kstnms = ["CO1V"; "BTAM"; "BVC2"];
    else
        trial_kstnms = ["BREF";"BTAM";"BVC2"];
    end
    S = loadWaveforms(day_,dayInc,trial_kstnms,chanList,"EC","",false,false);

    %%
    badI = isnat(pull(S,"ref")) | max(pull(S,"npts")) < MINNPTS;
    S(badI) = [];
    lS = length(S);
    if ~lS %min number of traces is 1
        fprintf("not enough data for day: %s\n",string(day_));
        continue;
    end

    Sf = syncWaveforms(S,false,true,true);
    Sf = filterWaveforms(detrendWaveforms(Sf),lfc,hfc);
    Sf = resampleWaveforms(Sf,newFs);
    for j = 1:length(old_cmps)
        old_ = old_cmps(j);
        new_ = rename_cmps(j);
        hh_kcmpnmI = strcmp(pull(Sf,"kcmpnm"),old_);
        if ~sum(hh_kcmpnmI)
            fprintf("no %s\n",old_);
            continue;
        end
        S_ = push(Sf(hh_kcmpnmI),"kcmpnm",repmat(new_,sum(hh_kcmpnmI,"all"),1),true);
        Sf(hh_kcmpnmI) = S_;
    end

    Sf = mergeWaveforms(Sf);
    Sf = syncWaveforms(Sf,false,true,true);
    Sf = nanGapWaveforms(Sf,0);
    Sf = padWaveforms(Sf);

    Ssncls = strcat(pull(Sf,"knetwk"),pull(Sf,"kstnm"),...
        pull(Sf,"khole"),pull(Sf,"kcmpnm"));
    [lia,locb] = ismember(Ssncls,Tsncls);
    dataExist = sum(lia);
    if ~dataExist
        continue;
    end

    templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,...
        verboseFlag,madThresh,custom_prefix,fileVersionNumber,CANUSEGPU);
end