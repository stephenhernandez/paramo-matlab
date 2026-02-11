function CotoLPSearch(dayStart,dayEnd,dayInc)
template_file_name = "CotoLPTemplate";
template_function_directory = fullfile("~","masa","template_search","template_functions");
load(fullfile(template_function_directory,template_file_name),"T");

Tsncls = strcat(pull(T(:,1),"knetwk"),pull(T(:,1),"kstnm"),...
    pull(T(:,1),"khole"),pull(T(:,1),"kcmpnm"));
lfc = 0.5;
hfc = 8;
newFs = 20;
madThresh = 7;
custom_prefix = "lp_cotopaxi";
fileVersionNumber = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

%%
writeFlag = true;
plotFlag = false;
verboseFlag = true;
CANUSEGPU = canUseGPU();
for i = 1:lDays
    day_ = dayVec(i);
    if day_ >= datetime(2025,01,23)
        trial_kstnms = ["CO1V";"BTAM";"BVC2"];
    else
        trial_kstnms = ["CO1V";"BREF";"BTAM";"BVC2"];
    end

    S = loadWaveforms(day_,dayInc,...
        trial_kstnms,...
        ["HHZ";"HHN";"HHE";"BHZ";"BHN";"BHE"],"EC","",false,false);

    badI = isnat(pull(S,'ref'));
    S(badI) = [];
    lS = length(S);
    if ~lS
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end
    Sf = syncWaveforms(S,false,true,true);
    Sf = filterWaveforms(detrendWaveforms(Sf),lfc,hfc);
    Sf = resampleWaveforms(Sf,newFs);
    Sf = nanGapWaveforms(Sf,0);
    Sf = padWaveforms(Sf);

    Ssncls = strcat(pull(Sf,"knetwk"),pull(Sf,"kstnm"),...
        pull(Sf,"khole"),pull(Sf,"kcmpnm"));
    [lia,locb] = ismember(Ssncls,Tsncls);
    if ~sum(lia)
        continue;
    end
    S_ = Sf(lia);
    templateSearch(S_,T(locb(lia),:),writeFlag,plotFlag,verboseFlag,...
        madThresh,custom_prefix,fileVersionNumber,CANUSEGPU);
end