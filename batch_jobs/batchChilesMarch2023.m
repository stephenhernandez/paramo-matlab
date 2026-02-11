function batchChilesMarch2023(dayStart,dayEnd,dayInc)
template_file_name = "CVCCN_CHL1_March2023Templates_N233";
template_function_directory = fullfile("~","masa","template_search","template_functions");
load(fullfile(template_function_directory,template_file_name),"T");

lfc = 2;
hfc = 8;
newFs = 32;
madThresh = 10;
fileVersionNumber = 1;
custom_prefix = "CHL1";
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

%%
writeFlag = true;
plotFlag = false;
verboseFlag = true;
CANUSEGPU = canUseGPU();

%%
for i = 1:lDays
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        "CHL1",["HHZ";"HHN";"HHE"],"EC");
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
    Sf = syncWaveforms(Sf,false,true,true);
    Sf = nanGapWaveforms(Sf,0);
    Sf = padWaveforms(Sf);

    templateSearch(Sf,T,writeFlag,plotFlag,verboseFlag,madThresh,...
        custom_prefix,fileVersionNumber,CANUSEGPU);
end