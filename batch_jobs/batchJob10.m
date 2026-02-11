function batchJob10(dayStart,dayEnd,dayInc)
% template_file_name = "ggp_no_pino";
% template_function_directory = fullfile("~","masa","template_search","template_functions");
% load(fullfile(template_function_directory,template_file_name),"T");
% Tsncls = strcat(pull(T(:,1),"knetwk"),pull(T(:,1),"kstnm"),...
%     pull(T(:,1),"khole"),pull(T(:,1),"kcmpnm"));
% 
% %%
% lfc = 2;
% hfc = 8;
% newFs = 64;
% thresh = 7;
% dayVec = (dayEnd:-dayInc:dayStart)';
% lDays = length(dayVec);
% 
% %%
% writeFlag = true;
% plotFlag = false;
% verboseFlag = true;
% MINNPTS = max(pull(T,"npts"),[],"all");
% custom_prefix = "ggp_no_pino";
% CANUSEGPU = canUseGPU();
% fileVersionNumber = 1;
% for i = 1:lDays
%     day_ = dayVec(i);
%     S = loadWaveforms(day_,dayInc,...
%         ["BTER";"GGPC";"GGPT";"YANA"],...
%         ["SHZ";"SHN";"SHE";"HHZ";"HHN";"HHE"],"EC","",false,false);
% 
%     badI = isnat(pull(S,'ref')) | max(pull(S,'npts')) < MINNPTS;
%     S(badI) = [];
%     lS = length(S);
%     if ~lS %min number of traces is 1
%         fprintf("not enough data for day: %s\n",string(day_));
%         continue;
%     end
% 
%     %%
%     Sf = syncWaveforms(S,false,true,true);
%     Sf = filterWaveforms(detrendWaveforms(Sf),lfc,hfc);
%     Sf = resampleWaveforms(Sf,newFs);
%     Sf = nanGapWaveforms(Sf,0);
%     Sf = padWaveforms(Sf);
% 
%     Ssncls = strcat(pull(Sf,"knetwk"),pull(Sf,"kstnm"),...
%         pull(Sf,"khole"),pull(Sf,"kcmpnm"));
%     [lia,locb] = ismember(Ssncls,Tsncls);
%     dataExist = sum(lia);
%     if ~dataExist
%         continue;
%     end
%     templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,verboseFlag,...
%         thresh,custom_prefix,fileVersionNumber,CANUSEGPU);
% end
% 
% %%
% clearvars -except dayStart dayEnd dayInc
template_function_directory = fullfile("~","masa","template_search","template_functions");
load(fullfile(template_function_directory,"GGPTemplateStructStack_7.mat"),"T");
Tsncls = strcat(pull(T(:,1),"knetwk"),pull(T(:,1),"kstnm"),...
    pull(T(:,1),"khole"),pull(T(:,1),"kcmpnm"));

%
%JUA2 on 16-Nov-2023 is multiplexed and has data from other sncls and other
%years....
%

%%
lfc = 2;
hfc = 16;
newFs = 100;
thresh = 7;

dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);

%%
writeFlag = true;
plotFlag = false;
verboseFlag = true;
MINNPTS = max(pull(T,"npts"),[],"all");
custom_prefix = "ggp";
CANUSEGPU = canUseGPU();
fileVersionNumber = 1;
for i = 1:lDays
    tic;
    day_ = dayVec(i);
    S = loadWaveforms(day_,dayInc,...
        ["BTER";"PINO";"YANA";"JUA2";"GGPC";"GGPT"],...
        ["HHZ";"HHN";"HHE";"SHZ";"SHN";"SHE"],"EC");

    badI = isnat(pull(S,"ref")) | max(pull(S,"npts")) < MINNPTS;
    S(badI) = [];
    lS = length(S);
    if ~lS %min number of traces is 1
        fprintf("not enough data for day: %s\n",day_);
        continue;
    end
    kstnms = pull(S,"kstnm");
    kcmpnms = pull(S,"kcmpnm");
    kstnmI = strcmp(kstnms,"PINO");
    sensorType = char(kcmpnms);
    sensorType = string(sensorType(:,1:2));

    %%
    HHI = kstnmI & strcmp(sensorType,"HH");
    lH = sum(HHI);
    if day_ >= datetime(2024,01,229)
        %delete PINO.SH? channels
        badI = kstnmI & strcmp(sensorType,"SH");
        if lH % flip polarity of PINO.HH? channels (if they exist)
            S_ = S(HHI);
            S_ = scaleWaveforms(S_,-1);
            lS_ = length(S_);
            for j = 1:lS_
                kcmpnm_ = char(S_(j).kcmpnm);
                S_(j).kcmpnm = string(['S',kcmpnm_(2:end)]);
            end
            S(HHI) = S_;
        end
    else
        badI = HHI;
    end
    S(badI) = [];

    kstnmI = strcmp(kstnms,"YANA");
    HHI = kstnmI & strcmp(kcmpnms,"HHZ");
    lH = sum(HHI);
    if day_ >= datetime(2025,01,176)
        badI = kstnmI & strcmp(kcmpnms,"SHZ");
        if lH % YANA.HHZ exists
            S_ = S(HHI);
            %S_ = scaleWaveforms(S_,-1); %no polarity flips for YANA
            lS_ = length(S_);
            for j = 1:lS_
                kcmpnm_ = char(S_(j).kcmpnm);
                S_(j).kcmpnm = string(['S',kcmpnm_(2:end)]);
            end
            S(HHI) = S_;
        end
    else
        badI = HHI;
    end
    S(badI) = [];

    kstnmI = strcmp(kstnms,"JUA2");
    HHI = kstnmI & strcmp(kcmpnms,"HHZ");
    lH = sum(HHI);
    if day_ >= datetime(2025,01,156)
        badI = kstnmI & strcmp(kcmpnms,"SHZ");
        if lH % JUA2.HHZ exists
            S_ = S(HHI);
            %S_ = scaleWaveforms(S_,-1); % polarity unknown for JUA2, confirm later
            lS_ = length(S_);
            for j = 1:lS_
                kcmpnm_ = char(S_(j).kcmpnm);
                S_(j).kcmpnm = string(['S',kcmpnm_(2:end)]);
            end
            S(HHI) = S_;
        end
    else
        badI = HHI;
    end
    S(badI) = [];

    %%
    Sf = detrendWaveforms(S);
    Sf = differentiateWaveforms(Sf);
    Sf = syncWaveforms(Sf,false,true,true);
    Sf = filterWaveforms(detrendWaveforms(Sf),lfc,hfc);
    Sf = resampleWaveforms(Sf,newFs);
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
        verboseFlag,thresh,custom_prefix,fileVersionNumber,CANUSEGPU);
end