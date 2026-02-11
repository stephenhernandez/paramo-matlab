function PlataRepeaterSearch_v2(dayStart,dayEnd,dayInc)
warning off MATLAB:polyfit:PolyNotUnique

template_file_name = "PlataRepeaterTemplates_v2";
template_function_directory = fullfile("~","masa","template_search","template_functions");
load(fullfile(template_function_directory,template_file_name),"TNew");

lfc = 1;
hfc = 4;
newFs = 20;
madThresh = 7;
fileVersionNumber = 1;
dayVec = (dayEnd:-dayInc:dayStart)';
lDays = length(dayVec);
npoles = 4;
custom_prefix = "isla_de_la_plata";

T = TNew;
Tsncls = strcat(pull(T(:,1),"knetwk"),pull(T(:,1),"kstnm"),...
    pull(T(:,1),"khole"),pull(T(:,1),"kcmpnm"));
Torig = T;
npts = pull(Torig,"npts");
MINNPTS = max(npts,[],"all");
chanList = ["BHZ";"BHN";"BHE";"HHZ";"HHN";"HHE"];
rename_cmps = ["BHZ";"BHN";"BHE"];
old_cmps = ["HHZ";"HHN";"HHE"];

%%
writeFlag = true;
plotFlag = false;
verboseFlag = true;
CANUSEGPU = canUseGPU();
DECONFLAG = true;
WAFLAG = false;
tic;
for i = 1:lDays
    day_ = dayVec(i);
    if day_ < datetime(2025,01,23)
        kstnms = ["BREF";"BNAS";"BVC2";"BTAM";"BMOR";"BMAS";"BPAT";"BULB";"BRUN";"BBIL"];
    else
        kstnms = ["BNAS";"BVC2";"BTAM";"BMOR";"BMAS";"BPAT";"BULB";"BRUN";"BBIL"];
    end
    S = loadWaveforms(day_,dayInc,kstnms,chanList,"EC"); toc;

    %%
    badI = isnat(pull(S,"ref")) | max(pull(S,"npts")) < MINNPTS;
    S(badI) = [];
    lS = length(S);
    if ~lS %min number of traces is 1
        fprintf("not enough data for day: %s\n",string(day_));
        continue;
    end

    %
    fprintf("starting deconvolution procedure....\n");
    Sf = resampleWaveforms(S,newFs); %<-- do not erase, makes transfer run faster by pre-downsampling
    Sf = transferWaveforms(Sf,lfc,hfc,npoles,newFs,"vel",DECONFLAG,WAFLAG,CANUSEGPU);
    Sf = scaleWaveforms(Sf,1e9); toc;
    fprintf("done with deconvolution....\n");
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

    fprintf("data exist (%d sncls) for day: %s\n",dataExist,day_);
    templateSearch(Sf(lia),T(locb(lia),:),writeFlag,plotFlag,verboseFlag,...
        madThresh,custom_prefix,fileVersionNumber,CANUSEGPU); toc;
    fprintf("\n");
end