clear; close all;
wdir = '~/research/now/nyiragongo/';

%%
tStart = datetime(2021,04,01);
tEnd = datetime(2021,05,26);

%%
kstnm = 'MBAR';
knetwk = 'II';
khole = '00';
kcmpnms = ["BHZ";"BH1";"BH2";"BDF"];
lcmps = length(kcmpnms);

%% BEWARE %%
%
% BH1 = 17.5 degrees from north
% BH2 = 106.2 degrees from north 
% rotate accordingly...
%
%%

%%
dayRange = (tStart:tEnd)';
ldays = length(dayRange);
for i = 1:ldays
    tStart_ = dayRange(i);
    [yyyy,mm,dd] = datevec(tStart_);
    
    %%
    yyyyStr = num2str(yyyy);
    dayStr = num2str(dd);
    monthStr = num2str(mm);
    
    %%
    if dd < 10
        dayStr = ['0',dayStr];
    end
    
    if mm < 10
        monthStr = ['0',monthStr];
    end
    
    %%
    for j = 1:lcmps
        kcmpnm_ = char(kcmpnms(j));
        
        %%
        S_ = iris2sh(irisFetch.Traces(knetwk,kstnm,khole,kcmpnm_,datestr(tStart_,31),datestr(tStart_+1,31)));
        if isnat(S_(1).ref)
            fprintf('no data for day: %s, cmp: %s\n',datestr(tStart_),kcmpnm_);
            continue;
        end
        
        fName = [knetwk,'.',kstnm,'.',khole,'.',kcmpnm_,'.',yyyyStr,'.',monthStr,'.',dayStr,'.SAC'];
        fullFileName = fullfile(wdir,fName);
        sacwrite(fullFileName,S_);
    end
end