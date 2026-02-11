cd ~/products/rsam/
clear; close all; clc;

%%
dur = 30;
files = dir('EC*1Hz*_SLIDEMED_*');
files = flipud(files);

%%
for i = 1:length(files)
    try
        load(files(i).name,'S');
    catch
        disp('couldnt load data');
    end
    
    %%
    kstnm = S.kstnm;
    kcmpnm = S.kcmpnm;
    knetwk = S.knetwk;
    khole = S.khole;
    
    %%
    if strcmp(kcmpnm,'BDF')
        lfc = 1/4;
        hfc = 1;
    else
        lfc = 1;
        hfc = 8;
    end
    
    %%
    nowTime = dn2dt(now) + hours(5);
    tStart = datetime(1997,01,01);
    tEnd = dateshift(nowTime,'start','day');
    meanFlag = false; %when false, we take the median amplitude in the window
    rmsFlag = true;
    
    %%
    try
        S = rmsGather(tStart,tEnd,dur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag);
        S.d = single(S.d);
    catch
        disp('couldnt generate new data');
        continue;
    end
    
    %%
    if strcmp(kcmpnm,"ENZ") || strcmp(kcmpnm,"HNZ")
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
    elseif strcmp(kcmpnm,"BDF")
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
    else
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
    end
    
    %%
    disp(fname);
    save(fname,'S');
    
    %% 0.6 - 1.2 filter
    if ~(strcmp(kcmpnm,"BDF") || strcmp(kcmpnm,"HNZ") || strcmp(kcmpnm,"ENZ"))...
            && (strcmp(kstnm,"PUYO") || strcmp(kstnm,"TAIS") || strcmp(kstnm,"TAMH")...
            || strcmp(kstnm,"PORT") || strcmp(kstnm,"PKYU") || strcmp(kstnm,"BBIL")...
            || strcmp(kstnm,"BRUN") || strcmp(kstnm,"BPAT") || strcmp(kstnm,"BMAS")...
            || strcmp(kstnm,"BULB"))
        
        %%
        lfc = 0.6;
        hfc = 1.2;
        nowTime = dn2dt(now) + hours(5);
        tEnd = dateshift(nowTime,'start','day');
        meanFlag = false; %when false, we take the median amplitude in the window
        rmsFlag = true;
        
        %%
        try
            S = rmsGather(tStart,tEnd,dur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
            S.d = single(S.d);
        catch
            disp('couldnt generate new data');
            continue;
        end
        
        %%
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
        
        %%
        disp(fname);
        save(fname,'S');
    end
    
    %% 4 Sec. - 1 Hz. filter
    if ~(strcmp(kcmpnm,"BDF") || strcmp(kcmpnm,"HNZ") || strcmp(kcmpnm,"ENZ"))...
            && (strcmp(kstnm,"CASC") || strcmp(kstnm,"SAGA") || strcmp(kstnm,"ANTI")...
            || strcmp(kstnm,"ANTC") || strcmp(kstnm,"ANTG") || strcmp(kstnm,"ANTM")...
            || strcmp(kstnm,"ANTS") || strcmp(kstnm,"BREF") || strcmp(kstnm,"BNAS")...
            || strcmp(kstnm,"BULB"))
        
        %%
        lfc = 0.25;
        hfc = 1;
        nowTime = dn2dt(now) + hours(5);
        tEnd = dateshift(nowTime,'start','day');
        meanFlag = false; %when false, we take the median amplitude in the window
        rmsFlag = true;
        
        %%
        try
            S = rmsGather(tStart,tEnd,dur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
            S.d = single(S.d);
        catch
            disp('couldnt generate new data');
            continue;
        end
        
        %%
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4Sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
        
        %%
        disp(fname);
        save(fname,'S');
    end
    
    %% 10 Hz.
    if strcmp(kcmpnm,"HHZ") || strcmp(kcmpnm,"BHZ") || strcmp(kcmpnm,"BLZ") || strcmp(kcmpnm,"SHZ")
        %%
        lfc = 10;
        hfc = -inf;
        nowTime = dn2dt(now) + hours(5);
        tEnd = dateshift(nowTime,'start','day');
        meanFlag = false; %when false, we take the median amplitude in the window
        rmsFlag = true;
        
        %%
        try
            S = rmsGather(tStart,tEnd,dur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
            S.d = single(S.d);
        catch
            disp('couldnt generate new data');
            continue;
        end
        
        %%
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_10Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
        
        %%
        disp(fname);
        save(fname,'S');
    end
    
    %% 0.5 - 4 Hz.
    if strcmp(kcmpnm,"HHZ") || strcmp(kcmpnm,"BHZ") || strcmp(kcmpnm,"BLZ")
        
        %%
        lfc = 0.5;
        hfc = 4;
        nowTime = dn2dt(now) + hours(5);
        tEnd = dateshift(nowTime,'start','day');
        meanFlag = false; %when false, we take the median amplitude in the window
        rmsFlag = true;
        
        %%
        try
            S = rmsGather(tStart,tEnd,dur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
            S.d = single(S.d);
        catch
            disp('couldnt generate new data');
            continue;
        end
        
        %%
        fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_2Sec4Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
        
        %%
        disp(fname);
        save(fname,'S');
    end
end

