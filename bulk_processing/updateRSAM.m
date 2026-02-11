function updateRSAM(kMaster,currentDay,verboseFlag)
if nargin < 1
    disp('doing all');
end

if nargin < 2
    currentDay = dateshift(dn2dt(now)+hours(5),'start','day'); %add 5 hours for UTC
end

if nargin < 3
    verboseFlag = false;
end

%%
nets = ["EC";"IU";"OP";"PE";"CM"];
lNets = length(nets);

chans = ["HHZ"; "SHZ"; "BHZ"; "HNZ"; "BLZ"; "ENZ"; "BDF"];
lChans = length(chans);

[yyyy,~,~] = datevec(currentDay);
doy = day(currentDay,'doy');

if doy < 10
    doyStr = ['00',num2str(doy)];
elseif doy < 100
    doyStr = ['0',num2str(doy)];
else
    doyStr = num2str(doy);
end



%%
secDur = 30;
if ~(exist('~/rawdata','dir') == 7)
    fprintf(2,'rawdata data directory is not mounted.\n');
    return;
end

cd('~/rawdata/')
cd(num2str(yyyy));
for i = 1:lNets
    net_ = nets(i);
    dirs = dir(nets(i));
    kstnmNames = string({dirs.name})';
    goodI = ([dirs.isdir])' & ~strcmp(kstnmNames,".") & ~strcmp(kstnmNames,"..");
    kstnmNames = kstnmNames(goodI);
    lNames = sum(goodI);
    if ~lNames
        if verboseFlag
            fprintf('there are no valid SNCLs for this network\n');
        end
        continue;
    end

    %%
    cd(net_);
    for j = 1:lNames
        for k = 1:lChans
            chanStr = strcat(chans(k),'.D');
            f = fullfile(pwd,kstnmNames(j),chanStr);
            if exist('kMaster','var')
                miniseedfiles = dir(strcat(f,'/',net_,'.',char(kMaster),'.*.',doyStr));
            else
                miniseedfiles = dir(strcat(f,'/',net_,'.*.',doyStr));
            end

            %%
            if isempty(miniseedfiles)
                continue;
            end

            %%
            lminiseedfile = length(miniseedfiles);
            for l = 1:lminiseedfile
                if verboseFlag
                    fprintf('processing: %s\n',datestr(currentDay));
                end
                mseedname = miniseedfiles(l).name;
                splits = split(mseedname,".");
                knetwk = string(splits(1));
                kstnm = string(splits(2));
                khole = string(splits(3));
                kcmpnm = string(splits(4));

                %% select only "01" BDF channel
                if strcmp(kcmpnm,"BDF") && (strcmp(khole,"02") || strcmp(khole,"03"))
                    continue;
                end

                %%
                if strcmp(kcmpnm,'BDF')
                    lfc = 1/4;
                    hfc = 1;
                else
                    lfc = 1;
                    hfc = 8;
                end

                %%
                if ~strcmp(kcmpnm,"BDF")
                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                else
                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                end

                %%
                nowTime = dn2dt(now) + hours(5);
                tStart = datetime(1997,01,01);
                tEnd = dateshift(nowTime,'start','day');
                meanFlag = false;   % when false, we take the median amplitude in the window instead of the mean
                rmsFlag = true;     % root mean square (as opposed to the mean of the absolute value)

                if ~exist(fname,'file')
                    try
                        S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                        %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        if verboseFlag
                            fprintf(1,'saving file %s\n',fname);
                        end
                        save(fname,'S');
                        clear S;
                    catch
                        warning('couldnt generate new file %s\n',fname);
                        continue;
                    end
                end

                %% update already existing
                try
                    if verboseFlag
                        fprintf(1,"loading %s\n",fname);
                    end
                    load(fname,'S');
                catch
                    warning('Couldnt load %s\n',fname);
                    continue;
                end

                %%
                lS = length(S);
                d = S(1).d;
                if isnat(S(1).ref) || isempty(d) || lS > 1
                    S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                    %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                    save(fname,'S');
                    clear S;
                    continue;
                end

                %%
                newStart = dateshift(S.ref + S.e,'start','day');
                if newStart < datetime(2018,01,01)
                    if verboseFlag
                        fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                    end
                    continue;
                end

                %%
                tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                %
                try
                    [M,status] = mergeWaveforms([S; S_]);
                catch
                    S_ = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S_,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                    [M,status] = mergeWaveforms([S; S_]);
                end
                if status
                    S = M(1);
                    save(fname,'S');
                end
                clear S_ M S;

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
                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                    if ~exist(fname,'file')
                        try
                            S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                            if verboseFlag
                                fprintf(1,'saving file %s\n',fname);
                            end
                            save(fname,'S');
                            clear S;
                        catch
                            warning('couldnt generate new file %s\n',fname);
                        end
                        continue;
                    end

                    %% update already existing
                    try
                        if verboseFlag
                            fprintf(1,"loading %s\n",fname);
                        end
                        load(fname,'S');
                    catch
                        warning('Couldnt load %s',fname);
                        continue;
                    end

                    %%
                    lS = length(S);
                    d = S(1).d;
                    if isnat(S(1).ref) || isempty(d) || lS > 1
                        S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                        %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        save(fname,'S');
                        clear S;
                        continue;
                    end

                    %%
                    newStart = dateshift(S.ref + S.e,'start','day');
                    if newStart < datetime(2018,01,01)
                        if verboseFlag
                            fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                        end
                        continue;
                    end

                    %%
                    tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                    S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                    try
                        [M,status] = mergeWaveforms([S; S_]);
                    catch
                        S_ = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S_,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        [M,status] = mergeWaveforms([S; S_]);
                    end
                    if status
                        S = M(1);
                        save(fname,'S');
                    end
                    clear S_ M S;
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

                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4Sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                    if ~exist(fname,'file')
                        try
                            S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                            if verboseFlag
                                fprintf(1,'saving file %s\n',fname);
                            end
                            save(fname,'S');
                            clear S;
                        catch
                            warning('couldnt generate new file %s\n',fname);
                        end
                        continue;
                    end

                    %% update already existing
                    try
                        if verboseFlag
                            fprintf(1,"loading %s\n",fname);
                        end
                        load(fname,'S');
                    catch
                        warning('Couldnt load %s',fname);
                        continue;
                    end

                    %%
                    lS = length(S);
                    d = S(1).d;
                    if isnat(S(1).ref) || isempty(d) || lS > 1
                        S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                        %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        save(fname,'S');
                        clear S;
                        continue;
                    end

                    newStart = dateshift(S.ref + S.e,'start','day');
                    if newStart < datetime(2018,01,01)
                        if verboseFlag
                            fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                        end
                        continue;
                    end

                    tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                    S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                    try
                        [M,status] = mergeWaveforms([S; S_]);
                    catch
                        S_ = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S_,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        [M,status] = mergeWaveforms([S; S_]);
                    end
                    if status
                        S = M(1);
                        save(fname,'S');
                    end
                    clear S_ M S;
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

                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_10Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                    if ~exist(fname,'file')
                        try
                            S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                            if verboseFlag
                                fprintf(1,'saving file %s\n',fname);
                            end
                            save(fname,'S');
                            clear S;
                        catch
                            warning('couldnt generate new file %s\n',fname);
                        end
                        continue;
                    end

                    %% update already existing
                    try
                        if verboseFlag
                            fprintf(1,"loading %s\n",fname);
                        end
                        load(fname,'S');
                    catch
                        warning('Couldnt load %s',fname);
                        continue;
                    end

                    %%
                    lS = length(S);
                    d = S(1).d;
                    if isnat(S(1).ref) || isempty(d) || lS > 1
                        S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                        %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        save(fname,'S');
                        clear S;
                        continue;
                    end

                    newStart = dateshift(S.ref + S.e,'start','day');
                    if newStart < datetime(2018,01,01)
                        if verboseFlag
                            fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                        end
                        continue;
                    end

                    %%
                    tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                    S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                    try
                        [M,status] = mergeWaveforms([S; S_]);
                    catch
                        S_ = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S_,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        [M,status] = mergeWaveforms([S; S_]);
                    end
                    if status
                        S = M(1);
                        save(fname,'S');
                    end
                    clear S_ M S;
                end

                %% 0.5 - 4 Hz.
                if strcmp(kcmpnm,"HHZ") || strcmp(kcmpnm,"BHZ") || strcmp(kcmpnm,"BLZ")
                    lfc = 0.5;
                    hfc = 4;
                    nowTime = dn2dt(now) + hours(5);
                    tEnd = dateshift(nowTime,'start','day');
                    meanFlag = false; %when false, we take the median amplitude in the window
                    rmsFlag = true;

                    fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_2Sec4Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                    if ~exist(fname,'file')
                        try
                            S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                            if verboseFlag
                                fprintf(1,'saving file %s\n',fname);
                            end
                            save(fname,'S');
                            clear S;
                        catch
                            warning('couldnt generate new file %s\n',fname);
                            continue;
                        end
                    end

                    %% update already existing
                    try
                        if verboseFlag
                            fprintf(1,"loading %s\n",fname);
                        end
                        load(fname,'S');
                    catch
                        warning('Couldnt load %s',fname);
                        continue;
                    end

                    %%
                    lS = length(S);
                    d = S(1).d;
                    if isnat(S(1).ref) || isempty(d) || lS > 1
                        S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                        %S = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        save(fname,'S');
                        clear S;
                        continue;
                    end

                    newStart = dateshift(S.ref + S.e,'start','day');
                    if newStart < datetime(2018,01,01)
                        if verboseFlag
                            fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                        end
                        continue;
                    end

                    %%
                    tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                    S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                    try
                        [M,status] = mergeWaveforms([S; S_]);
                    catch
                        S_ = rmfield(rmfield(rmfield(rmfield(rmfield(rmfield(S_,'poles'),'zeroes'),'A0'),'sensitivity'),'gain'),'constant');
                        [M,status] = mergeWaveforms([S; S_]);
                    end
                    if status
                        S = M(1);
                        save(fname,'S');
                    end
                    clear S_ M S;
                end
            end
        end
    end
    cd ..;
end
