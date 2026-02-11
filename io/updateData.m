function [S,Se] = updateData(S,MAXDUR,HORFLAG,verboseFlag)
PWD = string(pwd);
if nargin < 2
    MAXDUR = 8;
end

if nargin < 3
    HORFLAG = true;
end

if nargin < 4
    verboseFlag = false;
end

%%
nets = ["EC";"IU";"OP";"PE";"CM"];
nets = sort(nets);
lNets = length(nets);

if HORFLAG
    chans = ["HHZ"; "HHN"; "HHE"; "HH1"; "HH2"; ...
        "SHZ"; "SHN"; "SHE"; "SH1"; "SH2"; ...
        "BHZ"; "BHN"; "BHE"; "BH1"; "BH2"; ...
        "BLZ"; "BLN"; "BLE"; "BL1"; "BL2"; ...
        "BDF"];
else
    chans = ["HHZ"; "SHZ"; "BHZ"; "BLZ"; "BDF"];
end
chans = sort(chans);
%lChans = length(chans);

%%
origSNCLs = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
startTimes = pull(S,'ref');
endTimes = startTimes + pull(S,'e');

%%
MAXDUR = hours(MAXDUR);
fudgeFactor = 4;
nowTime = dn2dt(now)+hours(5);
idealStartTime = nowTime - MAXDUR;
currentDay = dateshift(nowTime,'start','day');

[yyyy,~,~] = datevec(currentDay);
doy = day(currentDay,'doy');
doyStr = padCalendarDay(doy);

%%
if ~(exist('~/rawdata','dir') == 7)
    fprintf(2,'rawdata data directory is not mounted.\n');
    return;
end

%%
DEBUGGING = false;

%%
nUpdated = 0;
Se = S;
cd('~/rawdata/')
cd(num2str(yyyy));
for i = 1:lNets
    net_ = nets(i);
    cd(net_);
    dirs = dir();
    
    %%
    kstnmNames = string({dirs.name})';
    goodI = ([dirs.isdir])' & ~strcmp(kstnmNames,".") & ~strcmp(kstnmNames,"..");
    kstnmNames = kstnmNames(goodI);
    lkstnmNames = sum(goodI);
    if ~lkstnmNames
        if verboseFlag
            fprintf('there are no valid SNCLs for this network\n');
        end
        continue;
    end
    
    %%
    for j = 1:lkstnmNames
        thisKSTNM = kstnmNames(j);
        cd(thisKSTNM);
        dirs = dir();
        
        %%
        chanDirNames = string({dirs.name})';
        goodI = ([dirs.isdir])' & ~strcmp(chanDirNames,".") & ~strcmp(chanDirNames,"..");
        chanDirNames = chanDirNames(goodI);
        lChanDirNames = sum(goodI);
        if ~lChanDirNames
            if verboseFlag
                fprintf('there are no valid SNCLs for this network\n');
            end
            continue;
        end
        
        %%
        for k = 1:lChanDirNames
            thisChanDirD = chanDirNames(k);
            thisChan = strsplit(thisChanDirD,".D");
            thisChan = thisChan(1);
            if ~ismember(thisChan,chans)
                continue;
            end
            
            %%
            f = fullfile(pwd,thisChanDirD);
            fileStr = strcat(net_,'.',thisKSTNM,'.*.',doyStr);
            searchStr = fullfile(f,fileStr);
            
            %%
            miniseedfile = dir(searchStr);
            if isempty(miniseedfile)
                continue;
            end
            
            %%
            lminiseedfile = length(miniseedfile);
            for l = 1:lminiseedfile
                if verboseFlag
                    disp(currentDay);
                end
                mseedname = miniseedfile(l).name;
                splits = split(mseedname,".");
                knetwk = splits(1);
                kstnm = splits(2);
                khole = splits(3);
                kcmpnm = splits(4);
                
                %%
                thisSNCL = string(strcat(knetwk,kstnm,khole,kcmpnm));
                
                %% select only "01" BDF channel
                if strcmp(kcmpnm,"BDF") && (strcmp(khole,"02") || strcmp(khole,"03"))
                    continue;
                end
                
                %%
                tic;
                mseedFileName = fullfile(f,mseedname);
                [lia,locb] = ismember(thisSNCL,origSNCLs);
                
                if ~lia
                    S_2 = readMiniSeed(mseedFileName,0,0,idealStartTime);
                    if isnat(S_2.ref)
                        %fprintf(1,'channel ALMOST added to master list: %s\n',thisSNCL);
                        continue;
                    end
                    
                    %%
                    S2_eTime = S_2.ref + S_2.e;
                    if S2_eTime < idealStartTime
                        continue;
                    end
                    S = [S; S_2]; %#ok<AGROW>
                    origSNCLs = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
                    endTimes = pull(S,'ref') + pull(S,'e');
                    fprintf(1,'new channel added to master list: %s\n',mseedname);
                    continue;
                end
                
                %%
                thisEnd = endTimes(locb);
                prevEndDay = dateshift(thisEnd,'start','day');
                if prevEndDay == currentDay
                    S_2 = readMiniSeed(mseedFileName,0,0,thisEnd);
                else
                    loopDays = (prevEndDay:currentDay)';
                    lloops = length(loopDays);
                    S_2 = populateWaveforms(lloops);
                    
                    %%
                    n = 0;
                    for kk = 1:lloops
                        loopDay = loopDays(kk);
                        doy2 = day(loopDay,'doy');
                        doy2str = padCalendarDay(doy2);
                        mseedFileNameTmp = char(mseedFileName);
                        mseedFileNameTmp(end-2:end) = [];
                        mseedFileNameTmp = [mseedFileNameTmp doy2str]; %#ok<AGROW>
                        
                        %%
                        S_2_ = readMiniSeed(mseedFileNameTmp,0,0,thisEnd);
                        if isnat(S_2_.ref)
                            break;
                        end
                        n = n+1;
                        S_2(n) = S_2_;
                        
                    end
                    
                    %%
                    if n ~= lloops
                        continue;
                    end
                    S_2 = mergeWaveforms(S_2);
                end
                
                %%
                if isnat(S_2.ref)
                    continue;
                end
                
                %%
                S1_delta = S_2.delta;
                tmpFs = 1/S1_delta;
                S2_eTime = S_2.ref + S_2.e;
                newSamples = round(tmpFs*seconds(S2_eTime - endTimes(locb)));
                if newSamples <= 0 || S2_eTime < idealStartTime - minutes(fudgeFactor)
                    % elapsedTime = toc;
                    % fprintf(1,'no update for %s, elapsed time: %f seconds\n',string(thisSNCL),elapsedTime);
                    continue;
                end
                
                %%
                S_1 = S(locb);
                try
                    S_1 = mergeWaveforms([S_1; S_2]);
                catch
                    if ~DEBUGGING
                        continue;
                    end
                    S = S_1;
                    Se = S_2;
                    fprintf(1,'tried to update: %s,%d\n',thisSNCL,locb);
                    fprintf(2,'ABORT\n');
                    return;
                end
                
                %%
                nUpdated = nUpdated + 1;
                Se(nUpdated,1) = S_2;
                S1_ref = S_1.ref;
                gapFlag = S_1.gapFlag;
                if gapFlag
                    %
                    fprintf('gaps exist for: %s\n',thisSNCL);
                    gapInfo = S_1.gapInfo;
                    old_npts = S_1.npts;
                    newStart = sum(gapInfo(end,:),2);
                    doCut = newStart < old_npts;
                    if doCut
                        dtmp = S_1.d;
                        dtmp = dtmp(newStart:end);
                        if ~isfinite(dtmp(1))
                            % do something
                        end
                        newRef = i2t(newStart,S1_ref,S1_delta);
                        [S_1,new_npts] = dealHeader(S_1,dtmp,tmpFs,newRef);
                        fprintf('i cut everything but last continguous section, there are now %d points\n',old_npts-newStart+1);
                        S1_ref = newRef;
                        newSamples = new_npts;
                    end
                end
                
                %%
                
                startIndex = t2i(idealStartTime,S1_ref,S1_delta);
                doCut = startIndex > tmpFs*60*fudgeFactor;
                if doCut
                    dtmp = S_1.d;
                    dtmp = dtmp(startIndex:end);
                    if ~isfinite(dtmp(1))
                        % do something
                    end
                    newRef = i2t(startIndex,S1_ref,S1_delta);
                    S_1 = dealHeader(S_1,dtmp,tmpFs,newRef);
                end
                S(locb) = S_1;
                elapsedTime = toc;
                
                if doCut
                    deletedSamples = startIndex - 1;
                    fprintf(1,'%d, %d, updated %s in %f sec., deleted %d samples, added %d samples\n',nUpdated,locb,thisSNCL,elapsedTime,deletedSamples,newSamples);
                else
                    fprintf(1,'%d, %d, updated %s in %f sec., added %d samples\n',nUpdated,locb,thisSNCL,elapsedTime,newSamples);
                end
            end
        end
        cd ..;
    end
    cd ..
end

%% delete channels that are too old
endTimes = pull(S,'ref') + pull(S,'e');
badI = endTimes < idealStartTime - minutes(fudgeFactor);
S(badI) = [];

%%
origSNCLs = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
[~,sI] = sort(origSNCLs);
S = S(sI);
Se = Se(1:nUpdated);

%%
cd(PWD);
