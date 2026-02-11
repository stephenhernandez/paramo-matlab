function [S,Se] = updateWaveforms(S,MAXDUR,HORFLAG,MINDUR,verboseFlag)
PWD = string(pwd);
if nargin < 2
    MAXDUR = 1600; %8*60;
end

if nargin < 3
    HORFLAG = true;
end

if nargin < 4
    MINDUR = 2;
end

if nargin < 5
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
        "HNZ"; "HNN"; "HNE"; "HN1"; "HN2"; ...
        "ENZ"; "ENN"; "ENE"; "EN1"; "EN2"; ...
        "BDF"];
else
    chans = ["HHZ"; "SHZ"; "BHZ"; "BLZ"; "ENZ"; "HNZ"; "BDF"];
end
chans = sort(chans);

%%
origSNCLs = string(strcat(pull(S,"knetwk"),pull(S,"kstnm"),...
    pull(S,"khole"),pull(S,"kcmpnm")));
startTimes = pull(S,"ref");
endTimes = startTimes + pull(S,"e");

%%
MAXDUR = seconds(MAXDUR);
fudgeFactor = MINDUR;
nowTime = datetime("now") + hours(5);
latestEndTime = nowTime - MAXDUR;
currentDay = dateshift(nowTime,'start','day');

[yyyy,~,~] = datevec(currentDay);
doy = day(currentDay,'doy');
doyStr = padCalendarDay(doy);

%%
if ~(exist('~/rawdata','dir') == 7)
    fprintf(2,'rawdata data directory is not mounted.\n');
    cd(PWD);
    return;
end

%%
DEBUGGING = false;

%%
nUpdated = 0;
Se = populateWaveforms(length(S));

%%
for i = 1:lNets
    cd('~/rawdata/')
    try
        cd(num2str(yyyy));
    catch
        fprintf("Year %d does not exist.\n",yyyy);
        continue;
    end

    %%
    net_ = nets(i);
    try
        cd(net_);
    catch
        fprintf("Network %s does not exist.\n",net_);
        continue;
    end
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
        if strcmp(thisKSTNM,"HNZ.D") || strcmp(thisKSTNM,"HNN.D") ...
                || strcmp(thisKSTNM,"HNE.D") || strcmp(thisKSTNM,"LOG.L")
            continue;
        end
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

                %%
                tic;
                mseedFileName = fullfile(f,mseedname);
                [lia,locb] = ismember(thisSNCL,origSNCLs);

                %%
                nowTime = datetime("now")+hours(5);      %test
                latestEndTime = nowTime - MAXDUR;  %test
                if ~lia
                    S_2 = readMiniSeed(mseedFileName,0,0,latestEndTime);

                    %% too many
                    lS_ = length(S_2);
                    if lS_ > 1
                        continue;
                    end

                    %% empty
                    ref_ = S_2.ref;
                    if isnat(ref_)
                        %fprintf(1,'channel bad, ALMOST added to master list: %s\n',thisSNCL);
                        continue;
                    end

                    %% too short
                    dur_ = S_2.e;
                    if dur_ < minutes(MINDUR)
                        continue;
                    end

                    %%
                    gf_ = S_2.gapFlag;
                    if gf_
                        gapInfo = S_2.gapInfo;
                        old_npts = S_2.npts;
                        newStart = sum(gapInfo(end,:),2);
                        CUTFLAG = newStart < old_npts;
                        if CUTFLAG
                            delta_ = S_2.delta;
                            dtmp = S_2.d;
                            dtmp = dtmp(newStart:end);
                            if ~isfinite(dtmp(1))
                                disp('do something');
                            end

                            %%
                            tmpFs = 1/delta_;
                            if tmpFs > 1
                                tmpFs = round(tmpFs);
                            else
                                delta_ = round(delta_);
                                tmpFs = 1./delta_;
                            end
                            newRef = i2t(newStart,ref_,delta_);
                            S_2 = dealHeader(S_2,dtmp,tmpFs,newRef);
                            dur_ = S_2.e;

                            %% too short
                            if dur_ < minutes(MINDUR)
                                continue;
                            end
                            ref_ = newRef;
                        end
                    end

                    %% too old
                    S2_eTime = ref_ + dur_;
                    if S2_eTime < latestEndTime
                        if verboseFlag
                            fprintf(1,'Channel is too old: %s\n',thisSNCL);
                        end
                        continue;
                    end


                    %% ends in future
                    newNow = datetime("now") + hours(5);
                    if S2_eTime >= newNow
                        if verboseFlag
                            fprintf(1,"Channel ends in FUTURE: %s\n",thisSNCL);
                        end
                        continue;
                    end

                    %% passed all the tests!
                    S = [S; S_2]; %#ok<AGROW>
                    origSNCLs = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
                    endTimes = pull(S,'ref') + pull(S,'e');


                    %%
                    nUpdated = nUpdated + 1;
                    Se(nUpdated,1) = S_2;
                    if verboseFlag
                        fprintf(1,'%d, new channel added to master list: %s\n',nUpdated,mseedname);
                    end
                    continue;
                end

                %%
                thisEnd = endTimes(locb);
                prevEndDay = dateshift(thisEnd,'start','day');
                if prevEndDay == currentDay
                    S_2 = readMiniSeed(mseedFileName,0,0,thisEnd);
                    lS_ = length(S_2);
                    if lS_ > 1
                        continue;
                    end
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
                        if kk == 1
                            S_2_ = readMiniSeed(mseedFileNameTmp,0,0,thisEnd);
                            if isnat(S_2_.ref)
                                continue;
                            end
                        else
                            S_2_ = readMiniSeed(mseedFileNameTmp,0,0);
                        end

                        n = n+1;
                        S_2(n) = S_2_;
                    end

                    %%
                    if n ~= lloops
                        continue;
                    end

                    try
                        S_2 = mergeWaveforms(S_2);
                    catch
                        if verboseFlag
                            fprintf(1,"couldnt merge these two separate days, skipping...\n");
                        end
                        continue;
                    end
                end

                %%
                if isnat(S_2.ref)
                    continue;
                end
                newNow = datetime("now") + hours(5);

                %%
                S1_delta = S_2.delta;
                tmpFs = 1/S1_delta;
                S2_eTime = S_2.ref + S_2.e;

                addedSamples = round(tmpFs*seconds(S2_eTime - endTimes(locb)));
                if addedSamples <= 0 || S2_eTime < latestEndTime - minutes(fudgeFactor) || S2_eTime >= newNow
                    % elapsedTime = toc;
                    if verboseFlag
                        fprintf(1,'no update for %s, elapsed time: %f seconds\n',string(thisSNCL));
                    end
                    continue;
                end

                if S2_eTime >= newNow
                    % elapsedTime = toc;
                    if verboseFlag
                        fprintf(1,'no update for %s, end time greater than now time: %s %s\n',string(thisSNCL),S2_eTime,newNow);
                    end
                    continue;
                end

                %%
                S_1 = S(locb);
                try
                    S_1 = mergeWaveforms([S_1; S_2]);
                catch
                    if ~DEBUGGING
                        if verboseFlag
                            fprintf(1,'-----------------merge broke-----------------\n');
                        end
                        S(locb) = [];
                        origSNCLs(locb) = [];
                        endTimes(locb) = [];
                        continue;
                    end
                    S = S_1;
                    Se = S_2;
                    if verboseFlag
                        fprintf(1,"tried to update: %s,%d\n",thisSNCL,locb);
                        fprintf(2,"ABORT\n");
                    end
                    cd(PWD);
                    return;
                end

                %%
                nUpdated = nUpdated + 1;
                Se(nUpdated,1) = S_2;
                S1_ref = S_1.ref;
                idealStartTime = S2_eTime - MAXDUR;

                gapFlag = S_1.gapFlag;
                if gapFlag
                    %
                    if verboseFlag
                        fprintf('gaps exist for: %s\n',thisSNCL);
                    end
                    gapInfo = S_1.gapInfo;
                    old_npts = S_1.npts;
                    newStart = sum(gapInfo(end,:),2);
                    CUTFLAG = newStart < old_npts; % || S1_ref < idealStartTime;
                    if CUTFLAG
                        dtmp = S_1.d;
                        dtmp = dtmp(newStart:end);
                        if ~isfinite(dtmp(1))
                            disp('do something');
                        end
                        newRef = i2t(newStart,S1_ref,S1_delta);
                        [S_1,new_npts] = dealHeader(S_1,dtmp,tmpFs,newRef);
                        Se(nUpdated,1) = S_1;
                        if verboseFlag
                            fprintf('removed everything except last contiguous section, now %d points\n',new_npts);
                        end
                        S1_ref = newRef;
                    end
                end

                %%
                totSamples = S_1.npts;
                startIndex = t2i(idealStartTime,S1_ref,S1_delta);
                CUTFLAG = startIndex > tmpFs*60*fudgeFactor & startIndex < totSamples;
                if CUTFLAG
                    dtmp = S_1.d;
                    dtmp = dtmp(startIndex:end);
                    newRef = i2t(startIndex,S1_ref,S1_delta);
                    [S_1,totSamples] = dealHeader(S_1,dtmp,tmpFs,newRef);
                    Se(nUpdated,1) = S_1;
                end

                %%
                try
                    S(locb) = S_1;
                catch
                    whos
                    for ii = 1:length(S_1)
                        disp(S_1(ii));
                    end
                    if verboseFlag
                        fprintf(1,'cant assign\n');
                    end
                    cd(PWD);
                    return;
                end
                elapsedTime = toc;

                if CUTFLAG
                    deletedSamples = startIndex - 1;
                    if verboseFlag
                        fprintf(1,'%d, %d, %s updated in: %f; deleted %d; added: %d; total %d\n',nUpdated,locb,thisSNCL,elapsedTime,deletedSamples,addedSamples,totSamples);
                    end
                else
                    if verboseFlag
                        fprintf(1,'%d, %d, %s updated in: %f; added: %d; total %d\n',nUpdated,locb,thisSNCL,elapsedTime,addedSamples,totSamples);
                    end
                end
            end
        end
        cd ..;
    end
    cd ..
end

%%
origSNCLs = string(strcat(pull(S,'knetwk'),pull(S,'kstnm'),pull(S,'khole'),pull(S,'kcmpnm')));
[origSNCLs,sI] = sort(origSNCLs);
S = S(sI);

%% delete channels that are too old or too short
dur = pull(S,'e');
endTimes = pull(S,'ref') + dur;
badI = endTimes < latestEndTime - minutes(fudgeFactor);
if verboseFlag
    fprintf('number of old channels to delete from master vector: %d\n',sum(badI));
end
badI = badI | dur < minutes(MINDUR);
sumbad = sum(badI);
S(badI) = []; %delete bad

if sumbad
    badI = find(badI);
    for i = 1:length(badI)
        if verboseFlag
            fprintf(1,'deleted: %s\n',origSNCLs(badI(i)));
        end
    end
end

%%
if verboseFlag
    fprintf(1,"Nupdated: %d; Deleted: %d\n",nUpdated,sumbad);
end
if ~nUpdated
    Se = populateWaveforms();
    if verboseFlag
        fprintf("<strong>::NOTHING UPDATED::</strong>\n");
    end
    cd(PWD);
    return;
end

%%
Se = Se(1:nUpdated);
updatedSNCLs = string(strcat(pull(Se,"knetwk"),pull(Se,"kstnm"),...
    pull(Se,"khole"),pull(Se,"kcmpnm")));
[updatedSNCLs,sI] = sort(updatedSNCLs);
Se = Se(sI);

dur = pull(Se,'e');
endTimes = pull(Se,'ref') + dur;
badI = endTimes < latestEndTime - minutes(fudgeFactor);
if verboseFlag
    fprintf("number of traces to delete that are too old: %d\n",sum(badI));
end
badI = badI | dur < minutes(MINDUR); %too old or too short?
sumbad = sum(badI);

%%
if ~sumbad
    fprintf("<strong>UPDATED: %d\n</strong>",nUpdated);
    cd(PWD);
    return;
end

if sumbad >= nUpdated
    Se = populateWaveforms();
    fprintf("<strong>strange...::NOTHING UPDATED::</strong>\n");
    cd(PWD);
    return;
end

%%
badI = find(badI);
for i = 1:length(badI)
    if verboseFlag
        fprintf(1,"deleted %s from Supdated vector: Ts: %s; Te: %s\n",...
            updatedSNCLs(badI(i)),Se(badI(i)).ref,Se(badI(i)).ref + Se(badI(i)).e);
    end
end

Se(badI) = [];

%%
fprintf("<strong>UPDATED: %d; DELETED: %d\n</strong>",nUpdated,sumbad);
cd(PWD);