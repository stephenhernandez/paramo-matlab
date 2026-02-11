function S = seedData(MAXDUR,HORFLAG,verboseFlag)
if nargin < 1
    MAXDUR = 1600;
end

if nargin < 2
    HORFLAG = false;
end

if nargin < 3
    verboseFlag = 0;
end

%%
nets = ["EC";"IU";"OP";"PE";"CM"];
nets = sort(nets);
lNets = length(nets);

%%
if HORFLAG
    chans = ["HHZ"; "HHN"; "HHE"; "HH1"; "HH2"; ...
        "SHZ"; "SHN"; "SHE"; "SH1"; "SH2"; ...
        "BHZ"; "BHN"; "BHE"; "BH1"; "BH2"; ...
        "BLZ"; "BLN"; "BLE"; "BL1"; "BL2"; ...
        "HNZ"; "HNN"; "HNE"; ...
        "ENZ"; "ENN"; "ENE"; ...
        "EHZ"; "EHN"; "EHE"; ...
        "BDF"; "HDF"];
else
    chans = ["HHZ";"SHZ";"BHZ";"BLZ";"BDF";"HDF";"ENZ";"EHZ";"HNZ"];
end
chans = sort(chans);

%%
nowTime = datetime("now")+hours(5);
currentDay = dateshift(nowTime,'start','day'); %add 5 hours for UTC

%%
[yyyy,~,~] = datevec(currentDay);
doy = day(currentDay,'doy');
doyStr = padCalendarDay(doy);

%%
rawDataDir = '~/rawdata';
existFlag = exist(rawDataDir,'dir');
if existFlag ~= 7
    fprintf(2,'%s: data directory is not mounted.\n',rawDataDir);
    return;
end

%%
maxN = 400;
S = populateWaveforms(maxN);

%%
n = 0;
yyyyRootDir = fullfile(rawDataDir,num2str(yyyy));
for i = 1:lNets
    net_ = nets(i);
    netRootDir = fullfile(yyyyRootDir,net_);
    if ~exist(netRootDir,'dir')
        continue
    end

    dirs = dir(netRootDir);
    ldirs = length(dirs);
    if ~ldirs
        continue;
    end

    %%
    kstnmNames = pull(dirs,'name',"");
    goodI = logical(pull(dirs,'isdir')) & ~strcmp(kstnmNames,".") & ~strcmp(kstnmNames,"..");
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

        kstnmRootDir = fullfile(netRootDir,thisKSTNM);
        if ~exist(kstnmRootDir,'dir')
            if verboseFlag
                fprintf("something wrong with: %s\n",kstnmRootDir);
            end
            continue
        end

        dirs = dir(kstnmRootDir);
        ldirs = length(dirs);
        if ~ldirs
            fprintf('%s is empty, continuing\n',kstnmRootDir)
            continue;
        end

        chanDirNames = pull(dirs,'name',"");
        goodI = logical(pull(dirs,'isdir')) & ~strcmp(chanDirNames,".") & ~strcmp(chanDirNames,"..");
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
            chanRootDir = fullfile(kstnmRootDir,thisChanDirD);
            fileStr = strcat(net_,'.',thisKSTNM,'.*.',doyStr);
            searchStr = fullfile(chanRootDir,fileStr);

            miniseedfile = dir(searchStr);
            if isempty(miniseedfile)
                continue;
            end

            %%
            lminiseedfile = length(miniseedfile);
            if ~lminiseedfile
                continue;
            end

            for l = 1:lminiseedfile
                if verboseFlag
                    disp(currentDay);
                end
                mseedname = miniseedfile(l).name;

                %%
                mseedFileName = fullfile(chanRootDir,mseedname);
                tic;
                S_ = readMiniSeed(mseedFileName,0,0,nowTime - seconds(MAXDUR));
                if isnat(S_.ref)
                    fprintf(2,'something wrong with: %s\n',string(mseedFileName));
                    continue;
                end

                %%
                newNow = datetime("now") + hours(5);
                eTime = S_.ref + S_.e;
                if eTime >= newNow
                    fprintf(1,'%s is bad, has end time in FUTURE!\n',mseedname);
                    continue;
                end
                
                n = n + 1;
                S(n) = S_;
                elapsedTime = toc;
                fprintf(1,'%d, %s: %f\n',n,string(mseedFileName),elapsedTime);
            end
        end
    end
end
S = S(1:n);
