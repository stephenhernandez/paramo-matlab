function fileList = existsToday(queryDay,HORFLAG,verboseFlag,chans,nets,locList)
cd ~
if nargin < 1
    queryDay = datetime('now')+hours(5);
end

if nargin < 2
    HORFLAG = false;
end

if nargin < 3
    verboseFlag = true;
end

if nargin < 4
    chans = [];
end

if nargin < 5
    nets = [];
end

if nargin < 6
    locList = [];
end

%%
if isempty(chans)
    chans = ["HHZ";"BHZ";"BDF";"HDF";"ENZ";"BLZ";"HNZ";"SHZ";"EHZ"];
end

if isempty(nets)
    nets = ["EC";"IU";"OP";"PE";"CM";"4B";"8G";"XE";"US";"G";"II";"IU";...
        "NA";"Y2";"XF";"XX";"PW";"VE"]; %["9D";"C1";"CX"];
end

if isempty(locList)
    locList = ["";"00";"01";"02";"03";"04";"05";"10";"20";"32"];
end

%%
nets = sort(nets);

%%
if HORFLAG
    horTypes = ["E";"N";"Z";"1";"2"];
    chans_ = [];
    for i = 1:length(chans)
        chan1 = chans(i);
        band_gain = char(chan1);
        band_gain = string(band_gain(1:2));
        chanHor = [];
        for j = 1:length(horTypes)
            chanHor = [chanHor; strcat(band_gain,horTypes(j))];
        end
        chans_ = [chans_; chanHor];
    end
    chans = unique([chans; chans_]);
end
chans = sort(chans);

%%
currentDay = dateshift(queryDay,'start','day'); %add 5 hours for UTC
maxN = 400;

%%
[yyyy,~,~] = datevec(currentDay);
yyyyStr = num2str(yyyy); %string(yyyy);
doy = day(currentDay,'doy');
doyStr = sprintf("%03d",doy);

%%
dirList = "~/rawdata_cotopaxi"; %["~/rawdata/";"~/rawdata_cotopaxi/";"~/rawdata_old/"];
lDirs = length(dirList);
fileList = [];
for i = 1:lDirs
    dataDir = dirList(i); %fullfile('~','rawdata');
    if ~(exist(dataDir,'dir') == 7)
        fprintf(2,'rawdata data directory is not mounted.\n');
        return;
    end
    %try
        fileList_ = existsToday_(dataDir,yyyyStr,doyStr,verboseFlag,maxN,chans,nets,locList,currentDay);
    %catch
    %    fprintf("soemthing went wrong\n");
    %end
    fileList = [fileList; fileList_];
end
fileList = unique(fileList);
end

%%
function fileList = existsToday_(dataDir,yyyyStr,doyStr,verboseFlag,maxN,chans,nets,locList,currentDay)
fileList = repmat("",maxN,1);
lNets = length(nets);
if verboseFlag
    disp(currentDay);
end
nn = 0;
cd(dataDir);
cd(yyyyStr);
for i = 1:lNets
    net_ = nets(i);
    if ~exist(net_,'dir')
        fprintf('no stations within %s network\n',net_);
        continue;
    end
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
                cd ..;
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
            %EC.BTER..HHZ.D.2022.044
            for m = 1:length(locList)
                hole_ = locList(m);
                fileStr = strcat(net_,'.',thisKSTNM,'.',hole_,'.',thisChan,'.D.',yyyyStr,'.',doyStr);
                searchStr = fullfile(f,fileStr);

                %%
                miniseedfile = dir(searchStr);
                if isempty(miniseedfile)
                    continue;
                end

                %%
                lminiseedfile = length(miniseedfile);
                for l = 1:lminiseedfile
                    mseedname = miniseedfile(l).name;
                    splits = split(mseedname,".");
                    khole = (splits(3));
                    kcmpnm = (splits(4));

                    %%
                    mySNCL_ = strcat(net_,".",thisKSTNM,".",khole,".",kcmpnm);
                    if verboseFlag
                        disp(mySNCL_)
                    end
                    nn = nn + 1;
                    fileList(nn) = mySNCL_;
                end
            end
        end
        cd ../;
    end
    cd ../;
end
fileList = fileList(1:nn);
end