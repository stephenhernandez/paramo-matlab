clear;
close all;
cd ~/Desktop/
lehighFlag = false;
sangayFlag = false;
uclTungurahuaBB = ~true;
uclTungurahuaSP = true;
dataDir = "rawdata_cotopaxi";
if lehighFlag
    snclorig = string(importdata('filelist.txt'));
    goodChans = ["BLZ";"HNZ";"ENZ"];
    goodNets = "EC";
    saveDir = fullfile('~','products','ForLehigh');
elseif sangayFlag
    snclorig = ["EC.PUYO..HHZ";"EC.TAIS..HHZ";"EC.BULB..BHZ";"EC.PORT..HHZ";...
        "EC.TAMH..HHZ";"EC.BPAT..BHZ";"EC.BBIL..BHZ";"EC.BMAS..BHZ";"EC.PKYU..HHZ"];
    goodChans = ["BHZ";"HHZ"];
    goodNets = "EC";
    saveDir = fullfile('~','products','SangayMachineLearningData_v2');
elseif uclTungurahuaBB
    snclorig = ["EC.BULB..BHZ";"EC.BPAT..BHZ";"EC.BRUN..BHZ";"EC.BBIL..BHZ";...
        "EC.BMAS..BHZ";"EC.POND..HHZ";"EC.BRTU..HHZ";"EC.BULB..HHZ";...
        "EC.BPAT..HHZ";"EC.BRUN..HHZ";"EC.BBIL..HHZ";"EC.BMAS..HHZ";];
    goodChans = ["BHZ";"HHZ";"SHZ"];
    goodNets = "EC";
    saveDir = fullfile('~','products','TungurahuaBroadbandUCL');
elseif uclTungurahuaSP
    snclorig = ["EC.RETU..SHZ";"EC.BIL2..SHZ";"EC.ARA2..SHZ";...
        "EC.RUN5..SHZ";"EC.JUI6..SHZ"];
    goodChans = ["BHZ";"HHZ";"SHZ"];
    goodNets = "EC";
    saveDir = fullfile('~','products','TungurahuaShortPeriodUCL');
end

%%
splitstring = split(snclorig,".");
knetwks = splitstring(:,1);
kstnms = splitstring(:,2);
kholes = splitstring(:,3);
kcmpnms = splitstring(:,4);

chanI = ismember(kcmpnms,goodChans) & ismember(knetwks,goodNets);
disp(snclorig(chanI))
mySNCLs = [knetwks(chanI) kstnms(chanI) kholes(chanI) kcmpnms(chanI)];

%%
dayStart = datetime(2005,12,30);
dayEnd = datetime(2016,12,31);
dayInc = 1;
dayVec = (dayStart:dayInc:dayEnd)';
lDays = length(dayVec);

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

for i = 1:lDays
    day_ = dayVec(i);
    cd(saveDir);
    [yyyy,mm,dd] = datevec(day_);
    yyyyStr = num2str(yyyy);
    doy = day(datetime(yyyy,mm,dd),'doy');
    doyStr = padCalendarDay(doy);
    if ~exist(yyyyStr,"dir")
        mkdir(yyyyStr);
    end
    cd(yyyyStr);
    for j = 1:sum(chanI)
        thisSNCL = mySNCLs(j,:);
        thisknetwk = thisSNCL(1);
        thiskstnm = thisSNCL(2);
        thiskhole = thisSNCL(3);
        thiskcmpnm = thisSNCL(4);
        otherCMPNMs = char(thiskcmpnm);
        baseChan = otherCMPNMs(1:2);
        otherCMPNMs = [strcat(baseChan,"E"); strcat(baseChan,"N"); strcat(baseChan,"Z")];
        if ~exist(thisknetwk,"dir")
            mkdir(thisknetwk);
        end
        cd(thisknetwk);

        if ~exist(thiskstnm,"dir")
            mkdir(thiskstnm);
        end
        cd(thiskstnm);

        for k = 1:length(otherCMPNMs)
            thiskcmpnm = otherCMPNMs(k);
            chanDir = strcat(thiskcmpnm,".D");
            if ~exist(chanDir,"dir")
                mkdir(chanDir);
            end
            cd(chanDir);

            rootDir = fullfile('~',dataDir,yyyyStr,thisknetwk,thiskstnm,...
                chanDir);
            mseedfilename = sprintf("%s.%s.%s.%s.D.%s.%s",...
                thisknetwk,thiskstnm,thiskhole,thiskcmpnm,yyyyStr,doyStr);
            fullfilename = fullfile(rootDir,mseedfilename);

            cmd = sprintf("cp %s .;",fullfilename);
            %disp(cmd);
            %try
            status = unix(cmd);
            %catch
            %    fprintf(2,"could not run cmd: %s",cmd);
            %end
            cd ../
        end
        cd ../..;
    end
    cd ..;
end
