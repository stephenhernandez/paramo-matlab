%function SangayExplosionsCutRSAM
clear; close all;
cd ~/products/rsam/

secDur = seconds(48*3600);
explosionTimes = [datetime(2020,09,20,09,40,00); ...
    datetime(2021,03,06,03,00,00); ...
    datetime(2021,03,11,09,00,00); ...
    datetime(2021,04,13,00,10,00); ...
    datetime(2021,05,07,14,00,00); ...
    datetime(2022,02,08,08,20,00); ...
    datetime(2022,08,12,15,00,00); ...
    datetime(2022,11,04,09,20,00)];

lTimes = length(explosionTimes);

matFiles = dir('*.mat');
lFiles = length(matFiles);

goodList = ["PUYO";"PORT";"PKYU";"TAIS";"TAMH";"SAGA"];

for i = 1:lTimes
    expTime = explosionTimes(i);
    tStart = expTime - secDur;
    tEnd = expTime + secDur;
    tMain = (tStart:minutes(1):tEnd)';
    n = 0;
    fmtStr = 'yyyy-mm-dd HH:MM:SS';
    dAll = NaN(length(tMain),length(goodList)+1);
    kstnms = [];
    varNames = {};
    for j = 1:lFiles
        fileName = matFiles(j).name;
        kstnm = upper(string(fileName(1:4)));

        if ~ismember(kstnm,goodList)
            continue;
        end

        if strcmp(kstnm,"SAGA")
            if strcmp(fileName(1:5),'saga1')
                kstnm = "SAGA_0.25-2";
            else
                kstnm = "SAGA_2-8";
            end
        end
        kstnms = [kstnms; upper(kstnm)];


        S = load(fileName);
        S = S.S;
        d = S.d;
        t = getTimeVec(S);
        t = t - minutes(1);
        tI = t >= tStart & t <= tEnd;
        tSlice = t(tI);
        dSlice = d(tI);
        nElems = sum(tI);
        if ~nElems
            continue;
        end
        n = n+1;
        isgood = ismember(tMain,tSlice);
        dAll(isgood,n) = dSlice; %(isgood);
        %dAll = [dAll dSlice];
    end
    outFile = sprintf("~/research/now/sangay/ExplosionSangay_%s_v2.txt",datestr(expTime,'yyyy-mmm-dd'));

    %figure();
    figure('units','normalized','outerposition',[0 0 1 1]);
    semilogy(tSlice,dAll,'.'); zoom on; grid on; 
    legend(kstnms,'Location','Best');
    writeFlag = true;
    if ~writeFlag
        continue;
    end
    t = datestr(tSlice,fmtStr);
    %TT = array2timetable(dAll,'RowTimes',t,'VariableNames',cellstr(kstnms));
    %writetimetable(TT,outFile);

    formatSpec = '%s,%f,%f,%f,%f,%f,%f,%f';
    str = compose(formatSpec,t,dAll);
    str = string(str);

    fileID = fopen(outFile,'w');
    fprintf(fileID,'%s\n',str);
    fclose(fileID);    
end

