clear; close all;

writeDirBase = fullfile('~','products','events','waveforms');
saveDirectory = fullfile('~','products','events','html');
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

cd(saveDirectory)
totDurMinutes = 6;

firstAttempt = false;
if ~firstAttempt
    T = readtable('eq.txt');
    strID = T.Var18;
    for i = 1:length(strID)
        newStr(i,1) = string(strrep(strID(i), char(39), ''));
    end
    ids = newStr;
    clear newStr strID;
    searchComponents = ["BHZ";"HHZ";"SHZ"]; %;"HNZ";"BLZ";"ENZ"
else
    load tmp.mat
    searchComponents = ["BHZ";"HHZ";"SHZ";"EHZ";"HNZ";"BLZ";"ENZ"];
end

for i = 1:length(ids)
    E(i,1) = readSCBulletin(ids(i));
end

t = pull(E,'t');
eqmag = pull(E,'mag');
[eqmag,sI] = sort(eqmag,'descend');
ids = ids(sI);
t = t(sI);

%%
close all;
HORFLAG = true;
verboseFlag = true;
for i = 1:length(t)
    eventID = ids(i);
    origt_ = t(i);
    eqmag_ = eqmag(i);

    [yyyy,mm,~] = datevec(origt_);
    yyyyStr = string(num2str(yyyy));
    if mm < 10
        mmStr = string(strcat('0',num2str(mm)));
    else
        mmStr = string(num2str(mm));
    end

    eqmag_ = round(eqmag_*100)/100;
    dirName = sprintf("%s_M%3.2f_%s",eventID,eqmag_,datestr(origt_,30));

    writeDir = fullfile(writeDirBase,yyyyStr,mmStr,dirName);
    if ~exist(writeDir,'dir')
        mkdir(writeDir);
    end

    tStart = origt_-minutes(1);
    tEnd = tStart + minutes(totDurMinutes);
    kstnmListFlag = ~true;
    if kstnmListFlag
        kstnmList = ["RVRD";"PTGL";"ESM1";"SNLR";"APED";"PAC1";"JAMA";...
            "FLFR";"OTAV";"PDNS";"MAG1";"BV15";"FLF1";"FLFR";"CUSE";...
            "LGCB"];
        Sall = extractWaveforms(tStart,tEnd,kstnmList,["HHZ";"HHN";"HHE";...
            "BHZ";"BHN";"BHE";"SHZ";"SHN";"SHE";"HNZ";"HNN";"HNE";...
            "ENZ";"ENN";"ENE";"EHZ";"EHN";"EHE"],...
            ["EC";"IU"],["";"01";"10"],true,true,'~/rawdata',1,false);
        Sall = Sall(:);
        refs = pull(Sall,'ref');
        badI = isnat(refs);
        Sall(badI) = [];
    else
        Sall = seedData2(dateshift(origt_,'start','day'),tStart,tEnd,...
            HORFLAG,verboseFlag,searchComponents);
    end

    badI = isnat(pull(Sall,'ref'));
    Sall(badI) = [];
    if isempty(Sall)
        continue;
    end

    %%
    for j = 1:length(Sall)
        S1 = Sall(j);
        knetwks_ = S1.knetwk;
        kstnms_ = S1.kstnm;
        kholes_ = S1.khole;
        kcmpnm_ = S1.kcmpnm;
        cutStart = S1.ref;

        fNameSac_ = sprintf("%s.%s.%s.%s_%s.SAC",knetwks_,kstnms_,kholes_,kcmpnm_,datestr(cutStart,30));
        fNameSAC = fullfile(writeDir,fNameSac_);

        try
            disp(fNameSac_);
            sacwrite(fNameSAC,S1);
        catch
            fprintf(2,'couldnt write sac: %s\n',fNameSac);
        end
    end
    clear Sall;
end
