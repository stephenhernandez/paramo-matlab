function guralp_rename(stnm,netwk,hole,stnumber,baseDir)
%cd ~/data/galapagos/service_runs/SERVICE3/logs/;
%files = dir('*.log');
% clear;
% close all;
% baseDir = '~/data/galapagos/service_runs/SN12_DFD_FEB2022';
% stnm = "SN12";
% netwk = "EC";
% hole = "";
% stnumber = "57";
% another example: guralp_rename("CBHM","EC","","57","~/data/galapagos/service_runs/November2025/CBHM_Nov2025/");

%%
cd(baseDir)
dirs = dir([char(stnumber),'*']);
%dirs = flipud(dirs);
ldirs = length(dirs);
for i = 1%ldirs
    currentDirName = dirs(i).name;
    cd(currentDirName);
    files = dir('*.msd');
    lFiles = length(files);
    if ~lFiles
        files = dir('*.txt');
        lFiles = length(files);
        if ~lFiles
            fprintf(['neither miniseed nor text files found, ...' ...
                'trying another directory\n']);
            cd(baseDir);
            continue;
        end
    end

    for j = 1:lFiles
        origFilename = files(j).name;
        filename = string(split(origFilename,"."));
        filenamesuffix = filename(2);

        filename = filename(1);
        origSplits = split(filename,"_");
        date_time = origSplits(1);

        date_time = char(date_time);
        yyyy = str2double(date_time(1:4));
        mm = str2double(date_time(5:6));
        dd = str2double(date_time(7:8));
        nowDay = datetime(yyyy,mm,dd);
        if strcmp(filenamesuffix,"msd")
            S = readMiniSeed(origFilename);
            cmpnm = S.kcmpnm;
            ref = S.ref;
            if isnat(ref) || strcmp(cmpnm,"")
                fprintf('abort\n');
                break;
            end
            S.kstnm = string(stnm);
            if strcmp(S.knetwk,"")
                S.knetwk = netwk;
            end
            if ~strcmp(S.khole,hole)
                hole = S.khole;
            end
        else
            cmpnm = "LOG";
        end

        %%
        newFileName = sprintf('%s.%s.%s.%s.D.%d.%03d',netwk,stnm,hole,cmpnm,...
            yyyy,day(datetime(yyyy,mm,dd),'doy'));
        fprintf("%s %s\n",currentDirName,newFileName);

        newDir = sprintf('%d/%s/%s/%s.D',yyyy,netwk,stnm,cmpnm);
        if ~exist(newDir,'dir')
            mkdir(newDir);
        end

        %%
        if exist("S",'var')
            if j > 1
                if nowDay == oldDay
                    try
                        S = mergeWaveforms([Sold; S]);
                    catch
                        continue;
                    end
                end
            end

            mkmseed(sprintf('%s.%s.%s.%s.D.%d.%03d',...
                S.knetwk,S.kstnm,S.khole,S.kcmpnm,...
                year(S.ref),day(S.ref,'doy')),S.d,datenum(S.ref),1./S.delta,11,4096);

            oldDay = nowDay;
            Sold = S;

            cmd = sprintf('\\mv -f %s %d/%s/%s/%s.D',newFileName,yyyy,netwk,stnm,cmpnm);
            unix(cmd);
            disp(cmd);
        else
            cmd = sprintf('\\cp -f %s %d/%s/%s/%s.D/%s',origFilename,yyyy,netwk,stnm,cmpnm,newFileName);
            unix(cmd);
            disp(cmd);
        end
    end
    cd ..
end