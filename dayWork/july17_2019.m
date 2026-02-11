%july17_2019
clear; close all; clc;
cd ~/research/favors/gponce/
e = importdata('events.txt');
e = string(e);
startTime = -10;
endTime = 5*60 - startTime;
khole = '';
le = length(e);
E = populateSeisCompStructure(le);
for i = 1:le
    eventName = [char(e(i)),'.txt'];
    yyyy = str2double(eventName(6:9));
    fName = ['~/phaseInformationSC3/',num2str(yyyy),'/',eventName];
    E(i) = readSCBulletin(fName);
    originTime = E(i).t;
    requestedDay = dateshift(originTime,'start','day');
    disp(fName);
    id = E(i).id;
    if ~exist(id,'dir')
        mkdir(id)
    end
    cd(id)
    Pphases = E(i).Pphases;
    for j = 1:length(Pphases)
        Pphases_ = Pphases(j);
        stnm = Pphases_.stnm;
        ntwk = Pphases_.ntwk;
        chan = Pphases_.chan;
        chan = chan(1:2);
        for k = ['E','N','Z']
            chanStr = string([chan,k]);
            S = loadWaveforms(originTime,1,string(stnm),chanStr,string(ntwk));
            if ~isnat(S.ref)
                Scut = cutWaveforms(S,originTime,startTime,endTime);
                [yyyy_,month_,day_,hr_,min_,sec_] = datevec(Scut.ref);
                msec_ = round((sec_ - floor(sec_))*1000);
                sec_ = floor(sec_);
                
                monthStr = num2str(month_);
                dayStr = num2str(day_);
                hrStr = num2str(hr_);
                minStr = num2str(min_);
                secStr = num2str(sec_);
                msecStr = num2str(msec_);
                
                if month_ < 10; monthStr = ['0',monthStr]; end
                if day_ < 10; dayStr = ['0',dayStr]; end
                if hr_ < 10; hrStr = ['0',hrStr]; end
                if min_ < 10; minStr = ['0',minStr]; end
                if sec_ < 10; secStr = ['0',secStr]; end
                if msec_ < 10; msecStr = ['00',msecStr]; end
                if msec_ < 100; msecStr = ['0',msecStr]; end
                
                fileName = [stnm,'.',ntwk,'.',char(khole),'.',...
                    char(chanStr),'.',num2str(yyyy_),'.',monthStr,'.',dayStr,'.',...
                    hrStr,'.',minStr,'.',secStr,'.',msecStr,'.SAC'];
                mksac(fileName,Scut.d,datenum(Scut.ref),'DELTA',Scut.delta,...
                    'KSTNM',Scut.kstnm,'KCMPNM',Scut.kcmpnm,...
                    'KNETWK',Scut.knetwk);
            end
        end
    end
    cd ..
end
