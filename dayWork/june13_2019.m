%june13_2019
clear; close all; clc;
cd ~/research/now/sn_sws/

data = load('orig_cat.txt'); %cut1.txt');
t = datetime(data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,7));
eqmag = data(:,11);
mI = eqmag >= 3;
t = t(mI);

%%
kstnms =    ["SN02","SN06"]; %,"SN05","SN04","SN11","SN12","SN13","SN14","SN07","FER2","PVIL","ALCE","VCH1"];
knetwks =   ["9D",  "9D"]; %,   "9D", "9D",  "9D",  "9D",  "9D",  "9D",  "EC",  "EC",  "EC",  "EC",  "EC"];
%kstnms = ["CEAZ","FER1"];

for k = 1:length(kstnms)
    kstnm = kstnms(k); %"VCH1";
    knetwk = knetwks(k); %"EC";
    khole = "";
    if strcmp (kstnm,"SN02") || strcmp(kstnm,"SN06") || strcmp(kstnm,"SN09") || strcmp(kstnm,"SN10")
        E = extractSacFromList(t-seconds(15),45,kstnm,"HH2",knetwk);
        lE = length(E);
        for j = 1:lE
            E(j).kcmpnm = "HHE";
        end
        N = extractSacFromList(t-seconds(15),45,kstnm,"HH1",knetwk);
        lN = length(N);
        for j = 1:lN
            N(j).kcmpnm = "HHN";
        end
    else
        E = extractSacFromList(t-seconds(15),45,kstnm,"HHE",knetwk);
        N = extractSacFromList(t-seconds(15),45,kstnm,"HHN",knetwk);
    end
    Z = extractSacFromList(t-seconds(15),45,kstnm,"HHZ",knetwk);
    
    %%
    cd ~/research/now/sn_sws/sac/
    lT = length(t);
    A = [E,N,Z];
    for i = 1:lT
        ref = E(i).ref;
        [yyyy_,month_,day_,hr_,min_,sec_] = datevec(ref);
        msec_ = round((sec_ - floor(sec_))*1000);
        sec_ = floor(sec_);
        
        monthStr = num2str(month_);
        dayStr = num2str(day_);
        hrStr = num2str(hr_);
        minStr = num2str(min_);
        secStr = num2str(sec_);
        msecStr = num2str(msec_);
        
        if month_ < 10; monthStr = ['0',monthStr]; end
        if day_ < 1; dayStr = ['0',dayStr]; end
        if hr_ < 10; hrStr = ['0',hrStr]; end
        if min_ < 10; minStr = ['0',minStr]; end
        if sec_ < 10; secStr = ['0',secStr]; end
        if msec_ < 10; msecStr = ['00',msecStr]; end
        if msec_ < 100; msecStr = ['0',msecStr]; end
        
        [tyyyy_,tmonth_,tday_,thr_,tmin_,tsec_] = datevec(t(i));
        tmsec_ = round((tsec_ - floor(tsec_))*1000);
        tsec_ = floor(tsec_);
        
        tmonthStr = num2str(tmonth_);
        tdayStr = num2str(tday_);
        thrStr = num2str(thr_);
        tminStr = num2str(tmin_);
        tsecStr = num2str(tsec_);
        tmsecStr = num2str(tmsec_);
        
        if tmonth_ < 10; tmonthStr = ['0',tmonthStr]; end
        if tday_ < 1; tdayStr = ['0',tdayStr]; end
        if thr_ < 10; thrStr = ['0',thrStr]; end
        if tmin_ < 10; tminStr = ['0',tminStr]; end
        if tsec_ < 10; tsecStr = ['0',tsecStr]; end
        if tmsec_ < 10; tmsecStr = ['00',tmsecStr]; end
        if tmsec_ < 100; tmsecStr = ['0',tmsecStr]; end
        
        A_ = A(i,:);
        A_ = synchSacData(A_);
        newDir = ['Event_',num2str(tyyyy_),'.',tmonthStr,'.',tdayStr,'.',thrStr,...
            '.',tminStr,'.',tsecStr,'.',tmsecStr];
        
        disp(newDir)
        if ~exist(newDir,'dir')
            mkdir(newDir)
        end
        cd(newDir)
        
        for j = 1:3
            chanStr = A_(1,j).kcmpnm;
            fileName = [char(kstnm),'.',char(knetwk),'.',char(khole),'.',...
                char(chanStr),'.',num2str(yyyy_),'.',monthStr,'.',dayStr,'.',...
                hrStr,'.',minStr,'.',secStr,'.',msecStr,'.SAC'];
            mksac(fileName,A_(1,j).d,datenum(A_(1,j).ref),'DELTA',A_(1,j).delta,...
                'KSTNM',A_(1,j).kstnm,'KCMPNM',A_(1,j).kcmpnm,...
                'KNETWK',A_(1,j).knetwk);
        end
        cd ..
    end
end
