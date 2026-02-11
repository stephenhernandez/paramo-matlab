clear; close all; clc;
cd ~/products/rsam/

tStart = datetime(2018,11,01);
tEnd = datetime(2023,04,28);
secDur = 60;
meanFlag = true;
rmsFlag = false;

kstnms = ["PUYO";"TAIS";"PORT";"TAMH";"PKYU";"SAGA";"SAGA"];
for i = 1:length(kstnms)
    kstnm_ = kstnms(i);
    if ~strcmp(kstnm_,"SAGA")
        S = rmsGather(tStart,tEnd,secDur,0.6,1.2,kstnm_,"HHZ",meanFlag,rmsFlag,"EC");
        save(lower(kstnm_),'S');
    else
        if i == 6
            S = rmsGather(tStart,tEnd,secDur,0.25,2,kstnm_,"HHZ",meanFlag,rmsFlag,"EC");
            save('saga1','S');
        else
            S = rmsGather(tStart,tEnd,secDur,2,8,kstnm_,"HHZ",meanFlag,rmsFlag,"EC");
            save('saga2','S');
        end
    end
end