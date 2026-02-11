S = rmsGather(datetime(2009,08,01),datetime(2022,08,20),60,10,-inf,...
    "SAGA","HHZ",false,true,"EC","",false,false,[],true);

cd ~/products/rsam/temporary/
clearvars -except S

save('EC.SAGA..HHZ_10Hz_60DUR_MedRMSAmp_preFiltFalse.mat','S');