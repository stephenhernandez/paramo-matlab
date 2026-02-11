clear; close all; clc; 
cd ~/products/rsam;

S = rmsGather(datetime(2021,01,01),dn2dt(ceil(now)),60,10,-inf,"SAGA","HHZ",true,false,"EC",""); % true-false is mean of abs value
save('EC.SAGA..HHZ_10Hz_60DUR_MeanAbsAmp_preFiltFalse_2019.mat','S');