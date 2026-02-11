clear;
kstnm = "BREF";
kcmpnm = "BHZ";
knetwk = "EC";
khole = "";
fname1 = fullfile("~/products/rsam/",strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4sec2Hz_60DUR_RootMeanSquare_preFiltFalse.mat"));
fname2 = fullfile("~/products/rsam/",strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_2Hz8Hz_60DUR_RootMeanSquare_preFiltFalse.mat"));

secDur = 60;
tStart = datetime(2006,05,01);
tEnd = datetime(2022,10,22);

try
    RMS_4sec2Hz = rmsGather(tStart,tEnd,secDur,0.25,2,kstnm,kcmpnm,true,true,knetwk);
    save(fname1,"RMS_4sec2Hz");
    save("tmp1.mat","fname1");
catch
    fprintf("4Sec2Hz failed for some reason\n");
end

try
    RMS_2Hz8Hz = rmsGather(tStart,tEnd,secDur,2,8,kstnm,kcmpnm,true,true,knetwk);
    save(fname2,"RMS_2Hz8Hz");
    save('tmp2.mat',"fname2");
catch
    fprintf("2Hz8Hz failed for some reason\n");
end