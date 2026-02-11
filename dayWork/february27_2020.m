%february27_2020
clear; close all; clc;
cd ~/products/rsam/
load EC.FER1..BHZ_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat

%%
Ns = [1440,60,10];
ds = 10;
t = getTimeVec(S);
t = datenum(t);
t = downsample(t,ds);
t = datestr(dn2dt(t),31);
TT = string(t);

%%
for j = 1:length(Ns)
    N_ = Ns(j);
    d = int32(pull(medfiltWaveforms(S,N_)));
    d = downsample(d,ds); %TT = [string(t),num2str(d)];
    TT = [TT,num2str(d)];
end
writematrix(["t" "daily" "hourly" "10min"; TT],'~/public_html/ecuadorian-seismicity/js/test2.csv');