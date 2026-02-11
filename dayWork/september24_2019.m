clear; clear all; close all; clc;
lfc = 1;
hfc = 4;
cwiStart = 10;
cwiDur = 70;
allFlag = 0;
newFs = 32;
secDur = 2^12;
maxLag = 128;
dailyFlag = true;

cd ~/research/now/dV/tungurahua/brtu
[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHZ","EC",""],["BRTU","HHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHZHHZ_32sps_OneBit_4096SecDur','-v7.3');

[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHN","EC",""],["BRTU","HHN","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHNHHN_32sps_OneBit_4096SecDur','-v7.3');

[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHE","EC",""],["BRTU","HHE","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHEHHE_32sps_OneBit_4096SecDur','-v7.3');

[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHZ","EC",""],["BRTU","HHE","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHZHHE_32sps_OneBit_4096SecDur','-v7.3');

[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHZ","EC",""],["BRTU","HHN","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHZHHN_32sps_OneBit_4096SecDur','-v7.3');

[ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag,dayStack] = dV(datetime(2019,05,01),datetime(2019,07,01),datetime(2019,05,26),["BRTU","HHE","EC",""],["BRTU","HHN","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
save('BRTU_May2019June2019_1Hz4Hz_10s80s_HHEHHN_32sps_OneBit_4096SecDur','-v7.3');


% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,10,01),datetime(2020,01,01),["BMAS","BHZ","EC",""],["BMAS","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,true);
% cd ~/research/now/dV/tungurahua/bmas
% save('BMAS_2006_2019_1Hz4Hz_10s80s','-v7.3');

% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["BPAT","BHZ","EC",""],["BPAT","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('BPAT_2006_2019_1Hz4Hz_10s80s','-v7.3');
% 
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["BRUN","BHZ","EC",""],["BRUN","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('BRUN_2006_2019_1Hz4Hz_10s80s','-v7.3');
% 
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["BULB","BHZ","EC",""],["BULB","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('BULB_2006_2019_1Hz4Hz_10s80s','-v7.3');
% 
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["BBIL","BHZ","EC",""],["BBIL","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('BBIL_2006_2019_1Hz4Hz_10s80s','-v7.3');
% 
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["POND","HHZ","EC",""],["POND","HHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('POND_2006_2019_1Hz4Hz_10s80s','-v7.3');
% 
% [ccMaster,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags,dVsymStack,ccCaus,ccAcaus,ccSymmetric,refTrace,allFlag] = dV(datetime(2006,01,01),datetime(2019,09,01),datetime(2020,01,01),["BMAS","BHZ","EC",""],["BMAS","BHZ","EC",""],lfc,hfc,cwiStart,cwiDur,allFlag,newFs,secDur,maxLag,dailyFlag);
% cd ~/research/now/dV/tungurahua/bmas
% save('BMAS_2006_2019_1Hz4Hz_10s80s','-v7.3');