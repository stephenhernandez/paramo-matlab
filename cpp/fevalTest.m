% cd ~/matlab/
% PWD = pwd;
% PWD = strcat(PWD,'/working/');
% 
% setenv('TZ','America/Guayaquil');
% datetime.setDefaultFormats('default','dd-MMM-uuuu HH:mm:ss.SSS');
% javaaddpath([PWD,'/IRIS-WS-2.0.15.jar']);
% javaaddpath([PWD,'/matTaup/lib/matTaup.jar']);
% format long g;
% addpath(genpath(PWD),'-end');
% rng('shuffle');
% 
% set(groot,'defaultAxesFontSize',18);
% set(groot,'defaultColorbarFontSize',12);
% set(groot,'defaultLineLineWidth',0.1);
% set(groot,'defaultScatterLineWidth',0.1);
% set(groot,'defaultStemLineWidth',0.1);
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaultTextInterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% set(groot,'defaultColorbarTickLabelInterpreter','latex');
% set(groot,'defaultTextarrowshapeInterpreter','latex');
% set(groot,'defaultTextboxshapeInterpreter','latex');
% set(groot,'defaultLineMarkerSize',8);
% set(groot,'defaultFigurePaperPositionMode','auto');
% set(groot,'defaultFigureRenderer','Painters');
% set(groot,'defaultFigureColor',[1 1 1]);

cd ~/;

%%
%january29_2020();

%%
clear; close all; clc;
cd ~/products/rsam/

% n = 0;
% sensors = string(importdata('potentialSensors.txt'));
% for i = 1:length(sensors)
%     sensor_ = sensors(i);
%     splits = split(sensor_,".");
%     knetwk = splits(1);
%     kstnm = splits(2);
%     khole = splits(3);
%     kcmpnm = splits(4);
%     
%     if strcmp(kcmpnm,"ENZ") || strcmp(kcmpnm,"HNZ")
%         fname2 = strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat");
%     else
%         fname2 = strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat");
%     end
%     
%     if ~exist(fname2,'file')
%         disp('here i would run rmsGather');
%         disp(fname2);
%         S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,1,8,kstnm,kcmpnm,0,1,knetwk,khole,0,1);
%         save(fname2);
%         
%         n = n+1;
%         disp(n);
%     else
%         disp(' '); %oops, file already exists');
%     end
% end

%% format
% 'EC.FENY..HNZ_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat'
% S = rmsGather(datetime(2006,01,01),dn2dt(floor(now)),60,1,8,"FENY","HNZ",0,1,"EC","",0,1);

%%
files = dir('*.mat');
for i = 1:length(files)
    fname = files(i).name;
    
    %%
    disp(fname)
    try
        load(fname);
    catch
        warning('Couldnt load %s',fname);
    end
    
    %%
    newStart = dateshift(S.ref + S.e,'start','day');
    if newStart >= datetime(2018,01,01)
        tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
        S_ = rmsGather(newStart,tEnd,60,1,8,S(1).kstnm,S(1).kcmpnm,0,1,S(1).knetwk,S(1).khole,0,1);
        [M,status] = mergeWaveforms([S; S_]);
        if status
            S = M(1);
            save(fname,'S');
        end
        clear S_;
    else
        disp('too many gaps')
    end
end

%%
% % % %%
% % % %cd('~/products/rsam');
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"LAV4","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"CONE","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"REVN","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"REVS","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"GGPC","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"GGPT","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"JUI6","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"RUN5","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"CHL1","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"CHL2","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"LNGL","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"TERV","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12)),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"HPAL","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"GGP","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"PULU","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"CAYR","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"ANGU","SHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %S = rmsGather(datetime(2006,01,01),datetime(2019,12,31),60,0.75,12,"CAYA","HHZ",false,false,"EC","",false,true);
% % % %clearvars -except S
% % % %S.d = single(S.d);
% % % %save(strcat(S.knetwk,".",S.kstnm,".",S.khole,".",S.kcmpnm,"_",num2str(0.75),"_",num2str(12),".mat"),'S','-v7.3');
% % % %
% % % %%
% % % %clear
% % % %cd ~/products/events/html/
% % % %fNames = importdata('~/importFileNames.txt');
% % % %fNames = string(fNames);
% % % %for i = 1:length(fNames)
% % % %    disp(i);
% % % %    load(fNames(i));
% % % %    id_ = regexp(fNames(i),'\/','split');
% % % %    id(i) = id_(1);
% % % %    W = table2struct(W);
% % % %    W = differentiateWaveforms(W);
% % % %    W = filterWaveforms(W,0.75,12,4,0.02);
% % % %    plotRecordSection(W,-inf,-inf,false,true,id(i),'Z',false);
% % % %end
% % % %
% % % %%%
% % % %clear
% % % %cd ~/products/events/html/
% % % %fNames = importdata('~/ispt_events.txt');
% % % %fNames = string(fNames);
% % % %for i = 1:length(fNames)
% % % %    disp(i);
% % % %    [Sall,E,goodI] = loadEventWaveforms(fNames(i),-5,45,true,true,true,true);
% % % %end
% % % %
% % % %%
% % % % clear; close all; clc; 
% % % % lfc = 1;
% % % % hfc = 4;
% % % % cwiStart = 10;
% % % % cwiDur = 80;
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2009,05,31),datetime(2010,05,31),datetime(2010,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_1','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2010,05,31),datetime(2011,05,31),datetime(2011,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_2','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2011,05,31),datetime(2012,05,31),datetime(2012,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_3','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2012,05,31),datetime(2013,05,31),datetime(2013,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_4','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2013,05,31),datetime(2014,05,31),datetime(2014,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_5','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2014,05,31),datetime(2015,05,31),datetime(2015,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_6','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2015,05,31),datetime(2016,05,31),datetime(2016,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_7','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2016,05,31),datetime(2017,05,31),datetime(2017,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_8','-v7.3');
% % % 
% % % [dayStack,dVstack,dVcaus,dVacaus,dVsymmetric,t,lags] = dV(datetime(2017,05,31),datetime(2018,05,31),datetime(2018,05,31),["RETU","SHZ","EC",""],["RETU","SHZ","EC",""],lfc,hfc,cwiStart,cwiDur,0);
% % % cd ~/research/now/tungurahua/dv/
% % % save('ictpRun_9','-v7.3');
% % % 
% % % %%
% % % clear; close all; clc;
% % % [t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
% % %     timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
% % %     locMethod,earthModel,eqType,creationTime] = readCat1(datetime(2018,01,01),datetime(2020,01,01),2);
% % % roundDaysStart = unique(dateshift(t,'start','day'));
% % % for i = 1:length(roundDaysStart)
% % %     n(i) = sum(t >= roundDaysStart(i) & t < roundDaysStart(i)+1);
% % % end
% % % weightedN = datenum(roundDaysStart-min(roundDaysStart));
% % % 
% % % [~,bI] = sort(weightedN,'descend');
% % % bI(1:10)
% % % roundDaysStart(bI(1:10))
% % % 
% % % for i = 1:10000
% % %     disp(i);
% % %     tI = t >= roundDaysStart(bI(i)) & t < roundDaysStart(bI(i))+1;
% % %     ids_ = ids(tI);
% % %     disp(ids_)
% % %     [Sall,E,goodI] = loadEventWaveforms(ids_,-10,50,true,true,true,true);
% % %     plotEvents(E,false);
% % % end
% % % 
% % % %%
% % % % clear; close all; clc; [t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,nPhases,nMLv,...
% % % % timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
% % % % locMethod,earthModel,eqType,creationTime] = readCat1();
% % % % [eqmag,magI] = sort(eqmag,'descend');
% % % % ids = ids(magI);
% % % % eI = eqmag >= 3;
% % % % ids = ids(eI);
% % % % E = populateSeisCompStructure(sum(eI));
% % % % parfor i = 1:sum(eI)
% % % % disp(i); E(i) = readSCBulletin(ids(i));
% % % % end
% % % % close all; 
% % % % 
% % % % for i = 1:length(E)
% % % % plotEvents(E(i),false);
% % % % end
