%%
cd ~/products/rsam/

AGETHRESHOLD = 180; % oldest age in days

%now1 = now;
%now2 = now;
previousUTCDay = dateshift(dn2dt(now)+hours(5),'start','day');


while true
    currentUTCDay = dateshift(dn2dt(now)+hours(5),'start','day');
    oldestThreshold = currentUTCDay - AGETHRESHOLD;
    [yyyy,~,~] = datevec(currentUTCDay);
    doy = day(currentUTCDay,'doy');
    
    if doy < 10
        doyStr = ['00',num2str(doy)];
    elseif doy < 100
        doyStr = ['0',num2str(doy)];
    else
        doyStr = num2str(doy);
    end
    
    %%
    files = dir('*.mat');
    lFiles = length(files);
    Sold = populateWaveforms(lFiles);
    for i = 1:lFiles
        fname = files(i).name;
        
        %%
        disp(fname)
        try
            Sold(i) = load(fname,'S');
        catch
            warning('Couldnt load %s',fname);
        end
        
        %% at this point, the file has been read, and the data are in memory
        newStart = dateshift(S.ref + S.e,'start','day'); % assume UTC, what day do i want?
        
        %% 
        if newStart >= oldestThreshold
            tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
            S_ = rmsGather(newStart,tEnd,60,1,8,S(1).kstnm,S(1).kcmpnm,0,1,S(1).knetwk,S(1).khole,0,1);
            [M,status] = mergeWaveforms([Sold(i); S_]);
            if status
                disp('good');
                S = M(1);
                save(fname,'S');
            end
        end
    end
    
    
    cd ~/products/events/html/
    fNames = importdata('~/importFileNames.txt');
    fNames = string(fNames);
    for i = 1:length(fNames)
        disp(i);
        load(fNames(i));
        id_ = regexp(fNames(i),'\/','split');
        id(i) = id_(1);
        W = table2struct(W);
        W = differentiateWaveforms(W);
        W = filterWaveforms(W,0.75,12,4,0.02);
        plotRecordSection(W,-inf,-inf,false,true,id(i),'Z',false);
    end
    previousUTCDay = currentUTCDay;
end

%%
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
% % % clear; close all; clc; 
% % % lfc = 1;
% % % hfc = 4;
% % % cwiStart = 10;
% % % cwiDur = 80;
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
