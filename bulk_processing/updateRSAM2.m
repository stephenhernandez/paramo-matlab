%function updateRSAM()
clear; close all; clc;
nets = "EC";
lNets = length(nets);

chans = ["HHZ"; "SHZ"; "BHZ"; "HNZ"; "BLZ"; "ENZ"; "BDF"];
lChans = length(chans);
meanFlag = false; %when false, we take the median amplitude in the window
rmsFlag = true;

%%
currentDay = dateshift(dn2dt(now)+hours(5),'start','day'); %add 5 hours for UTC
disp(currentDay);

[yyyy,~,~] = datevec(currentDay);
doy = day(currentDay,'doy');

if doy < 10
    doyStr = ['00',num2str(doy)];
elseif doy < 100
    doyStr = ['0',num2str(doy)];
else
    doyStr = num2str(doy);
end

%%
invN = 0;
secDur = 30;
if exist('~/rawdata','dir') == 7
    cd('~/rawdata/')
    cd(num2str(yyyy));
    for i = 1:lNets
        dirs = dir(nets(i));
        dirNames = string({dirs.name})';
        goodI = ([dirs.isdir])' & ~strcmp(dirNames,".") & ~strcmp(dirNames,"..");
        dirNames = dirNames(goodI);
        lNames = sum(goodI);
        if lNames
            net_ = nets(i);
            cd(net_);
            for j = 1:lNames
                for k = 1:lChans
                    chanStr = strcat(chans(k),'.D');
                    f = fullfile(pwd,dirNames(j),chanStr);
                    miniseedfile = dir(strcat(f,'/',net_,'.*.',doyStr)); %this line is potentially slow
                    if ~isempty(miniseedfile)
                        lminiseedfile = length(miniseedfile);
                        for l = 1:lminiseedfile
                            mseedname = miniseedfile(l).name;
                            mseedFullfile = fullfile(miniseedfile(l).folder,mseedname);
                            %S_ = readMiniSeed(mseedFullfile);
                            splits = split(mseedname,".");
                            knetwk = splits(1);
                            kstnm = splits(2);
                            khole = splits(3);
                            kcmpnm = splits(4);
                            
                            %% select only "01" BDF channel
                            if strcmp(kcmpnm,"BDF") && (strcmp(khole,"02") || strcmp(khole,"03"))
                                continue;
                            end
                            invN = invN + 1;
                            disp(invN);
                            
                            %%
                            if strcmp(kcmpnm,'BDF')
                                lfc = 1/4;
                                hfc = 1;
                            else
                                lfc = 1;
                                hfc = 8;
                            end
                            
                            %%
                            %                             if ~strcmp(kcmpnm,"BDF")
                            %                                 fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                            %                             else
                            %                                 fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4sec1Hz_30DUR_MedAmpRMS_preFiltFalse.mat"));
                            %                             end
                            
                            searchQuery = strcat('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_*"));
                            files = dir(searchQuery);
                            lfiles = length(files);
                            
                            if lfiles
                                for m = 1:lfiles
                                    fname_ = fullfile(files(m).folder,files(m).name);
                                    %disp(fname_);
                                    load(fname_,'lastPoint');
                                    if exist('lastPoint','var')
                                        %disp(lastPoint);
                                        clearvars lastPoint
                                    else
                                        load(fname_,'S');
                                        lastPoint = S.ref + S.e;
                                        %disp(lastPoint);
                                        save(fname_,'lastPoint','-append');
                                        clearvars lastPoint
                                    end
                                end
                            end
                            
                            %%
                            %                             nowTime = dn2dt(now) + hours(5);
                            %                             tStart = datetime(1997,01,01);
                            %                             tEnd = dateshift(nowTime,'start','day');
                            
                            
                            %S(n) = rsamWaveforms(S_,dur,lfc,hfc,npoles,meanFlag,rmsFlag,diffFlag,0,0,medPreFiltFlag);
                            %S =  rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %                             if ~exist(fname,'file')
                            %                                 S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %                                 save(fname,'S');
                            %                                 clear S;
                            %                             else
                            %                                 %% update already existing
                            %                                 try
                            %                                     load(fname,'S');
                            %                                 catch
                            %                                     warning('Couldnt load %s',fname);
                            %                                     continue;
                            %                                 end
                            %
                            %                                 %%
                            %                                 lS = length(S);
                            %                                 d = S(1).d;
                            %                                 if isnat(S(1).ref) || isempty(d) || lS > 1
                            %                                     S = rmsGather(tStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %                                     save(fname,'S');
                            %                                     clear S;
                            %                                 else
                            %                                     newStart = dateshift(S.ref + S.e,'start','day');
                            %                                     if newStart >= datetime(2019,01,01)
                            %                                         tEnd = dateshift(dn2dt(now)+hours(5),'start','day');
                            %                                         S_ = rmsGather(newStart,tEnd,secDur,lfc,hfc,kstnm,kcmpnm,meanFlag,rmsFlag,knetwk,khole);
                            %
                            %                                         [M,status] = mergeWaveforms([S; S_]);
                            %                                         if status
                            %                                             S = M(1);
                            %                                             save(fname,'S');
                            %                                         end
                            %                                         clear S_ M S;
                            %                                     else
                            %                                         fprintf('end time is %s, will not attempt to update %s\n',datestr(S.ref+S.e),fname);
                            %                                     end
                            %                                 end
                            %                             end
                            
                        end
                    end
                end
            end
        else
            disp('there are no valid directories for this network');
        end
    end
else
    disp('rawdata data directory is not mounted.');
end
