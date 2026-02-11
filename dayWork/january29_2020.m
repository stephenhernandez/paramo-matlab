% january29_2020
clear; close all; clc;
nets = "EC";
lNets = length(nets);

chans = ["HHZ"; "SHZ"; "BHZ"; "HNZ"; "BLZ"; "ENZ"; "BDF"];
lChans = length(chans);

%%
lfc = [1/8; 1/4; 1/2; 1; 2; 4; 8];
bwOctaves = 2.^[1; 2; 3; 4; 5; 6; 7; 8];
newFs = 100;

%%
currentDay = dateshift(dn2dt(now)+hours(5),'start','day'); %add 5 hours for UTC
disp(currentDay);
%dn2dt(now); %dn2dt(floor(now));
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
                    miniseedfile = dir(strcat(f,'/',net_,'.*.',doyStr));
                    if ~isempty(miniseedfile)
                        lminiseedfile = length(miniseedfile);
                        for l = 1:lminiseedfile
                            mseedname = miniseedfile(l).name;
                            
                            splits = split(mseedname,".");
                            knetwk = splits(1);
                            kstnm = splits(2);
                            khole = splits(3);
                            kcmpnm = splits(4);
                            
                            if strcmp(kcmpnm,"BDF") && (strcmp(khole,"02") || strcmp(khole,"03"))
                                continue;
                            end
                            
                            if strcmp(kcmpnm,"ENZ") || strcmp(kcmpnm,"HNZ")
                                fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_60DUR_ACC_RMS_SLIDEMED_MedPreFilt3Points.mat"));
                            elseif strcmp(kcmpnm,"BDF")
                                fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_4sec1Hz_60DUR_BDF_RMS_SLIDEMED_MedPreFilt3Points.mat"));
                            else
                                fname = fullfile('~/products/rsam/',strcat(knetwk,".",kstnm,".",khole,".",kcmpnm,"_1Hz8Hz_60DUR_VEL_RMS_SLIDEMED_MedPreFilt3Points.mat"));
                            end
                            
                            if ~exist(fname,'file')
                                %% run rmsGather
                                if strcmp(kcmpnm,"BDF")
                                    S = rmsGather(datetime(2006,01,01),dateshift(dn2dt(now)+hours(5),'end','day'),60,1/4,1,kstnm,kcmpnm,0,1,knetwk,khole,0,1);
                                else
                                    S = rmsGather(datetime(2006,01,01),dateshift(dn2dt(now)+hours(5),'end','day'),60,1,8,kstnm,kcmpnm,0,1,knetwk,khole,0,1);
                                end
                                save(fname,'S');
                            else
                                %% update already existing
                                try
                                    load(fname,'S');
                                catch
                                    warning('Couldnt load %s',fname);
                                end
                                
                                %%
                                newStart = dateshift(S.ref + S.e,'start','day');
                                if newStart >= datetime(2020,01,01)
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
