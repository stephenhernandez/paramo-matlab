clear; close all; clc;

%%
cd ~/products/rsam/
files = dir('EC.*');
lFiles = length(files);
for i = 1:lFiles
    str = split(files(i).name,'.');
    sensorName(i,1) = string(str{2});
end

%%
alphaLevel = 0.25;

%%
uniqSensors = unique(sensorName);
lSensors = length(uniqSensors);

%%
for i = 152%:lSensors
    uniqSensor_ = uniqSensors(i);
    disp(uniqSensor_);
    
    %%
    pzStart = dn2dt(now) + hours(5);
    
    %%
    cd ~/response;
    unixCommand = ['bash getResponseInfo.sh ',char(uniqSensor_),' EC'];
    status = unix(unixCommand);
    if exist(strcat(char(uniqSensor_),'_EC'),'dir')
        cd(strcat(char(uniqSensor_),'_EC'));
        pzFiles = dir('SAC_PZs_EC*');
        lPZs = length(pzFiles);
        for k = 1:lPZs
            pz_ = pzFiles(k).name;
            pzSplit = split(pz_,'_');
            kcmpnm = pzSplit{5};
            if strcmp(kcmpnm(3),'Z')
                [~,~,~,~,~,~,~,pzStart_] = read_sac_pole_zero(pz_);
                pzStart = min([pzStart pzStart_]);
            end
        end
    end
    cd ~/products/rsam/
    sensorFiles = dir(strcat('EC.',char(uniqSensor_),'.*'));
    nFiles = length(sensorFiles);
    
    %%
    wStart = dn2dt(now) + hours(5);
    wEnd = wStart;
    
    %%
    W = populateWaveforms(nFiles);
    for j = 1:nFiles
        load(sensorFiles(j).name,'S');
        if ~isnat(S.ref)
            W(j) = S;
            wStart = min([wStart S.ref]);
        end
    end
    
    %%
    if ~isnat(W(1).ref)
        cutStart = min([pzStart wStart]);
        refs = pull(W,'ref');
        badRefs = isnat(refs);
        W(badRefs) = [];
        
        %%
        close all;
        [~,ax] = plotWaveforms(medfiltWaveforms(W,101),...
            [],[],[],[],[],[],[],...
            'on');
        xlim([cutStart wEnd]);
        
        %%
        nFiles = length(W);
        for j = 1:nFiles
            thisStart = W(j).ref;
            if ~isnat(thisStart)
                ax(j).YScale = 'log';
                ylim_ = ax(j).YLim;
                ylim_(1) = max([ylim_(1) 1e-1]);
                ax(j).YLim = ylim_;
                if cutStart < thisStart
                    hold(ax(j),'on');
                    aPlot = area(ax(j),[cutStart thisStart],[ylim_(2) ylim_(2)]);
                    aPlot.LineStyle = 'none';
                    aPlot.FaceAlpha = alphaLevel;
                end
            end
        end
        
        %%
        %suptitle();
        %%
%         fname1 = fullfile('~','products','station_health',[char(uniqSensor_),'_EC.jpg']);
%         print('-djpeg',fname1);
%         pause(5);
%         close all;
    end
end
