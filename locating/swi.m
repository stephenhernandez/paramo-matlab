%function swi
%stein-wysession inversion

clear;
close all;
clc;

singleTrue = true;
load ~/igdata/soam_noec.mat
load ~/igdata/ec_boundaries.mat

if singleTrue
    plotFlag = true;
    iasp91 = load('~/igdata/iasp91_matTauP_2.mat');
    ttp = iasp91.ttp;
    tts = iasp91.tts;
    depths = iasp91.depths;
    dists = iasp91.dists;
    clear iasp91;

    [dTdZ,dTdR] = gradient(ttp,1,0.5);
    [dTsdZ,dTsdR] = gradient(tts,1,0.5);
    [Xorig,Yorig] = meshgrid(depths,dists);
    load ~/igdata/soam.mat
    diasFlag = false;
    %cd ~/regions_old//pedernales/sc3_catalog/
    %eventID = 'igepn2016elzh.txt';
    %eventID = 'igepn2016vkti.txt';
    %eventID = 'igepn2016hpdo.txt';
    %eventID = 'igepn2016hnmu.txt';
    %eventID = 'igepn2016xkjs.txt';
    %eventID = 'igepn2016yyrs.txt';
    %eventID = 'igepn2016yvmv.txt';
    %eventID = 'igepn2016htqj.txt';
    %eventID = 'igepn2017auxv.txt';
    %eventID = 'igepn2016yydj.txt';
    %eventID = 'igepn2016yjjn.txt';
    %eventID = 'igepn2017cdxo.txt';
    %eventID = 'igepn2017howb.txt';
    %eventID = 'igepn2016nedp.txt';
    %eventID = 'igepn2017hztk.txt';
    %eventID = 'igepn2017hsar.txt';
    %eventID = 'igepn2014pbsq.txt';
    %eventID = 'igepn2012ggmt.txt';
    %eventID = 'igepn2016qnts.txt';
    %eventID = 'igepn2017kdtr.txt'; %pichincha event, may 2017
    
    %eventID = 'igepn2022bmbh.txt'; %19h53
    %eventID = "igepn2022blp2.txt"; %13h40
    eventID = "igepn2022fzqp.txt"; %13h40
    
    %eventID = 'igepn2017wovf.txt'; %
    %eventID = 'igepn2018efav.txt'; %

    %eventID = 'igepn2018mkta.txt';
    %eventID = 'igepn2018fmoy.txt';
    %eventID = 'igepn2015bpwq.txt';
    %eventID = 'igepn2018cyni.txt';
    %eventID = 'igepn2015benbernard.txt';
    %eventID = 'igepn2015kdqd.txt';
    %eventID = 'igepn2015qctt.txt';
    %eventID = 'igepn2019clmg.txt'; %event from guayaquil, february 4 2019
    %eventID = 'dias2018mgea.bulletin'; diasFlag = true;
    %eventID = 'igepn2017wnrz.txt';  %
    %eventID = 'igepn2021wyvu.txt';
    [Mbest,minRMS_,minIter,mag,NobsP,origModel]  = eql(eventID,ttp,dTdR,dTdZ,tts,dTsdR,dTsdZ,Xorig,Yorig,plotFlag,diasFlag);
    %save('dumdum')
else
    cd ~/igdata/
    load iasp91_matTauP_2.mat
    [dTdZ,dTdR] = gradient(ttp,1,0.5);
    [dTsdZ,dTsdR] = gradient(tts,1,0.5);
    [Xorig,Yorig] = meshgrid(depths,dists);
    load ~/igdata/soam.mat;
    plotFlag = false;
    magThresh = 3;
    maxMag = 7;
    diasFlag = false;

    [~,~,t,eqlat,eqlon,eqdepth,eqmag,id,stderr,azgap,nPhases,nMLv,...
        timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
        locMethod,earthModel,eqType,creationTime] = readCat1(datetime(2021,10,01),datetime(2030,12,31),magThresh,maxMag);

    eqmag = round(eqmag*100)/100;
    eqmagI = find(eqmag >= magThresh & eqdepth < 200);
    eqmag = eqmag(eqmagI);
    data = [eqmag eqmagI];
    data = sortrows(data,'descend');
    sI = data(:,2);
    cd ~/phaseInformationSC5/2022/;
    %files = dir('igepn2016*.txt');
    %files = importdata('quick_esmeraldas_list.txt');
    %files = importdata('events_gtM5_2016.txt');
    %files = flipud(files);
    %
    files = id(sI);
    Nfiles = length(id);
    Mcurrent = zeros(4,Nfiles);
    Morig = Mcurrent;
    minRMS = zeros(Nfiles,1);
    minIter = minRMS;
    minRMSOrig = minRMS;
    mags = minRMS;
    Nobs = minRMS;
    nnn = 0;
    for i = 1%:Nfiles
        eid = files{i};%(i).name;
        eventID = [eid,'.txt'];
        yyyy = eventID(6:9);
        cd(strcat('~/phaseInformationSC3/',yyyy,'/'));
        %if ~exist(eventID,'file')
        %    disp([eventID,' does not exist, attempting to retrieve...']);
        %    %unix(['~/scripts/get_sc3_phase_information.sh ',eid]);
        %    %disp('done with retrieval attempt');
        %end

        dirdir = dir(eventID);
        if ~exist(eventID,'file')
            fprintf('does not exist, skipping...\n');
            continue;
        end

        if dirdir.bytes == 0
            fprintf('file is empty, skipping...\n');
            continue
        end

        nnn = nnn + 1;
        %eventID = files(i).name; %strcat(ids{i},'.txt');
        %[Mcurrent_,minRMS_,minIter_,mag_,Nobs_,Morig_] = isdd_(eventID);
        [Mcurrent_,minRMS_,minIter_,mag_,Nobs_,Morig_,minRMSOrig_] = eql(eventID,ttp,dTdR,dTdZ,tts,dTsdR,dTsdZ,Xorig,Yorig,plotFlag,diasFlag);
        disp(['FILE NUMBER: ',num2str(nnn)]);
        Mcurrent(:,nnn) = Mcurrent_;
        Morig(:,nnn) = Morig_;
        minRMS(nnn) = minRMS_;
        minIter(nnn) = minIter_;
        mags(nnn) = mag_;
        Nobs(nnn) = Nobs_;
        minRMSOrig(nnn) = minRMSOrig_;
    end

    ecCoast = load('~/igdata/ecuador_coast_islands');
    ecCoastLon = ecCoast.lon;
    ecCoastLat = ecCoast.lat;

    disp('saving data')
    cd ~/igdata
    save('bulk_relocations_2022')

    %%
    close all;
    rmsI = Nobs >= 3 & minRMS > 0 & minRMS < 2 & minRMSOrig > 0 & Mcurrent(4,:)' >= 0 & Mcurrent(3,:)' >= 0 & Mcurrent(2,:)' >= -6 & Mcurrent(1,:)' >= -94;
    %rmsI = Nobs >= 1 & minRMS >= 0 & Mcurrent(4,:)' >= 0 & Mcurrent(3,:)' >= 0 & Mcurrent(2,:)' >= -6 & Mcurrent(1,:)' >= -92;
    rmsI = find(rmsI);

    figure(1);
    ha(1) = subplot(311);
    scatter(Mcurrent(1,rmsI),-Mcurrent(3,rmsI),[],minRMS(rmsI),'filled','marker','s'); colorbar; grid on; zoom on;
    ha(2) = subplot(312);
    scatter(Morig(1,rmsI),-Morig(3,rmsI),[],minRMS(rmsI),'marker','d'); hold on; colorbar; grid on; zoom on;
    figure(1);
    ha(3) = subplot(313); hold on;
    for i = 1:length(rmsI)
        plot([Mcurrent(1,rmsI(i)) Morig(1,rmsI(i))],-[Mcurrent(3,rmsI(i)) Morig(3,rmsI(i))],'k'); hold on;
    end
    scatter(Mcurrent(1,rmsI),-Mcurrent(3,rmsI),[],minRMS(rmsI),'filled','marker','s');
    scatter(Morig(1,rmsI),-Morig(3,rmsI),[],minRMS(rmsI),'marker','d'); hold on;
    colorbar; grid on;
    linkaxes(ha,'xy');

    figure(100);
    plot(Morig(1,rmsI),-Morig(3,rmsI),'ko'); hold on;
    scatter(Mcurrent(1,rmsI),-Mcurrent(3,rmsI),[],minRMSOrig(rmsI)-minRMS(rmsI),'filled','marker','s');
    colorbar;
    grid on; zoom on;

    %color by depth
    figure(2); hold on;
    plot(lonEC,latEC,'k','linewidth',2); hold on;
    plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
    geoshow('~/igdata/fallas/fallas2016.shp','Color','r');
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
    grid on; zoom on;
    scatter(Mcurrent(1,rmsI),Mcurrent(2,rmsI),100,Mcurrent(3,rmsI),'filled','marker','s');
    colorbar;
    axis equal;
    axis([-83 -76 -4 4]);
    title('velocity')
    caxis([-1 25])
    grid on; zoom on;

    %color by depth
    figure(3);
    plot(lonEC,latEC,'k','linewidth',2); hold on;
    plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
    geoshow('~/igdata/fallas/fallas2016.shp','Color','r');
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
    grid on; zoom on;
    for i = 1:length(rmsI)
        plot([Mcurrent(1,rmsI(i)) Morig(1,rmsI(i))],[Mcurrent(2,rmsI(i)) Morig(2,rmsI(i))]); hold on;
    end
    scatter(Mcurrent(1,rmsI),Mcurrent(2,rmsI),[],Mcurrent(3,rmsI),'filled','marker','s');
    scatter(Morig(1,rmsI),Morig(2,rmsI),100,Morig(3,rmsI),'marker','d');
    colorbar;
    axis equal;
    axis([-83 -76 -4 4]);
    title('depth')

    %color by rms
    figure(4);
    plot(lonEC,latEC,'k','linewidth',2); hold on;
    plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
    geoshow('~/igdata/fallas/fallas2016.shp','Color','r');
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
    grid on; zoom on;
    for i = 1:length(rmsI)
        plot([Mcurrent(1,rmsI(i)) Morig(1,rmsI(i))],[Mcurrent(2,rmsI(i)) Morig(2,rmsI(i))]); hold on;
    end
    scatter(Mcurrent(1,rmsI),Mcurrent(2,rmsI),[],minRMS(rmsI),'filled','marker','s');
    scatter(Morig(1,rmsI),Morig(2,rmsI),100,minRMS(rmsI),'marker','d');
    colorbar;
    axis equal;
    axis([-83 -76 -4 4]);
    title('rms')

    %color by rms
    figure(5);
    plot(lonEC,latEC,'k','linewidth',2); hold on;
    plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
    geoshow('~/igdata/fallas/fallas2016.shp','Color','r');
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
    grid on; zoom on;
    for i = 1:length(rmsI)
        plot([Mcurrent(1,rmsI(i)) Morig(1,rmsI(i))],[Mcurrent(2,rmsI(i)) Morig(2,rmsI(i))]); hold on;
    end
    scatter(Mcurrent(1,rmsI),Mcurrent(2,rmsI),[],mags(rmsI),'filled','marker','s');
    scatter(Morig(1,rmsI),Morig(2,rmsI),100,mags(rmsI),'marker','d');
    colorbar;
    axis equal;
    axis([-83 -76 -4 4]);
    title('magnitude');

    figure;
    scatter(Mcurrent(4,rmsI),-Mcurrent(3,rmsI),[],mags(rmsI),'filled');
    colorbar; caxis([4 8]);
    grid on; zoom on;

    figure;
    scatter(Mcurrent(4,rmsI),-Mcurrent(3,rmsI),[],minRMS(rmsI),'filled');
    colorbar;
    grid on; zoom on;

    figure; plot(Nobs(rmsI),'o')
    grid on; zoom on;

    %
    figure;
    ampFact = -1;
    plot(lonEC,latEC,'k','linewidth',2); hold on;
    plot(lon_noec,lat_noec,'linewidth',2,'Color',[0.5 0.5 0.5]);
    geoshow('~/igdata/fallas/fallas2016.shp','Color','r');
    geoshow('~/igdata/ZonaUrbana/ZonaUrbana.shp');
    grid on; zoom on;
    quiver(Mcurrent(1,rmsI),Mcurrent(2,rmsI),ampFact*(Mcurrent(1,rmsI)-Morig(1,rmsI)),ampFact*(Mcurrent(2,rmsI)-Morig(2,rmsI)),'linewidth',1);
    scatter(Mcurrent(1,rmsI),Mcurrent(2,rmsI),[],mags(rmsI),'filled','marker','o');
    colorbar;
    caxis([3.5 5.5])
    axis equal;
    axis([-83 -76 -4 4]);
end
