function [strikeO,dipO,rakeO,conf0] = SC2HASH(E,plotFlag,printFlag)
if nargin < 2; plotFlag = false; end
if nargin < 3; printFlag = false; end

baseDir = "~/masa/old/research/now/focmechs/";
%clear; close all; clc;
%E = readSCBulletin('~/phaseInformationSC3/2019/igepn2019kjzn.txt');
%E = readSCBulletin("igepn2019kjzn");
%E = readSCBulletin("igepn2025dtkg");
%E = readSCBulletin("igepn2025dzol");
%E = readSCBulletin("igepn2025ceos");
%plotFlag = false; printFlag = false;
%SC2HASH(E);

%%
t = pull(E,'t');
lat = pull(E,'lat');
lon = pull(E,'lon');
depth = pull(E,'depth');
depth = depth + 1.1;
laterr = pull(E,'laterr');
lonerr = pull(E,'lonerr');
verErr = pull(E,'deptherr');
mag = pull(E,'mag');
horizErr = 1*sqrt(laterr.^2 + lonerr.^2);

%%                                                                   %lp
formatSpec1 = '%04d%02d%02d%02d%02d%5.2f%02d%s%5.2f%03d%s%5.2f%5.2f%s%02d%s%5.2f%s%5.2f%s%4.2f%s%s';
formatSpec2 = '%4s %2s  %3s %s %s';
formatSpec2a = '%72s';
formatSpec3 = '%4s %3s %s %9.5f %10.5f %5d %s %2s'; %a4,1x,a3,33x,f9.5,1x,f10.5,1x,i5,23x,a2
formatSpec4 = '%4s %d        %d';

%%
N = length(t);
latSign = ' ';
lonSign = ' ';
padder = ' ';
padder1 = repmat(padder,1,49);
padder2 = padder;
padder3 = repmat(padder,1,40);
padder4 = repmat(padder,1,6);
padder5 = repmat(padder,1,31);
padder6 = repmat(padder,1,21);
padderN = repmat(padder,1,4);

%%
strikeO = NaN(N,1);
dipO = strikeO;
rakeO = strikeO;
conf0 = strikeO;

for i = 1:N
    id = E(i).id;
    cd(baseDir);
    idDirDNE = ~exist(id,'dir');
    if idDirDNE
        mkdir(id);
    end
    cd(id);
    !rm *.out*
    outfile1 = sprintf("%s.phase",id);
    outfile2 = strcat("%s.stations",id);
    outfile3 = strcat("%s.reverse",id);

    phases = E(i).Pphases;
    polarities = pull(phases,'polarity');
    pI = polarities ~= 0;
    lPols = sum(pI);
    pols = "D";

    if lPols < 3
        fprintf("not enough polarities. cannot process: %d\n",i);
        continue;
    end

    phases = phases(pI);
    polarities = pull(phases,'polarity');
    pols = repmat(pols,lPols,1);
    dI = polarities < 0;
    pols(dI) = "U";
    inpFile = strcat(id,'.inp');

    disp([E(i).lat E(i).lon])
    [yyyy_,month_,day_,hour_,minute_,seconds_] = datevec(t(i));
    if lat(i) > 0
        latDeg = floor(lat(i));
        latMinute = abs(lat(i) - latDeg)*60;
    else
        latDeg = abs(ceil(lat(i)));
        latMinute = 60*abs(lat(i) - ceil(lat(i)));
        latSign = 'S';
    end
    lonDeg = abs(ceil(lon(i)));
    lonMinute = abs(lon(i) + lonDeg)*60;

    %%
    yrStr = num2str(yyyy_);
    sec2 = seconds_*100;

    %%
    mnStr = sprintf("%02d",month_);
    dyStr = sprintf("%02d",day_);
    hrStr = sprintf("%02d",hour_);
    minStr = sprintf("%02d",minute_);
    secStr = sprintf("%04d",floor(sec2));
    idStr = strcat(yrStr,mnStr,dyStr,hrStr,minStr,secStr);

    %%
    str1 = string(compose(formatSpec1,yyyy_,month_,day_,hour_,minute_,seconds_,...
        latDeg,latSign,latMinute,lonDeg,lonSign,lonMinute,19,padderN,lPols,...
        padder1,horizErr(i),padder2,verErr(i),padder3,mag(i),padder4,idStr));

    %%
    pStnms = string(pull(phases,'stnm'));
    pNet = pull(phases,'ntwk');
    pChan = pull(phases,'chan');
    pImp = repmat('I',lPols,1);
    [pStnmSort,sI] = sort(pStnms);
    pNetSort = pNet(sI);
    pChanSort = pChan(sI);
    for k = 1:length(pChanSort)
        p_ = pChanSort(k);
        if strcmp(p_,"SH?")
            pChanSort(k) = "SHZ";
        end
    end
    pols = pols(sI);

    %% create polarity file
    str2 = string(compose(formatSpec2,pStnmSort,pNetSort,pChanSort,pImp,pols));
    str2a = string(compose(formatSpec2a,idStr));
    strMaster = [string(str1); string(str2); string(str2a)];

    %%
    fileID = fopen(outfile1,"w");
    fprintf(fileID,"%s\n",strMaster);
    fclose(fileID);

    %% create station list (has to be in alphabetical order)
    [stlat,stlon,stelev] = metaDataFromStationList(pStnmSort,pNetSort,pChanSort);
    %table(pStnmSort,pNetSort,pChanSort,stlat,stlon,stelev)

    % optional, using velocity model that start at -2 km
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,zeros(size(stelev)),padder6,pNetSort);
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,-stelev,padder6,pNetSort);  %PVIL 97
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,stelev-2000,padder6,pNetSort); %PVIL 97
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,stelev,padder6,pNetSort);

    % optional, use velocity model that starts at zero
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,zeros(size(stelev)),padder6,pNetSort); %PVIL 105
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,-stelev,padder6,pNetSort);  %PVIL 105
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,stelev-2000,padder6,pNetSort); %PVIL 105
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,2000-stelev,padder6,pNetSort); %PVIL 105

    % optional, use velocity model that starts at -1.1
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,zeros(size(stelev)),padder6,pNetSort); %PVIL 101
    str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,round(stelev),padder6,pNetSort);  %PVIL 101
    %str3 = compose(formatSpec3,pStnmSort,pChanSort,padder5,stlat,stlon,stelev-2000,padder6,pNetSort); %PVIL 101

    %%these results indicate that HASH does not care about the station
    %%elevations when computing take off angles, azimuth

    %%
    fileID = fopen(outfile2,'w');
    fprintf(fileID,'%s\n',string(str3));
    fclose(fileID);

    %% create reveral file (in alphabetical order)
    str4 = compose(formatSpec4,pStnmSort,zeros(lPols,1),zeros(lPols,1));
    fileID = fopen(outfile3,'w');
    fprintf(fileID,'%s\n',string(str4));
    fclose(fileID);

    %%
    !ln -s ../vz.* .
    inpFile = char(inpFile);
    dlmwrite(inpFile,char(outfile2),'delimiter','');                % station list file
    dlmwrite(inpFile,char(outfile3),'-append','delimiter','');      % reversal file
    dlmwrite(inpFile,char(outfile1),'-append','delimiter','');      % input (observed polarities) file
    dlmwrite(inpFile,'test2.out','-append','delimiter','');         % output file name for focal mechanisms
    dlmwrite(inpFile,'test2.out2','-append','delimiter','');        % output file name for acceptable planes
    dlmwrite(inpFile,8,'-append','delimiter','');                   % mininum number of polarities (e.g., 8)
    dlmwrite(inpFile,350,'-append','delimiter','');                 % maximum azimuthal gap (e.g., 90)
    dlmwrite(inpFile,170,'-append','delimiter','');                 % maximum takeoff angle gap (e.g., 60)
    dlmwrite(inpFile,1,'-append','delimiter','');                   % grid angle for focal mech search, in degrees(min 5)
    dlmwrite(inpFile,99,'-append','delimiter','');                  % number of trials (e.g., 30)
    dlmwrite(inpFile,500,'-append','delimiter','');                 % maxout for focal mech. output (e.g., 500)
    dlmwrite(inpFile,0.1,'-append','delimiter','');                 % fraction of picks presumed bad (e.g., 0.10)
    dlmwrite(inpFile,1000,'-append','delimiter','');                % maximum allowed source-station distance,in km (e.g., 120)
    dlmwrite(inpFile,45,'-append','delimiter','');                  % angle for computing mechanisms probability, in degrees (e.g., 45)
    dlmwrite(inpFile,0.01,'-append','delimiter','');                % probability threshold for multiples (e.g., 0.1)
    dlmwrite(inpFile,7,'-append','delimiter','');                   % number of velocity models (max 10)
    dlmwrite(inpFile,'vz.sierranegra3','-append','delimiter','');
    dlmwrite(inpFile,'vz.vb1','-append','delimiter','');
    dlmwrite(inpFile,'vz.quito','-append','delimiter','');
    dlmwrite(inpFile,'vz.socal','-append','delimiter','');
    dlmwrite(inpFile,'vz.north','-append','delimiter','');
    dlmwrite(inpFile,'vz.lab1','-append','delimiter','');
    dlmwrite(inpFile,'vz.sgm1','-append','delimiter','');
    %dlmwrite(inpFile,' ','-append','delimiter',''); %<-- experimental, 13 JUN 2022
    %! ~/soft/HASH_v1.2/hash_driver2 < igepn2019kjzn.inp ;

    runCmd = sprintf("~/soft/HASH_v1.2/hash_driver2 < %s",inpFile);
    disp(runCmd);
    fprintf('running command: %s\n',runCmd);
    unix(runCmd);
    !wc test2.out*;
    !cat test2.out

    %cd ../
    %% here write a little routine that determines whether a solution was found, then extract that solution and plot it
    %mt = sdr2mt(27,82,133); %338,66,89); %25,74,135); %142,29,66); %31,73,134); %142,28,67);
    %focalmech(mt,0,0,1); axis equal;

    %%
    pwd
    fid = fopen('test2.out','r');
    %5            %10            %15            %20            %25         %29
    C1 = textscan(fid,'%d %d %d %d %d %d %f %s %f %s %f %f %f %s %f %f %s %f %d %d %d %d %d %d %d %d %d %d %s %d %d %s');

    fclose(fid);
    lc = length(C1{1});
    if lc
        strike = double(C1{22});
        dip = double(C1{23});
        rake = double(C1{24});
        percentLikely = C1{30};
        conf0 = max(percentLikely);
        bI = find(percentLikely == conf0);
        bI = bI(1);

        for j = bI %1:lc
            close all;
            strike_ = strike(j);
            dip_ = dip(j);
            rake_ = rake(j);
            strikeO(i) = strike_;
            dipO(i) = dip_;
            rakeO(i) = rake_;
            mt = sdr2mt(strike_,dip_,rake_);
            figure('units','normalized','outerposition',[0 0 1 1]);
            focalmech(mt,0,0,1);
            axis equal;
            ax = gca;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            if plotFlag
                close all;
                %%
                ax = plotEvents(E(i),true);
                hold(ax,'on')
                ax(2) = axes;
                focalmech(mt,E(i).lon,E(i).lat,0.01);
                %focalmech(ax(2),mt,E(i).lon,E(i).lat,0.01);
                %focalmech(mt,-79.76+mtShifter,1.07,0.015*sqrt(mlv));


                %
                ax(2).Visible = 'off';
                axis(ax(2),'equal');
                linkaxes(ax,'xy');
                axis(ax,[E(i).lon - 0.1,E(i).lon + 0.1,E(i).lat - 0.1,E(i).lat + 0.1]);

                %%
                cd(baseDir);
                cd(id)
                fname = strcat("mechanism_",num2str(1));
                title(ax(1),{strcat(datestr(t),", mag: ",num2str(mag),...
                    ", lat: ",num2str(lat),", lon: ",num2str(lon),...
                    ", depth: ",num2str(depth));strcat(num2str(strike_),...
                    ", ",num2str(dip_),", ",num2str(rake_))},'interpreter','latex');
                if printFlag
                    print('-djpeg',fname)
                end
            end
        end
    end
    cd(baseDir);
end