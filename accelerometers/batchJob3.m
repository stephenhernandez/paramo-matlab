%function batchJob3(ids)
%if nargin < 1
%%
clear; close all; clc;
startup(true);
minMag = 3.8;
maxMag = 8;

% [evDescription,evStatus,t,eqlat,eqlon,eqdepth,eqmag,ids,stderr,azgap,phaseTot,nMLv,...
% timerr,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,magType,meanMethod,...
% locMethod,earthModel,creationTime,agencyID,evMode,scbullmag,authorID,evType] = readCat1();
try
    [~,~,t,~,~,~,eqmag,ids,stderr,azgap,nPhases,nMLv,...
        ~,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,~,~,...
        ~,~,~,~,~,scbullmag] = readCat1();
catch ME
    warning(ME.message);
    try
        [~,~,t,~,~,~,eqmag,ids,stderr,azgap,nPhases,~,...
            ~,eqlaterr,eqlonerr,eqdeptherr,eqmagerr,~,~,...
            ~,~,~,~,~,scbullmag] = readCat1();
    catch ME
        rethrow(ME);
    end
end

eqmag = round(eqmag*10)/10;
scbullmag = round(scbullmag*10)/10;

sError = sqrt(eqlaterr.^2 + eqlonerr.^2 + eqdeptherr.^2);
eI = eqmag <= maxMag & eqmag >= minMag & t > datetime(2011,01,01) & stderr > 0 & stderr < 2 & ...
    eqmagerr >= 0 & eqmagerr < 1 & sError <= 100 & azgap <= 360 & nMLv >= 4 & nPhases >= 4;

%%
eqmag = eqmag(eI);
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

[eqmag,eI] = sort(eqmag,"descend");
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

scbullmag(scbullmag < -100) = eqmag(scbullmag < -100);
eqmag = scbullmag;

eI = eqmag >= minMag & t > datetime(2011,01,01);

eqmag = eqmag(eI);
ids = ids(eI);

[~,eI] = sort(eqmag,"descend");
ids = ids(eI);
t = t(eI);
scbullmag = scbullmag(eI);

%%
npoles = 4;
hfc = 40;
mainDir = fullfile("~","masa","strong_motion/"); %tharp %'~/products/EcuadorStrongMotion/';
noiseWin = 2; %minutes
totalDur = 8; %minutes
for i = 1:length(ids)
    eventID_ = ids(i);
    yyyy = char(eventID_);
    yyyy = yyyy(6:9);
    files = dir(strcat(mainDir,yyyy,"/*/",eventID_,...
        "*/AmplitudeVsDistance_ACC.png"));
    lfiles = length(files);

    if lfiles > 0
        fprintf('already processed %s, mag: %f, skipping...\n',...
            eventID_,eqmag(i));
        continue;
    end
    fprintf('will process: <strong>%s</strong>\n',eventID_);
    if strcmp(mainDir,'~/masa/strong_motion/')
        try
            %writeStrongMotionRecords(eventID_,noiseWin,totalDur,1/40,hfc,npoles,'acc');
            writeStrongMotionRecordsNoPicks(eventID_,noiseWin,totalDur,1/40,hfc,npoles,'acc');
        catch
            fprintf(2,'something went wrong\n');
            continue;
        end

        try
            %writeStrongMotionRecords(eventID_,noiseWin,totalDur,1/10,hfc,npoles,'vel');
            writeStrongMotionRecordsNoPicks(eventID_,noiseWin,totalDur,1/10,hfc,npoles,'vel');
        catch
            fprintf(2,'something went wrong\n');
            continue;
        end

        try
            %writeStrongMotionRecords(eventID_,noiseWin,totalDur,1/5,hfc,npoles,'disp');
            writeStrongMotionRecordsNoPicks(eventID_,noiseWin,totalDur,1/5,hfc,npoles,'disp');
        catch
            fprintf(2,'something went wrong\n');
            continue;
        end

        try
            writeStrongMotionRecordsNoPicks(eventID_,noiseWin,totalDur,1/2,10,npoles,'WA');
        catch
            fprintf(2,'something went wrong\n');
            continue;
        end
    else
        try
            %writeStrongMotionRecords(eventID_,1,6,1/5,hfc,npoles,'disp');
            %writeStrongMotionRecordsNoPicks(eventID_,1,10,1/2,10,npoles,'WA');
        catch
            fprintf(2,'something went wrong\n');
            continue;
        end
    end
end