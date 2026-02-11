%% Common commands
% n = count(conn,collectionName);
% n = insert(conn,collectionName,documents)
% n = remove(conn,collectionName,mongoquery)
% n = update(conn,collectionName,findquery,updatequery)

% createCollection(conn,collectionName);
% dropCollection(conn,collectionName);
% totaldocs2 = find(conn,collectionName,'Skip',N-1,'Projection','{"t":1.0}');

% mongoquery = '{"t":{$gte:1609459200000}}';
% totaldocs = count(conn,collectionName,'Query',mongoquery);

%% Potential collections to create and populate
% SNCLCollectionName1 = strcat("RSAM_",sncl_);
% SNCLCollectionName2 = strcat("HEALTH_",sncl_);
% collectionName3 = strcat("NOISE_",sncl_);
% collectionName4 = strcat("ITERATION_STATISTICS_",sncl_);

%% TO ADD:
% global catalog from usgs
% local catalog from sc3

%%
clear; close all;
server = "127.0.0.1";
port = 27017;
%dbname = "IGEPN_19AUG2022";
%dbname = "IGEPN_23NOV2022";
%dbname = "IGEPN_10APR2023";
dbname = "IGEPN_07JAN2025";
conn = mongoc(server,port,dbname);

%% filter bank
lfcs = [1/50; 0.2; 0.25; 0.3; 0.5; 0.6;  1; 2;  5;   10];
hfcs = [1/20; 0.8;    2; 0.5;   4; 1.2; 16; 5; 15; -inf];
meanFlags = [true; true; false; false];
rmsFlags = [true; false; true; false];

llfc = length(lfcs);
lFlags = length(meanFlags);
MAXCOMBOS = llfc*lFlags;
MAXINSERTSIZE = 10;
HORFLAG = true;
units = 'vel';
scaleFactor = 1e9;
rsamDur = 60;
npoles = 4;
epochRef = datetime(1970,01,01);

%
writeToDB = true;
tw = 0.008;
secDur = 1200;       %seconds
MINDUR = 4;         %minutes
newFs = 50;
minFs = 40;
horFlag = true;
verboseFlag = true;

%%
%S = seedData(secDur,horFlag);
nowTime = datetime("now")+hours(5); %dn2dt(now) + hours(5);
currentDay = dateshift(nowTime,"start","day");
S = seedData2(currentDay,nowTime-seconds(2*secDur),nowTime,HORFLAG,verboseFlag);
refs = pull(S,'ref');
badRefs = isnat(refs);
if ~any(~badRefs) %if none are good...
    S = seedData(secDur,true);
end

nIter = 0;
iterationCollection = "ITERATION_STATISTICS";
SAG1InfrasoundCollection1 = "INFRASOUND_SAG1_MAIN";
SAG1InfrasoundCollection2 = "INFRASOUND_SAG1_LAST";

%%
SAG1InfrasoundArray = ["ECSAG101BDF";...
    "ECSAG102BDF";...
    "ECSAG103BDF";...
    "ECSAG104BDF";...
    "ECSAG105BDF"];

%%
while nIter < 5
    iterationStartTic = tic;
    nowTime = datetime("now")+hours(5); %dn2dt(now) + hours(5);
    currentDay = dateshift(nowTime,"start","day");
    if ~isopen(conn)
        fprintf("Connection not open, trying again: %s\n",nowTime);
        pause(5);
        conn = mongoc(server,port,dbname);
        continue;
    end

    %
    collectionNames = conn.CollectionNames;

    %
    if ~exist('S','var')
        S = seedData2(currentDay,nowTime-seconds(secDur),nowTime,HORFLAG,verboseFlag);
        %S = seedData(secDur,horFlag,MINDUR);
        continue;
    elseif length(S) < 2
        S = seedData2(currentDay,nowTime-seconds(secDur),nowTime,HORFLAG,verboseFlag);
        %S = seedData(secDur,horFlag,MINDUR);
        continue;
    end

    %
    if nIter == 0
        Supdated = S;
        deltas = pull(Supdated,'delta');
        dI = deltas > 1/minFs;
        Supdated(dI) = [];
        nUpdated = length(Supdated);
    else
        [S,Supdated] = updateWaveforms(S,secDur,horFlag,MINDUR);
        deltas = pull(Supdated,'delta');
        dI = deltas > 1/minFs;
        Supdated(dI) = [];
        nUpdated = length(Supdated);
    end

    if nUpdated < 2
        iterationEndToc = toc(iterationStartTic);
        if iterationEndToc < 40
            fprintf('waiting a bit..., be patient\n');
            pause(20)
        end
        continue;
    end

    %
    Supdated = resampleWaveforms(...
        detrendWaveforms(...
        differentiateWaveforms(Supdated)),newFs);

    e = pull(Supdated,'e');
    refs = pull(Supdated,'ref');
    npts = pull(Supdated,'npts');
    deltas = pull(Supdated,'delta');
    knetwks = pull(Supdated,'knetwk');
    kstnms = pull(Supdated,'kstnm');
    kholes = pull(Supdated,'khole');
    kcmpnms = pull(Supdated,'kcmpnm');

    tEnd = refs + e;

    updatedSNCLs = strcat(knetwks,kstnms,kholes,kcmpnms);
    lUpdatedSNCLs = length(updatedSNCLs);
    if ~lUpdatedSNCLs
        continue;
    end

    Fs = 1./deltas;
    FsI = Fs >= 1;
    Fs(FsI) = round(Fs(FsI));

    mySNCLs = [kstnms kcmpnms knetwks kholes];
    responseStructure = singleSNCLFreqResponse(mySNCLs,currentDay,currentDay+1,npts,Fs,units);
    uniqFs = unique(Fs);

    %read info related to current and previous iterations
    collectionNames = createNewCollection(conn,iterationCollection,collectionNames,verboseFlag);
    cumulativeIterations = count(conn,iterationCollection);
    cumulativeIterations = cumulativeIterations + 1;

    %
    nWritten = 0;
    for i = 1:lUpdatedSNCLs
        sncl_ = updatedSNCLs(i);
        Fs_ = Fs(i);
        if Fs_ < minFs
            fprintf('%s has too low of a sample rate, continuing...\n',sncl_);
            continue;
        end
        e_ = minutes(e(i));
        if e_ < MINDUR
            fprintf('%s duration too short, continuing...\n',sncl_);
            continue;
        end

        npts_ = npts(i);
        S_ = Supdated(i);
        kcmpnm_ = kcmpnms(i);
        R = responseStructure(i);

        ref_ = refs(i);
        tEnd_ = tEnd(i);

        delta_ = 1/Fs_;
        Tresp = R.Tstart;
        Hresp = R.H;
        dOrig = S_.d;           % time-domain
        dOrig = taper(detrend(dOrig),tw);
        Dorig = fft(dOrig);     % frequency-domain

        SNCLCollectionName1 = strcat("RSAM_",sncl_);
        collectionNames = createNewCollection(conn,SNCLCollectionName1,collectionNames,verboseFlag);

        nDocs = count(conn,SNCLCollectionName1);
        if ~nDocs
            lastRecordedObservation = epochRef;
        else
            lastRow = find(conn,SNCLCollectionName1,'Skip',nDocs-1,'Projection','{"t":1.0}');
            lastRecordedObservation = lastRow.t;
            if isempty(lastRecordedObservation)
                if nDocs > 1
                    % DANGEROUS! delete last row, drop collection, create anew, reinsert everything leftover...
                    collAll = find(conn,SNCLCollectionName1);
                    collAll(end) = [];
                    collAll = rmfield(collAll,'_id');
                    dropCollection(conn,SNCLCollectionName1); %<-- very dangerous
                    createCollection(conn,SNCLCollectionName1);
                    nReinserted = insert(conn,SNCLCollectionName1,struct2table(collAll));
                    if verboseFlag
                        fprintf(1,'error detected, collection dropped, %d rows RE-inserted...\n',nReinserted);
                    end
                    continue;
                else
                    fprintf(2,"MAJOR ERROR\n");
                    return;
                end
            end
            lastRecordedObservation = epochRef + milliseconds(lastRecordedObservation); %convert to matlab datetime
        end


        % decide whether enough new data to proceed
        lengthNew = minutes(tEnd_ - lastRecordedObservation); %assume no gaps
        if lengthNew < MINDUR
            continue;
        end

        % loop through all filters
        nFilterCombos = 0;
        for j = 1:llfc
            lfc = lfcs(j);
            hfc = hfcs(j);

            % apply filters
            if ~isnat(Tresp) && ~strcmp(kcmpnm_,"BDF") && ...
                    ~strcmp(kcmpnm_,"SHZ")
                Hbu = freqOperator(npts_,lfc,hfc,Fs_,npoles);
                D = Dorig.*Hresp.*Hbu;
                df = scaleFactor*ifft(D,'symmetric');
            else
                if verboseFlag
                    %fprintf('no response info for: %s, proceed to use original data...\n',sncl_);
                end
                Hd = zpkOperator(lfc,hfc,Fs_,npoles);
                df = filter(Hd,dOrig);
            end
            df = detrend(cumsum(df)*delta_);

            % RSAM Section
            k = 1;
            meanFlag_ = meanFlags(k);
            rmsFlag_ = rmsFlags(k);
            technique = rsamTechniqueStr(meanFlag_,rmsFlag_);
            tmpFieldName = rsamFieldNames(lfc,hfc,technique);

            [ampVec,winlen] = amplitudeVector(df,Fs_,rsamDur,meanFlag_,rmsFlag_);
            newRef = dateshift(ref_,'end','minute');
            iStart = t2i(newRef,ref_,delta_);
            ampVec = ampVec(iStart:winlen:end);
            lAmpVec = length(ampVec);
            tVec = newRef + seconds(rsamDur*(0:lAmpVec-1)');
            tI = tVec > lastRecordedObservation;
            tI(end) = false; %in theory, will help prevent data modified by taper from entering database
            nNewRecords = sum(tI);
            if nNewRecords < 1
                if verboseFlag
                    fprintf('%s: not enough data to insert new records, continuing...\n',sncl_);
                end
                continue;
            end

            ampVec = ampVec(tI);
            tVec = tVec(tI);

            latency_ = minutes(nowTime - tVec);
            tVec = milliseconds(tVec - epochRef); %convert to javascript format
            if j == 1
                RSAM = NaN(nNewRecords,MAXCOMBOS);
            end
            nFilterCombos = nFilterCombos + 1;
            RSAM(:,nFilterCombos) = ampVec;

            for k = 2:lFlags
                meanFlag_ = meanFlags(k);
                rmsFlag_ = rmsFlags(k);
                ampVec = amplitudeVector(df,Fs_,rsamDur,meanFlag_,rmsFlag_);
                ampVec = ampVec(iStart:winlen:end); % at one point i forgot this line...
                technique = rsamTechniqueStr(meanFlag_,rmsFlag_);
                tmpFieldName = rsamFieldNames(lfc,hfc,technique);
                ampVec = ampVec(tI);
                nFilterCombos = nFilterCombos + 1;
                RSAM(:,nFilterCombos) = ampVec;
            end
        end

        %
        if writeToDB
            try
                if nNewRecords
                    for jj = 1:nNewRecords
                        t = tVec(jj); %in milliseconds since 1970...
                        latency = latency_(jj);
                        rsam = RSAM(jj,:);
                        nInserted = insert(conn,...
                            SNCLCollectionName1,table(latency,t,rsam));
                    end
                    if verboseFlag
                        fprintf('%d new records inserted into %s\n',nNewRecords,sncl_);
                    end
                end
                nWritten = nWritten + 1;
            catch ME
                fprintf(2,'couldnt write: %s\n',sncl_);
                warning(ME.message);
                nIter = nIter+1;
                iterationEndToc = toc(iterationStartTic);
                fprintf('Iteration: <strong>%d</strong>, Duration: <strong>%f</strong>, Total SNCLs Written: <strong>%d</strong>, End Time: %s\n',...
                    nIter,iterationEndToc,nWritten,nowTime);

                dur = iterationEndToc;
                t = milliseconds(nowTime - epochRef);
                nInserted = insert(conn,iterationCollection,table(cumulativeIterations,t,nWritten,dur));
                if iterationEndToc < 20
                    pause(10)
                end
                continue;
            end
        end
    end

    nIter = nIter+1;
    iterationEndToc = toc(iterationStartTic);
    fprintf('Iteration: <strong>%d</strong>, Duration: <strong>%f</strong>, Total SNCLs Written: <strong>%d</strong>, End Time: %s\n',...
        nIter,iterationEndToc,nWritten,nowTime);

    if writeToDB
        dur = iterationEndToc;
        t = milliseconds(nowTime - epochRef);
        itsInserted = insert(conn,iterationCollection,table(cumulativeIterations,t,nWritten,dur));
    end

    knetwks = pull(S,'knetwk');
    kstnms = pull(S,'kstnm');
    kholes = pull(S,'khole');
    kcmpnms = pull(S,'kcmpnm');
    allMySNCLs = [kstnms kcmpnms knetwks kholes];
    allMySNCLs = strcat(allMySNCLs(:,3),allMySNCLs(:,1),allMySNCLs(:,4),allMySNCLs(:,2));

    % update SAG1 Infrasound Catalog
    collectionNames = createNewCollection(conn,SAG1InfrasoundCollection1,collectionNames,verboseFlag);
    collectionNames = createNewCollection(conn,SAG1InfrasoundCollection2,collectionNames,verboseFlag);
    updateSAG1InfrasoundCatalog(S,SAG1InfrasoundArray,allMySNCLs,updatedSNCLs,writeToDB,...
        conn,SAG1InfrasoundCollection1,SAG1InfrasoundCollection2);
    if iterationEndToc < 10
        pause(20);
    end

    % nowTime2 = datetime("now") + hours(5);
    % localDay = dateshift(nowTime2,'start','day');
    % crossedMidnight = nowTime2 - dateshift(nowTime,'start','day') > hours(0);
    % if crossedMidnight
    %     %%
    %     downloadTropomi();
    %     createTropomiFrames(true);
    %     createTropomiFrames(false);
    %     stackTropomi();
    %     close all;
    %     %%
    %     updateCotopaxiTremorCatalogBinaryClassifier();
    %
    %     %%
    %     dayStart = localDay - 1;
    %     dayEnd = localDay;
    %     dayInc = 1;
    %
    %     reventadorRepeaterSearch(dayStart,dayEnd,dayInc);
    %     batchChilesMarch2023(dayStart,dayEnd,dayInc);
    %     batchJob10(dayStart,dayEnd,dayInc); % guagua
    %     PlataRepeaterSearch_v2(dayStart,dayEnd,dayInc);
    %     CotoLPSearch(dayStart,dayEnd,dayInc);
    %     SAG1_infrasound_array_processing_v2(dayStart,dayEnd,dayInc);
    %     batchJob14(dayStart,dayEnd,dayInc); % coto vlps
    %     runSubspaceDetector(dayStart,dayEnd,dayInc);
    % end
end
close(conn);