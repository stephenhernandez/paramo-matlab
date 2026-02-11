function updateMongoRSAM(N)
if nargin < 1
    N = 1;
end
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
% clear; close all;
for n = 1:N
    server = "127.0.0.1";
    port = 27017;
    dbname = "IGEPN_19AUG2022";
    conn = mongoc(server,port,dbname);
    verboseFlag = true;

    % filter bank
    lfcs = [1/50; 0.2; 0.25; 0.3; 0.5; 0.6;  1; 2;  5;   10];
    hfcs = [1/20; 0.8;    2; 0.5;   4; 1.2; 16; 5; 15; -inf];
    meanFlags = [true; true; false; false];
    rmsFlags = [true; false; true; false];

    llfc = length(lfcs);
    units = 'vel';
    scaleFactor = 1e9;
    rsamDur = 60;
    npoles = 4;
    EPOCHREF = datetime(1970,01,01);

    %
    writeToDB = true;
    twMax = 0.008;
    MINDUR = 4;         %minutes
    newFs = 50;
    minFs = 40;
    horFlag = false;

    nowTime = dn2dt(now) + hours(5);
    currentDay = dateshift(nowTime,"start","day");
    if ~isopen(conn)
        pause(10);
        conn = mongoc(server,port,dbname);
        pause(10);
        fprintf("Connection wasnt open, had to try again: %s\n",datestr(nowTime));
    end

    collectionNames = conn.CollectionNames;
    fileList = existsToday(false,true);
    lTodayFiles = length(fileList);
    for kk = 1:lTodayFiles
        thisFile = fileList(kk);
        splits = split(thisFile,".");
        knetwk = splits(1);
        kstnm = splits(2);
        khole = splits(3);
        kcmpnm = splits(4);

        sncl_ = strcat(knetwk,kstnm,khole,kcmpnm);
        SNCLCollectionName1 = strcat("RSAM_",sncl_);

        lia = ismember(SNCLCollectionName1,collectionNames);
        sumlia = sum(lia);
        if ~sumlia
            fprintf('%s as collection does not exist, skipping....\n',...
                SNCLCollectionName1);
            continue;
        end

        nDocs = count(conn,SNCLCollectionName1);
        if ~nDocs
            fprintf('%s as collection does not have records, skipping....\n',...
                SNCLCollectionName1);
            continue;
        end

        lastRow = find(conn,SNCLCollectionName1,'Skip',nDocs-1,'Projection','{"t":1.0}');
        lastRecordedObservation = lastRow.t;
        if isempty(lastRecordedObservation)
            fprintf('last row of %s is empty, skipping....\n',...
                SNCLCollectionName1);
            continue;
        end
        lastRecordedObservation = EPOCHREF + milliseconds(lastRecordedObservation);
        latency = ceil(days(nowTime - lastRecordedObservation));
        if latency > 10
            fprintf('%s: Last update too long ago, memory issues, skipping....\n',...
                SNCLCollectionName1);
            continue;
        end

        cutStart = max([dateshift(lastRecordedObservation,'start','day') ...
            lastRecordedObservation - minutes(5)]);
        S = extractWaveforms(cutStart,floor(seconds(nowTime-cutStart)),kstnm,kcmpnm,knetwk,khole);

        refs = pull(S,'ref');
        badRefs = isnat(refs);
        if ~any(~badRefs) %if none are good...
            fprintf('%s: no data, skipping....\n',...
                SNCLCollectionName1);
            continue;
        end

        nowTime = dn2dt(now) + hours(5);
        currentDay = dateshift(nowTime,"start","day");
        if ~isopen(conn)
            pause(10);
            conn = mongoc(server,port,dbname);
            pause(10);
            fprintf("Connection wasnt open, had to try again: %s\n",datestr(nowTime));
            continue;
        end

        %%
        if ~exist('S','var')
            continue;
        end

        %%
        Supdated = resampleWaveforms(...
            detrendWaveforms(...
            differentiateWaveforms(S)),newFs); %prefilt

        Supdated = filterWaveforms(Supdated,min(0.5*lfcs));
        Supdated = nanGapWaveforms(Supdated,0);

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

        %%
        nWritten = 0;
        for i = 1:lUpdatedSNCLs
            sncl_ = updatedSNCLs(i);
            Fs_ = Fs(i);
            if Fs_ < minFs %&& ( hfc >= 10 ||  lfc >= 10 )
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
            tw = min([twMax 720/npts]);
            dOrig = taper(detrend(dOrig),tw);
            Dorig = fft(dOrig);     % frequency-domain

            SNCLCollectionName1 = strcat("RSAM_",sncl_);
            lia = ismember(SNCLCollectionName1,collectionNames);
            sumlia = sum(lia);
            if ~sumlia
                if verboseFlag
                    fprintf('%s as collection does not exist, creating....\n',...
                        SNCLCollectionName1);
                end

                try
                    createCollection(conn,SNCLCollectionName1);
                    collectionNames = conn.CollectionNames;
                catch ME
                    warning(ME);
                    continue;
                end
            end

            nDocs = count(conn,SNCLCollectionName1);
            if ~nDocs
                lastRecordedObservation = EPOCHREF;
            else
                lastRow = find(conn,SNCLCollectionName1,'Skip',nDocs-1,'Projection','{"t":1.0}');
                lastRecordedObservation = lastRow.t;
                if isempty(lastRecordedObservation)
                    if nDocs > 1
                        %% DANGEROUS! delete last row, drop collection, create anew, reinsert everything leftover...
                        collAll = find(conn,SNCLCollectionName1);
                        collAll(end) = [];
                        collAll = rmfield(collAll,'_id');
                        dropCollection(conn,SNCLCollectionName1);
                        createCollection(conn,SNCLCollectionName1);
                        nReinserted = insert(conn,SNCLCollectionName1,struct2table(collAll));
                        if verboseFlag
                            fprintf(1,'error detected, collection dropped, %d rows RE-inserted...\n',nReinserted);
                        end
                        continue;
                    else
                        fprintf(2,"MAJOR ERROR\n");
                        close(conn);
                        return;
                    end
                end
                lastRecordedObservation = EPOCHREF + milliseconds(lastRecordedObservation);
            end

            %%


            %% decide whether enough new data to proceed
            lengthNew = minutes(tEnd_ - lastRecordedObservation); %assume no gaps
            if lengthNew < MINDUR
                continue;
            end

            %% Unique RSAM struct to this SNCL
            RSAMstruct = populateRSAM(lfcs,hfcs,meanFlags,rmsFlags);

            %% loop through all filters
            for j = 1:llfc
                lfc = lfcs(j);
                hfc = hfcs(j);

                %% apply filters
                if ~isnat(Tresp) && ~strcmp(kcmpnm_,"BDF") && ...
                        ~strcmp(kcmpnm_,"SHZ")
                    Hbu = freqOperator(npts_,lfc,hfc,Fs_,npoles);
                    D = Dorig.*Hresp.*Hbu;
                    df = scaleFactor*ifft(D,'symmetric');
                else
                    if verboseFlag
                        fprintf('no response info for: %s, proceed to use original data...\n',...
                            sncl_);
                    end
                    Hd = zpkOperator(lfc,hfc,Fs_,npoles);
                    df = filter(Hd,dOrig);
                end
                df = detrend(cumsum(df)*delta_);

                %% RSAM Section
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
                if nNewRecords < 2
                    if verboseFlag
                        fprintf('%s: not enough data to insert new records, continuing...\n',sncl_);
                    end
                    continue;
                end

                tVec = milliseconds(tVec(tI) - EPOCHREF);
                RSAMstruct.t = tVec;
                ampVec = ampVec(tI);
                RSAMstruct.(tmpFieldName) = single(ampVec);

                lFlags = length(meanFlags);
                for k = 2:lFlags
                    meanFlag_ = meanFlags(k);
                    rmsFlag_ = rmsFlags(k);
                    ampVec = amplitudeVector(df,Fs_,rsamDur,meanFlag_,rmsFlag_);
                    ampVec = ampVec(iStart:winlen:end); % at one point i forgot this line...
                    technique = rsamTechniqueStr(meanFlag_,rmsFlag_);
                    tmpFieldName = rsamFieldNames(lfc,hfc,technique);
                    ampVec = ampVec(tI);
                    RSAMstruct.(tmpFieldName) = single(ampVec);
                end
            end

            %%
            if writeToDB
                try
                    nWritten = nWritten + 1;
                    nInserted = insert(conn,SNCLCollectionName1,struct2table(RSAMstruct));
                    if verboseFlag
                        fprintf('%d new records inserted into %s\n',nInserted,sncl_);
                    end
                catch ME
                    warning(ME.message);
                    continue;
                end
            end
        end
    end
    close(conn);
end
