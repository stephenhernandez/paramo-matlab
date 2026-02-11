%function [tnew,pks,pks2,nMatch] = PlataRepeaterSearch_v1()
clear; close all; %clc;
warning off signal:findpeaks:largeMinPeakHeight

cd ~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/

kstnms = ["BREF";"BNAS";"BVC2";"BTAM";"BMOR";"BMAS";"BPAT";"BULB";"BRUN";"BBIL"];
lkstnms = length(kstnms);
kcmpnms = ["BHZ";"BHN";"BHE"];
lcmpnms = length(kcmpnms);

lfc = 1;
hfc = 4;
newFs = 16;
%tw = 0.002;
secDur = 100;
tw = round(0.5*newFs*secDur);
minTraces = 8;
verboseFlag = false;
weightedFlag = true;
maxLagSecs = 10;

%tStart = datetime(2006,01,180);
%tStart = datetime(2009,11,01);
%tEnd = datetime(2022,05,20);
%tStart = datetime(2016,05,05);
%tEnd = datetime(2016,05,05);
%tStart = datetime(2016,04,23);
%tEnd = datetime(2016,04,23);
tStart = datetime(2016,05,07);
tEnd = datetime(2016,05,07);

dayVec = (tStart:tEnd)';
ldays = length(dayVec);

close all;
families = dir('fam*.txt');
lfams = length(families);

%
popFlag = true;
if popFlag %do I populate waveform structs?
    tic;
    for i = 1:lfams
        famTimes = importdata(families(i).name);
        clear t;
        rowHeaders = famTimes.rowheaders;
        mins_ = famTimes.data(:,1);
        secs_ = famTimes.data(:,2);
        nTimes = length(rowHeaders);
        t = NaT(nTimes,1);

        for j = 1:nTimes
            t_ = strcat(string(rowHeaders(j)),":",string(mins_(j)),":",string(secs_(j)));
            try
                t(j) = dn2dt(datenum(char(t_),'yyyy-mm-ddTHH:MM:SS.FFF'));
            catch
                t(j) = dn2dt(datenum(char(t_),'yyyy-mm-ddTHH:MM:SS'));
            end
        end
        t = sort(t);

        S = extractWaveforms(t+seconds(0),seconds(secDur),kstnms,kcmpnms,"EC","",true,verboseFlag,2,true,[lfc,hfc,false,false]);
        sizeS = size(S);
        S = resampleWaveforms(detrendWaveforms(S),newFs);
        S = reshape(S,sizeS);

        natters = reshape(pull(S,'ref'),sizeS);
        goodTemplateSearchIndex = find(sum(~isnat(natters),2) > 0)';
        if isempty(goodTemplateSearchIndex)
            fprintf('insufficient sandro times to extract a stacked template, no Templates!\n');
            continue;
        end

        npts = pull(S,'npts');
        npts(~isfinite(npts)) = [];
        npts_ = max(unique(npts));

        lgood = length(goodTemplateSearchIndex);

        ddd = zeros(npts_*lcmpnms*lkstnms,lgood);
        SIMain = false(size(ddd));

        j_ = 1;
        for j = goodTemplateSearchIndex
            S__ = S(j,:)';
            S__ = double(pull(nanGapWaveforms(detrendWaveforms(nanGapWaveforms(S__,NaN)),0)));
            SI = S__ == 0 | ~isfinite(S__);
            S__(SI) = 0;
            S__ = normalizeWaveforms(S__);
            SI = S__ == 0;
            S__(SI) = NaN;
            SI = SI(:);
            SIMain(:,j_) = SI;
            ddd(:,j_) = S__(:);
            j_ = j_ + 1;
        end

        SIMainOrig = logical(SIMain);
        ddd(SIMainOrig) = 0;
        ddd = normalizeWaveforms(ddd);
        dddOrig = ddd;

        % at this point, for family 1, i pruned 2 of the 15 events (now 13)
        [maxccp,plags] = doccFreqCircShift(dddOrig,true,[],[],maxLagSecs*newFs);
        [shifted_data,~,~,~,raw_shifts] = apply_vdcc(ddd,[maxccp plags],weightedFlag);
        SIMain = apply_shifts(SIMainOrig,raw_shifts);
        SIMain = logical(SIMain);
        shifted_data(SIMain) = NaN;

        shifted_data = median(shifted_data,2,'omitnan');
        shifted_data(~isfinite(shifted_data)) = 0;
        shifted_data = reshape(shifted_data,size(S__));
        snrs = zeros(nTimes,1);
        for j = 1:nTimes
            S_ = S(j,:)';
            badI = isnat(pull(S_,'ref'));
            if sum(badI) < lcmpnms*lkstnms
                S_(badI) = [];
                ddumb = double(pull(S_));
                snrs(j) = sum(peak2rms(ddumb),'omitnan');
            end
        end

        [maxSnr,maxI] = max(snrs);
        if maxSnr == 0
            fprintf('no suitable templates\n');
            continue;
        end

        fprintf('max snr is %f and was found at index %d and time %s\n',...
            maxSnr,maxI,datestr(t(maxI)));
        S_ = S(sum(~isnat(natters),2) == max(sum(~isnat(natters),2)),:);
        S_ = S_(1,:)';
        referenceTime = S_(1).ref;

        j_ = 0;
        lS_ = length(S_);
        T = populateWaveforms(lS_);
        for j = 1:lS_
            shiftedData = shifted_data(:,j);
            snrTmp = rms(shiftedData)./mad(shiftedData,1);
            if rssq(shiftedData) > 0 && snrTmp >= 1.7
                j_ = j_ + 1;
                T(j_) = dealHeader(S_(j),shifted_data(:,j),newFs,referenceTime);
            end
        end
        T = syncWaveforms(T(1:j_));
        if i < 10
            save(strcat('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/optimal_templates_family',['0' num2str(i)]),'T');
        else
            save(strcat('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/optimal_templates_family',num2str(i)),'T');
        end
    end
    toc;
end

%%
close all;
thresh = 8;
maxN = 1e4;

famFiles = dir('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/optimal_templates_family*.mat');
lfams = length(famFiles);

tnew = NaT(maxN,lfams);
rawCC = NaN(maxN,lfams);
madScore = rawCC;
nMatch = zeros(1,lfams);

%loop through days
for j = 1:ldays
    tic;
    day_ = dayVec(j);
    fprintf('<strong>%s</strong>\n',datestr(day_));
    L2 = loadWaveforms(day_,1,kstnms,kcmpnms,"EC","",true,verboseFlag); %load long time series
    toc;

    badI = isnat(pull(L2,'ref'));
    if ~sum(~badI)
        fprintf('no data for day: %s\n',datestr(day_));
        continue;
    end

    L2(badI) = [];
    L2 = syncWaveforms(...
        detrendWaveforms(...
        filterWaveforms(...
        detrendWaveforms(...
        intWaveforms(...
        detrendWaveforms(...
        nanGapWaveforms(...
        resampleWaveforms(...
        detrendWaveforms(...
        differentiateWaveforms(L2)),newFs),NaN)))),lfc,hfc)));

    sncl2 = repmat("",length(L2),1);
    for ii = 1:length(L2)
        sncl2(ii) = strcat(L2(ii).knetwk,L2(ii).kstnm,L2(ii).khole,L2(ii).kcmpnm);
    end

    tLong = getTimeVec(L2);
    L2Orig = L2;
    L2 = double(pull(L2));

    [~,ltraceL] = size(L2);
    if ltraceL < minTraces
        fprintf('not enough traces for day: %s\n',datestr(day_));
        continue;
    end

    for i = 1:ltraceL
        L2Orig(i) = dealHeader(L2Orig(i),L2(:,i),newFs,L2Orig(1).ref);
    end
    L2 = nanGapWaveforms(L2Orig,0);
    L2 = double(pull(L2)); %these traces will be recycled for each multiplet family

    tmpFlag = false;
    if tmpFlag
        axMain = [];
    end
    fprintf("done loading all data for day: <strong>%s</strong>\n",datestr(day_));
    toc;
    for i = 1:lfams
        familyFile = famFiles(i).name;
        if ~exist(familyFile,'file')
            fprintf('file for family %d does not exist\n',i);
            continue;
        end
        load(familyFile,'T');

        %synch S_ and L2
        Torig = T;
        sncl1 = repmat("",length(Torig),1);
        for ii = 1:length(Torig)
            sncl1(ii) = strcat(Torig(ii).knetwk,Torig(ii).kstnm,...
                Torig(ii).khole,Torig(ii).kcmpnm);
        end

        [lia,locb] = ismember(sncl1,sncl2);
        if ~sum(lia) || sum(lia) < minTraces
            fprintf('not enough traces in common\n');
            continue;
        end

        L = L2(:,locb(lia)); %use only L-traces corresponding to THIS multiplet
        [~,ltraceL] = size(L);
        if ltraceL < minTraces
            fprintf('not enough traces for day: %s\n',datestr(day_));
            continue;
        end

        T = syncWaveforms(Torig(lia));
        T = double(pull(T));
        [winlen,ltraceT] = size(T);
        T = normalizeWaveforms(flipud(T));

        box1 = ones(winlen,1);
        normers = sqrt(abs(fftfilt(box1,L.^2)));
        badNormersI = normers < 1 | ~isfinite(normers);
        normers(badNormersI) = NaN;

        cc = fftfilt(T,L);
        ccnormOrig = cc./normers;
        ccnorm = mean(ccnormOrig,2,'omitnan');

        missingSNCLs = sum(badNormersI,2);
        uniqMissingSNCLs = unique(missingSNCLs)';
        weights = ones(ltraceT+1,1);
        for ii = uniqMissingSNCLs
            wi = missingSNCLs == ii;
            we = 1./mad(ccnorm(wi),1);
            if isfinite(we)
                weights(ii+1) = we;
                ccnorm(wi) = we.*ccnorm(wi);
            end
        end

        mad_ = mad(ccnorm,1);
        [pks_,locs_] = findpeaks(ccnorm,'MINPEAKDISTANCE',0.2*winlen,'MINPEAKHEIGHT',thresh*mad_);

        locs_ = locs_-winlen+1;
        locsI = locs_ > 0;
        pks_ = pks_(locsI);
        locs_ = locs_(locsI);
        if ~isfinite(pks_)
            fprintf('no events found: %s\n',datestr(day_));
            continue;
        end

        ll = length(pks_);
        if ~ll
            fprintf('no events found on day: %s\n',datestr(day_));
            continue;
        end

        tnew(nMatch(i)+1:nMatch(i)+ll,i) = tLong(locs_);
        rawCC(nMatch(i)+1:nMatch(i)+ll,i) = pks_./weights(missingSNCLs(locs_+winlen-1)+1);
        madScore(nMatch(i)+1:nMatch(i)+ll,i) = pks_;
        nMatch(i) = nMatch(i) + ll;
        fprintf('family: <strong>%d</strong>,number of events found: <strong>%d</strong>, on day: <strong>%s</strong>, with <strong>%d</strong> cumulative events so far\n',i,ll,datestr(day_),nMatch(i));
        %save('~/research/now/plata_repeaters/SANDRO/OTAV_EV/temporaryVariableStore2010','tnew','pks','pks2','nMatch');

        if tmpFlag
            figure('units','normalized','outerposition',[0 0 1 1]);
            plot(tLong,ccnorm); zoom on;
            ax_ = gca;
            hold on;
            plot(ax_,tLong(locs_+winlen-1),pks_,'p');
            title(familyFile);
            axMain = [axMain; ax_];
            linkaxes(axMain,'x');
            fprintf("done processing family: %d\n",i);
        end
    end
    toc;
end
