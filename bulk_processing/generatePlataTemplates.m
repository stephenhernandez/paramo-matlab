%function [tnew,pks,pks2,nMatch] = PlataRepeaterSearch_v1()
clear; close all; %clc;
warning off signal:findpeaks:largeMinPeakHeight

saveFlag = true;
cd ~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/
experimentFlag = true;
if experimentFlag
%     fid = fopen('targets_la_plata_v2.txt');
%     C1 = textscan(fid,'%d %f %f %f %d-%d-%d %d:%d:%f %f');
%     fclose(fid);
% 
%     eqmag = C1{11};
%     tTarget = datetime(C1{5},C1{6},C1{7},C1{8},C1{9},C1{10});
%     targetNumber = C1{1};
%     eqlon = C1{2};
%     eqlat = C1{3};
%     eqdepth = C1{4};
fid = fopen('allSandro.txt');
C1 = textscan(fid,'%4d %02d %02d %02d %02d %f %02d');
fclose(fid);

tTarget = datetime(cat(2,C1{1:6}));
targetNumber = C1{7};
else
    fid = fopen('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/fam08_full_OTAV.txt');
    C1 = textscan(fid,'%f-%f-%fT%f:%f:%f');
    fclose(fid);
    tTarget = datetime(cat(2,C1{1:6}));
    targetNumber = ones(length(tTarget),1)*8; %C1(:,9);
    %     outFile = "~/jica_plata_catalog_v2.txt";
    %     C1 = importdata(outFile);
    %     tTarget = datetime(C1(:,1:6));
    %     targetNumber = C1(:,9);
end

kstnms = ["BREF";"BNAS";"BVC2";"BTAM";"BMOR";"BMAS";"BPAT";"BULB";"BRUN";"BBIL"];
lkstnms = length(kstnms);
kcmpnms = ["BHZ";"BHN";"BHE"];
lcmpnms = length(kcmpnms);

lfc = 1;
hfc = 4;
newFs = 20;
secDur = 150;
tw = round(0.5*newFs*secDur);
minTraces = 6;
verboseFlag = false;
weightedFlag = true;
maxLagSecs = 2; % max lag allowed when trying to align waveforms

% load target_11;
% tTarget = t11;
% targetNumber = 11*ones(size(tTarget));
close all;
uniqTargets = unique(targetNumber); %length(families);
lTargets = length(uniqTargets);
transferFlag = true;
%%
tic;
for i = 1:lTargets                               
    targetNum_ = uniqTargets(i); %targetNumber(i);
    fI = targetNumber == targetNum_; %importdata(families(i).name);
    %     clear t;
    %     rowHeaders = famTimes.rowheaders;
    %     mins_ = famTimes.data(:,1);
    %     secs_ = famTimes.data(:,2);
    %     nTimes = length(rowHeaders);
    %     t = NaT(nTimes,1);
    %
    %     for j = 1:nTimes
    %         t_ = strcat(string(rowHeaders(j)),":",string(mins_(j)),":",string(secs_(j)));
    %         try
    %             t(j) = dn2dt(datenum(char(t_),'yyyy-mm-ddTHH:MM:SS.FFF'));
    %         catch
    %             t(j) = dn2dt(datenum(char(t_),'yyyy-mm-ddTHH:MM:SS'));
    %         end
    %     end
    t = tTarget(fI)-seconds(20); %sort(t);
    nTimes = length(t);

    S = extractWaveforms(t+seconds(0),seconds(secDur),kstnms,kcmpnms,"EC","",true,verboseFlag,2,true,[lfc,hfc,false,transferFlag]);
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
    if saveFlag
        if i < 10
            save(strcat('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/optimal_templates_family_v2_',['0' num2str(i)]),'T');
        else
            save(strcat('~/research/now/plata/plata_repeaters/SANDRO/OTAV_EV/optimal_templates_family_v2_',num2str(i)),'T');
        end
    end
end
toc;
