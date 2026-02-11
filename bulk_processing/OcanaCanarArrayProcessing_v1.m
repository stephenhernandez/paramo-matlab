clear; close all
cd ~/data/nodes/JG_Canar/
Zfiles = dir('*/*Z.miniseed');
names = pull(Zfiles,'name',"");
split_names = split(names,".");
tDayTmp = datetime(str2double(split_names(:,5)),str2double(split_names(:,6)),str2double(split_names(:,7)));

uniqDays = unique(tDayTmp);
[yyyy,mm,dd] = datevec(uniqDays);

knetwk = "SS";
khole = "SW";

nUsed = [];
tMain = [];
ampMain = [];
bazMain = [];
velMain = [];
meanCCsMain = [];
medCCsMain = [];
for i = length(uniqDays)-1
    searchStr = sprintf('%04d.%02d.%02d',yyyy(i),mm(i),dd(i));
    disp(searchStr);
    I = contains(names,searchStr);
    dayFileNames = names(I);
    S = populateWaveforms(sum(I));

    for j = 1:sum(I)
        kstnm = char(dayFileNames(j));
        kstnm = string(kstnm(6:9));
        S(j) = readMiniSeed(fullfile(kstnm,dayFileNames(j)));
    end
    knetwk = pull(S,'knetwk');
    kstnm2 = pull(S,'kstnm');
    khole = pull(S,'khole');
    kcmpnm = pull(S,'kcmpnm');
    allSNCLs = strcat(knetwk,kstnm2,khole,kcmpnm);
    uniqSNCLs = unique(allSNCLs);
    lUniqSNCLs = length(uniqSNCLs);
    n = 0;
    S2 = populateWaveforms(lUniqSNCLs);
    for j = 1:lUniqSNCLs
        n = n + 1;
        uniqSNCLs_ = uniqSNCLs(j);
        lia = ismember(allSNCLs,uniqSNCLs_);
        if sum(lia) == 1
            S2(n) = S(lia);
            continue;
        end

        fprintf('merge: %s\n',uniqSNCLs_);
        S_ = S(lia);
        S_ = mergeWaveforms(S_);
        S2(n) = S_;
    end

    %%
    % [~,nUsed_,tMain_,ampMain_,bazMain_,velMain_,meanCCsMain_,medCCsMain_,Sf] = ...
    %     process_OcanaCanar_v1((S2));
    [~,nUsed_,tMain_,ampMain_,bazMain_,velMain_,meanCCsMain_,medCCsMain_,Sf] = ...
        process_OcanaCanar_v1(differentiateWaveforms(S2));
    nUsed = [nUsed; nUsed_];
    tMain = [tMain; tMain_];
    ampMain = [ampMain; ampMain_];
    bazMain = [bazMain; bazMain_];
    velMain = [velMain; velMain_];
    meanCCsMain = [meanCCsMain; meanCCsMain_];
    medCCsMain = [medCCsMain; medCCsMain_];
end

%%
tMainOrig = tMain;
medCCsMainOrig = medCCsMain;
nUsedOrig = nUsed;
velMainOrig = velMain;
bazMainOrig = bazMain;
ampMainOrig = ampMain;
meanCCsMainOrig = meanCCsMain;

% minSeparation = 10;
% [tMain,ampMain,nUsed,velMain,bazMain,medCCsMain,meanCCsMain] = ...
%     filterCatalog(tMainOrig,ampMainOrig,minSeparation,nUsedOrig,velMainOrig,bazMainOrig,medCCsMainOrig,meanCCsMainOrig);

minSeparation = 40;
[tMain,medCCsMain,nUsed,velMain,bazMain,ampMain,meanCCsMain] = ...
    filterCatalog(tMain,medCCsMain,minSeparation,nUsed,velMain,bazMain,ampMain,meanCCsMain);

%%
close all;
figure();
ax(1) = subplot(311);
semilogy(tMain,velMain*1e-3,'.'); grid on; title('velocity');

ax(2,1) = subplot(312);
plot(tMain,bazMain,'.'); grid on; title('back-azimuth');

ax(3,1) = subplot(313);
semilogy(tMain,ampMain,'.'); grid on; title('amplitude');
zoom on;
linkaxes(ax,'x');







