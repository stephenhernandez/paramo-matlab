function importTemplateMatchData(dataFileName)

%dataFileName = 'dayTemplateSearch_2022.10.06.txt';
formatSpec = '%d %d %d %d %d %f %f %f %f %s %d %f %f';
fig = fopen(dataFileName);
C = textscan(fid,formatSpec);
fclose(fid);

evidMain = C{10};
tMain = datetime(C{1},C{2},C{3},C{4},C{5},C{6});

dMag = C{7};
ampMain = C{8};
ccMain = C{9};
templateNumber = C{10};
madMain = C{11};
magMain = C{12};

close all
[t,cc,removeIndices] = removeRepeatedMatches(tMain,ccMain,20);
figure(1); hold on; plot(t,1:length(t),'.'); zoom on;
amp_ = ampMain; evid_ = evidMain;
for i = 1:length(removeIndices)
    rI = removeIndices{i};
    amp_(rI) = [];
    evid_(rI) = [];
end
figure(); semilogy(t,amp_,'.'); zoom on;
figure(); plot(t,cc,'.'); zoom on;
