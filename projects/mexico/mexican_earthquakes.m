clear; close all; clc;
% load('~/igdata/globalCatalog','t','eqlat','eqlon','eqdepth','eqmag','magType','status','code');
% mI = eqmag >= 7 & t >= datetime(1900,01,01) & ...
%     eqlat >= 13 & eqlat <= 20 & eqlon <= -90 & eqlon >= -107;
% N = sum(mI);

% de 1900 en adelante, han habido 60 eventos >= M7
N = 60;

%%
tStart = datetime(1900,01,01);
tEnd = datetime(2023,01,01);

Nsim = 1e3;
imax = days(tEnd - tStart);
R = randi(imax,[N, Nsim]);

%%
tSynth = days(R) + tStart;
maxRepeats = NaN(Nsim,1);
N3 = maxRepeats;

[~,mm,dd] = datevec(tSynth);
synthCatIgnoreYear = erase(strcat(reshape(string(num2str(mm(:))),N,Nsim),repmat("_",size(R)),reshape(string(num2str(dd(:))),N,Nsim))," ");
for i = 1:Nsim
    gc = groupcounts(synthCatIgnoreYear(:,i));
    gc3 = sum(gc == 3);
    maxgc = max(gc);
    maxRepeats(i) = maxgc;
    N3(i) = gc3;
end

%%
close all;
figure(); plot(maxRepeats,'o'); zoom on;
disp(sum(maxRepeats == 3)/length(maxRepeats))
xlabel('simulation number'); ylabel('max number of repeats for a single day');