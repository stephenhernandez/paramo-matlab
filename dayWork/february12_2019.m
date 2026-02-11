%function february12_2019()
% this code is intended to extract long-period or low frequency transients
% from the months of june-july in 2017 at Reventador.
% What I need to do is filter the data in a low frequency band.
% then apply an sta/lta filter to extract pseudo-p times.
% from there, i extract the data using extractSacFromList().

%%
clear; close all; clc;
lfc = 1/4;
hfc = 1/2;
npoles = 3;
tw = 0.02;
winLen = 30;
mph = 3;
plotFlag = false;
envFiltFlag = true;
envHfc = 2/winLen;

%%
snr = [];
tabs = [];
days = datetime(2017,01,01):datetime(2018,12,31);
for i = 1:length(days)
    disp(days(i));
    S = loadSacData(days(i),1,"REVN","HHZ");
    if isnat(S.ref)
        continue;
    else
        S = filterSacData(S,lfc,hfc,npoles,tw);
        %S = intSacData(S);
        %[locs,snr_] = kurtosisPicker(S,winLen,mph,[],plotFlag,envFiltFlag,envHfc);
        [locs,snr_] = stalta(S,winLen,winLen,mph,[],true,plotFlag,envFiltFlag,envHfc);
        t = getTimeVec(S);
        tabs = [tabs; t(locs)]; %#ok<AGROW>
        snr = [snr; snr_]; %#ok<AGROW>
    end
end

%%
close all;
figure();
plot(tabs,(1:length(tabs))','o');

%
%ax = plotSacData(S);

figure();
pI = true(size(tabs));
histogram(tabs(pI),dn2dt(floor(datenum(min(tabs)))):hours(1):dn2dt(ceil(datenum(max(tabs)))));
%unique(dateshift(sort(tabs),'start','day')));
xlim([datetime(2017,06,21) datetime(2017,07,01)]);

%%
[LF,goodI_LF] = extractSacFromList(tabs-seconds(10),tabs+minutes(1),"REVN","HHZ");

%%
% % RETU picker
% S = loadSacData(datetime(2015,04,10),1,"RETU","SHZ");
% S = filterSacData(S,2,6,4,0.02);
% close all; ax = plotSacData([S;intSacData(S)]);
% [locs,snr,staOlta,sosSTA,sosSTA2] = kurtosisPicker(differentiateSacData(S),10,1.5,[],true,true,1/2);
% t = getTimeVec(S);
% tabs = t(locs);
% figure(); plot(tabs,1:length(tabs),'o');
