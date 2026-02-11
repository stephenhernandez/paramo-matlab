% cotopaxi_trends
clear; clc;
%[tsheet,type,Tsp,amp,mag,codaDuration,period,refStation] = loadVolcanicEvents('guagua',6); 
%tStart = datetime(2000,01,01);

[tsheet,type,Tsp,amp,mag,codaDuration,period,refStation] = loadVolcanicEvents('cotopaxi',6);
tStart = datetime(2014,07,01);
ampThresh = 2e3;

%%
tI = tsheet >= tStart;
tsheet(~tI) = [];
type(~tI) = [];
Tsp(~tI) = [];
amp(~tI) = [];
mag(~tI) = [];
codaDuration(~tI) = [];
period(~tI) = [];
refStation(~tI) = [];

%%
[tsheet,amp,removeIndices] = removeRepeatedMatches(tsheet,amp,20);
for kk = 1:length(removeIndices)
    rI = removeIndices{kk};
    type(rI) = [];
    Tsp(rI) = [];
    codaDuration(rI) = [];
    period(rI) = [];
    refStation(rI) = [];
end

%%
uniqueTypes = unique(type);
ltypes = length(uniqueTypes);
Ndays = length(dateshift(min(tsheet),'start','day'):dateshift(max(tsheet),'start','day'))-1;
N = NaN(Ndays,ltypes);
legString = [];
for k = 1:ltypes
    uniqueTypes_ = uniqueTypes(k);
    legString = [legString; string(uniqueTypes_)];
    tI2 = strcmp(uniqueTypes_,type);% & tsheet >= tStart;
    sumType = sum(tI2);
    [N_,edges] = histcounts(tsheet(tI2),dateshift(min(tsheet),'start','day'):dateshift(max(tsheet),'start','day'));
    N(:,k) = N_';
    disp(strcat(uniqueTypes_,': ',num2str(sumType),',mean:',num2str(nanmean(N_')),',ci:',num2str(prctile(N_',[5 95]))));
end
disp('-----');
disp(' ');
edges = edges(1:end-1)';

%%
close all; 
figure(); plot(edges,N(:,[2 7]),'.'); zoom on;
%figure(); plot(edges,N(:,[4 9]),'.'); zoom on; %guagua
legend('LP','VT');
grid on

%%
figure(); 
plot(tsheet(tI2),amp(tI2),'o'); zoom on;
ax = gca;
ax.YScale = 'log';


t1 = tI2 & tsheet >= tStart & amp >= ampThresh;
t2 = tI2 & tsheet >= tStart & amp < ampThresh;
tBigLP = tsheet(t1);
tSmallLP = tsheet(t2);

[NbigLP,edgesBigLP] = histcounts(tBigLP,dateshift(min(tBigLP),'start','day'):dateshift(max(tBigLP),'end','day'));
[NsmallLP,edgesSmallLP] = histcounts(tSmallLP,dateshift(min(tSmallLP),'start','day'):dateshift(max(tSmallLP),'end','day'));
NbigLP = NbigLP';
edgesBigLP = edgesBigLP(2:end)';
NsmallLP = NsmallLP';
edgesSmallLP = edgesSmallLP(2:end)';


figure(); 
plot(edgesBigLP,cumsum(NbigLP),'.'); 
ylabel('Eventos mayores de 2000 cuentas'); 
yyaxis right; 
plot(edgesSmallLP,cumsum(NsmallLP),'.'); 
zoom on;
yyaxis right
title('VT');
ylabel('Eventos de menos de 2000 cuentas');