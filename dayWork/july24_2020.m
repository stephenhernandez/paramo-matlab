clear; close all; clc;
cd ~/research/now/cotopaxi/;
load bref_data_28jun.mat;

%%
mI = maxAmpRMS >= 4e2;
ampFilt = maxAmpRMS(mI);
t = tref(mI);

%%
[N,edges] = histcounts(t,dateshift(min(t),'start','day'):hours(12):dateshift(max(t),'end','day'));
N = N';
edges = edges(1:end-1)';
[~,edges2] = histcounts(t,dateshift(min(t),'start','day'):hours(12):dateshift(max(t),'end','day'));
figure(); plot(edges,N,'o'); zoom on;
%%
[ccRate,lagsRate] = xcorr(detrend(N));
%%
figure(); plot(lagsRate/2,zpkFilter(ccRate,-inf,1/5,1,1,true)); zoom on; title('Rate AutoCorrelation');
%%
for i = 1:length(edges2)-1
    aI = t>=edges2(i) & t <= edges2(i+1); nBin = sum(aI);
    if nBin > 2
        thisVec = ampFilt(aI); energyMedSum(i) = sum(thisVec(3:end));
    elseif nBin > 0
        thisVec = ampFilt(aI); energyMedSum(i) = sum(thisVec);
    else
        energyMedSum(i) = 0;
    end
end

%%
figure(); semilogy(edges,energyMedSum','.'); zoom on;

%%
[ccEnergy,lagsEnergy] = xcorr(detrend(energyMedSum'));
figure(); plot(lagsEnergy/2,zpkFilter(ccEnergy,-inf,1/5,1,1,true)); zoom on; title('Energy AutoCorrelation');

%%
Nwin = 100;
mamp = cutWindows(ampFilt,Nwin,Nwin-1,false);
[winMax,maxI] = max(mamp);
for i = 1:size(mamp,2)
    m_ = sort(mamp(:,i));
    %m_(maxI(i)) = [];
    mamp(1:Nwin-10,i) = m_(1:Nwin-10);
end
mamp = mamp(1:Nwin-10,:);

%%
figure(); semilogy(mean(mamp),'.'); zoom on;
sqrtvar = std(mamp)./sqrt(Nwin-10);
ub = mean(mamp) + 1.96*sqrtvar;
lb = mean(mamp) - 1.96*sqrtvar;

%%
figure(); semilogy(t(Nwin:end),mean(mamp),'.'); hold on; semilogy(t(Nwin:end),lb,'-','color',[0.5 0.5 0.5]); semilogy(t(Nwin:end),ub,'-','color',[0.5 0.5 0.5]); zoom on;
[cc,lags] = xcorr(detrend(mean(mamp)));
figure(); plot(lags,cc); zoom on;