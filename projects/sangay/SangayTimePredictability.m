%function
%SangayTimePredictability
% %clear; close all; clc;

% cumEnergy = [1.0371121e12; ...
%     1.3041990e12;...
%     1.914499e12;...
%     2.231565e12;...
%     2.5058e12];

cumEnergyOrig = cumsum(medfiltSH(energyMag,Nmed));

tStart = [datetime(2015,01,02,11,59,42); ...
    datetime(2016,03,05,07,59,12); ...
    datetime(2017,07,15,10,48,55); ...
    datetime(2018,08,11,06,50,31); ...
    datetime(2019,05,06,06,34,02)];
tStart = tStart + minutes(1);

%%
cumEnergy = [];
figure(7);
hold on;
for i = 1:length(tStart)
    fI = find(tOrig <= tStart(i),1,'last');
    cumEnergy_ = cumEnergyOrig(fI);
    plot(tStart(i),cumEnergy_,'p');
    cumEnergy = [cumEnergy; cumEnergy_];
end

%%
tStart = datenum(tStart);
p = polyfit(tStart,cumEnergy,1);
tq = (datetime(2013,01,01):dn2dt(ceil(now)+365))';
b = robustfit(tStart,cumEnergy);
b = flipud(b);
yq = polyval(b,datenum(tq));

hold on;
pp = plot(tq,yq,'-','linewidth',2);
pp.Color(4) = 0.75;

%% define repose times and fit robust line through them
tStart = [datetime(2015,01,02,11,59,42); ...
    datetime(2016,03,05,07,59,12); ...
    datetime(2017,07,15,10,48,55); ...
    datetime(2018,08,11,06,50,31); ...
    datetime(2019,05,06,06,34,02)];
tStart = tStart + minutes(1);

% reposeTimes = [datetime(2013,06,02,22,15,42) tStart(1); ...
%     datetime(2015,04,19,07,16,00) tStart(2); ...
%     datetime(2016,07,31,09,27,00) tStart(3); ...
%     datetime(2017,10,28,03,39,00) tStart(4); ...
%     datetime(2018,12,04,19,18,00) tStart(5); ...
%     datetime(2021,03,11,08,55,00) tOrig(end)];

reposeTimes = [datetime(2013,06,02,22,15,42) tStart(1); ...
    datetime(2015,04,19,07,16,00) tStart(2); ...
    datetime(2016,07,31,09,27,00) tStart(3); ...
    datetime(2017,10,28,03,39,00) tStart(4); ...
    datetime(2018,12,04,19,18,00) tStart(5); ...
    datetime(2021,03,11,08,55,00) datetime(2022,01,01)];

brepose = [];
predictedTime = [];
for i = 1:size(reposeTimes,1)
    treposeStart_ = reposeTimes(i,1);
    treposeEnd_ = reposeTimes(i,2);

    tI = tOrig >= treposeStart_ & tOrig <= treposeEnd_;
    b_ = robustfit(datenum(tOrig(tI)),cumEnergyOrig(tI));
    b_ = flipud(b_);
    disp(b_);

    tq = (treposeStart_-365:treposeEnd_+365)';
    yq = polyval(b_,datenum(tq));

    hold on;
    pp = plot(tq,yq,'-','linewidth',3);
    pp.Color(4) = 0.5;
    brepose = [brepose; b_'];

    predictedTime_ = (b(2) - b_(2))./(b_(1) - b(1));
    predictedTime_ = dn2dt(predictedTime_);
    predictedTime = [predictedTime; predictedTime_];
    fprintf('predicted_time: %s\n',datestr(predictedTime_));
end
zoom on;
grid on;

%
figure(); 
ntest = 35;
plot(tOrig(2:end),...
    3600*(10.^medfiltSH(m(2:end),ntest))./medfiltSH(seconds(diff(tOrig)),ntest),'.'); 
zoom on; grid on;

figure(8); 
ax = gca; ax.YScale = 'log';
zoom on; grid on;

%%
nDays = 7;
[rate,~,medianMagsFixedTimeWin] = ...
    t2r(tOrig,days(nDays),eqmag);
figure(); 
axx(1) = subplot(311);
semilogy(tOrig,rate/nDays,'.'); grid on; hold on;

axx(2) = subplot(312);
semilogy(tOrig,medianMagsFixedTimeWin,'.'); grid on; hold on;

axx(3) = subplot(313); 
semilogy(tOrig,(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); 
zoom on; grid on; hold on;

linkaxes(axx,'x');

nDays = 1;
[rate,~,medianMagsFixedTimeWin] = ...
    t2r(tOrig,days(nDays),eqmag);

semilogy(axx(1),tOrig,rate/nDays,'.'); grid on;

semilogy(axx(2),tOrig,medianMagsFixedTimeWin,'.'); grid on; 

semilogy(axx(3),tOrig,(10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays,'.'); 
zoom on; grid on;

%%
% averageDailyEnergy7Days = (10.^(1.5*medianMagsFixedTimeWin + 4.8)).*rate/nDays;
% tQuery = sort((datetime(2022,09,15,15,00,00):-7:datetime(2019,05,06))');
% Yq = interp1(datenum(tOrig),averageDailyEnergy7Days,datenum(tQuery));
% -1 + Yq(3:end)./Yq(1:end-2)
% -1 + Yq(2:end)./Yq(1:end-1)
