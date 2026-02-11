function [z2p,tcut,kstnms] = rsam2magnitude(nMinutes,tStart,tEnd)
if nargin < 1
    nMinutes = 10;
end

if nargin < 2
    tStart = datetime(2021,03,01);
end

if nargin < 3
    tEnd = datetime(2021,04,01);
end

%%
load('~/research/now/sangay/sangay_svd_basis_functions_withSNR_v8','kstnm');

%% preallocate
ntot = 1 + (seconds(tEnd - tStart)/30);
lk = length(kstnm);
z2p = NaN(ntot,lk);
kstnms = repmat("",lk,1);

%%
for i = 1:lk
    kstnm_ = kstnm(i);
    chanStr = "HHZ";
    if strcmp(kstnm_,"BULB") || strcmp(kstnm_,"BMAS") ||  strcmp(kstnm_,"BPAT") ||  strcmp(kstnm_,"BRUN")
        chanStr = "BHZ";
    end
    fname = strcat("EC.",kstnm_,"..",chanStr,"_0.6Hz1.2Hz_30DUR_MedAmpRMS_preFiltFalse.mat");
    fname = fullfile('~','products','rsam',char(fname));
    try
        fprintf('trying to read: %s\n',fname);
        load(fname,'S');
    catch
        disp('something went wrong,trying again');
        pause(5);
        try
            fprintf('trying to read: %s\n',fname);
            load(fname,'S');
        catch
            disp('ok, that didnt work either, skipping entirely...')
            continue;
        end
    end
    
    %%
    S = medfiltWaveforms(S,1+nMinutes*2);
    %S = convWaveforms(S,nMinutes*2);
    kstnms(i) = S.kstnm;
    t = getTimeVec(S);
    d = S.d;
    d(d == 0) = NaN;
    
    tday = dateshift(t,'start','day');
    tm = minutes(t - tday);
    tDayStart = tday(1);
    tDayEnd = tday(end)+1;
    
    %maxWeeks = floor(days(tDayEnd - tDayStart)/7);
    nWeeks = floor(length(tm)/2880/7);
    maxWeeks = nWeeks;
    
    %disp([nWeeks maxWeeks])
    %tm2 = reshape(tm(1:2880*7*nWeeks),[2880*7 nWeeks]);
    d2 = reshape(d(1:2880*7*nWeeks),[2880*7 nWeeks]);
    tthisweek = t(2880*7*(nWeeks-maxWeeks)+1:end);
    dmed = nanmedian(d2,2);
    bgBounds = repmat(dmed,1,2);
    for j = 1:size(d2,1)
        Y = prctile(d2(j,:),[10 90]);
        bgBounds(j,:) = Y;
    end
    
    %%
    dmed = [repmat(dmed,nWeeks,1); dmed(1:length(d)-2880*7*nWeeks)];
    %     figure(); pp = plot(tthisweek,d(2880*7*(nWeeks-maxWeeks)+1:end),'.-','linewidth',1); zoom on; grid on; pp.Color(4) = 0.5;
    %     hold on;
    %     ax = gca;
    %     ax.YScale = 'log';
    %     plot(tthisweek,dmed(1:length(tthisweek)),'-','linewidth',2); zoom on; grid on;
    %     bgBounds2 = [repmat(bgBounds(:,1),nWeeks,1); bgBounds(1:length(d)-2880*7*nWeeks,1)];
    %     bgBounds2(:,2) = [repmat(bgBounds(:,2),nWeeks,1); bgBounds(1:length(d)-2880*7*nWeeks,2)];
    %     plot(tthisweek,bgBounds2(1:length(tthisweek),1),'k-','linewidth',2); zoom on; grid on;
    %     plot(tthisweek,bgBounds2(1:length(tthisweek),2),'k-','linewidth',2); zoom on; grid on;
    %     xlim([datetime(2021,02,01) datetime(2021,04,14)]);
    %     %xlim([datetime(2016,01,01) datetime(2016,02,01)]);
    tI = tthisweek >= tStart & tthisweek <= tEnd;
    sumti = sum(tI);
    if ~sumti
        disp('no data');
        continue;
    end
    dcut = d(2880*7*(nWeeks-maxWeeks)+1:end);
    dcut = dcut(tI);
    
    %%
    lcut = sumti;
    if lcut > ntot
        lcut = ntot;
    end
    z2p(1:lcut,i) = dcut;
    
    %%
    if ~exist('tcut','var')
        tcut = tthisweek(tI);
        cutcurr = sumti;
    elseif sumti > cutcurr
        tcut = tthisweek(tI);
    end
end
