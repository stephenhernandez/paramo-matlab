function [Msangay,Err_sangay,mlv,z2p,dmlv] = sangayMagnitudeCalculation(tabs,z2p,plotFlag)
%%
if nargin < 3
    plotFlag = false;
end

%%
load('~/masa/old/research/now/sangay/sangay_magnitude_attenuation_parameters.mat','Mparams','kstnms','dM'); %<-- dont use this one even though its better  %116k (588/573, 12 big) %119 (20 big)
%load('~/research/now/sangay/sangay_magnitude_attenuation_parameters_2.mat','Mparams','kstnms','dM'); %<-- use this one %% 114k (569/583, 12 big) % 117k (16 big)

alpha = Mparams(1);
beta = Mparams(2);
gamma = Mparams(3);

%% remove saga from station list and dM list
Mkstnms = kstnms(2:end);
dM = dM(2:end);

%%
[stla,stlo] = metaDataFromStationList(Mkstnms); %no saga
refEllipse = referenceEllipsoid('wgs84');
dist = distance(stla,stlo,-2.00535,-78.341294,refEllipse)*1e-3;

%%
load('~/masa/old/research/now/sangay/raw2wa_sangayRegionalSensors','b');
b_kstnm = ["BPAT";...
    "BMAS";...
    "BULB";...
    "BRUN";...
    "PUYO";...
    "TAMH";...
    "PORT";...
    "PKYU";...
    "TAIS"];

%%
kstnmMaster = ["PUYO";...
    "BULB";...
    "TAIS";...
    "TAMH";...
    "BMAS";...
    "BPAT";...
    "PKYU";...
    "PORT";...
    "BRUN"];

%%
mlv = NaN(size(z2p));
dmlv = mlv;
lK = min([length(kstnmMaster) size(z2p,2)]);
for i = 1:lK
    kstnm_ = kstnmMaster(i);
    
    %% convert to pseudo wood-anderson
    [lia,locb] = ismember(kstnm_,b_kstnm);
    if ~lia
        fprintf(2,'Something went wrong with WA conversion, fix!\n');
        return;
    end
    
    %disp(['locb: ',num2str(locb)]);
    if strcmp("BULB",kstnm_)
        raw = z2p(:,i);
        bI = tabs >= datetime(2016,01,32) & tabs <= datetime(2017,01,91);
        raw(bI) = raw(bI)/4;
        z2p(:,i) = raw;
    end
    
    %%
    raw = z2p(:,i);
    rI = raw > 0.1 & isfinite(raw);
    b_ = b(:,locb);
    raw(rI) = 10.^(b_(1) + b_(2)*log10(raw(rI)));
    z2p(:,i) = raw;
    
    %% convert to custom local magnitude unique to Sangay
    [lia,locb] = ismember(kstnm_,Mkstnms);
    if ~lia
        fprintf(2,'Something went wrong with mag params, fix!\n');
        return;
    end
    a = z2p(:,i); %PSEUDO-WOODANDERSON
    aI = a > 1e-3 & isfinite(a); %<-- this step is important!
    mlv(aI,i) = log10(a(aI)) + alpha*log10(dist(locb)) + beta*dist(locb) + gamma - dM(locb);
    %dmlv(aI,i) = mlv(aI,i);
end

%%
Msangay = nanmedian(mlv,2);
dmlv = mlv;
dmlv = -(Msangay - dmlv);
% Err_sangay = std(mlv-nanmedian(mlv,2),0,2,'omitnan');
% Err_sangay = mad(mlv,1,2);
Err_sangay = std(mlv,0,2,'omitnan');

%%
%plotFlag = false;
if plotFlag
    figure('units','normalized','outerposition',[0 0 0.8 1]);
    ax_(1) = subplot(211);
    semilogy(tabs,z2p,'.'); zoom on; grid on;
    title('individual wa amps');
    ax_(2) = subplot(212);
    plot(tabs,mlv,'.'); zoom on; grid on;
    title('individual mags');
    linkaxes(ax_,'x');
    
    %%
    figure('units','normalized','outerposition',[0 0 0.8 1]);
    ax(1) = subplot(211);
    plot(tabs,Msangay,'.'); zoom on; grid on;
    ax(2) = subplot(212);
    plot(tabs,Err_sangay,'.'); zoom on; grid on;
    linkaxes(ax,'x');
end
