%july29b_2020
clear; close all; clc;

%%
cd ~/research/now/reventador/

%%
refEllipse = referenceEllipsoid('wgs84');
center_lon = -77.6581;
center_lat = -0.0810;

%%
%listOfSensors = ["REVN","REVS"];
listOfSensors = ["CASC","REVN","REVS"];
[stla,stlo,stel] = metaDataFromStationList(listOfSensors);
d_ = distance(stla,stlo,center_lat,center_lon,refEllipse)*1e-3;
lS = length(listOfSensors);

%%
load tMaster.mat

load CASC_HHZ_2.mat
timeRef = min([tMaster(1) t(1)]);
t1 = seconds(t - timeRef);
t2 = seconds(tMaster - timeRef);
lia = ismembertol(t1,t2,2,'DataScale',1);
p2p_CASC = p2p(lia);

%%
load REVN_HHZ_2.mat
timeRef = min([tMaster(1) t(1)]);
t1 = seconds(t - timeRef);
t2 = seconds(tMaster - timeRef);
lia = ismembertol(t1,t2,2,'DataScale',1);
p2p_REVN = p2p(lia);

%%
load REVS_HHZ_2.mat
timeRef = min([tMaster(1) t(1)]);
t1 = seconds(t - timeRef);
t2 = seconds(tMaster - timeRef);
lia = ismembertol(t1,t2,2,'DataScale',1);
p2p_REVS = p2p(lia);

%%
CASC_ResponseInfo(1).tstart = datetime(2012,07,12);
CASC_ResponseInfo(1).sensitivity = 3.141950e+08;

%%
REVN_ResponseInfo(1).tstart = datetime(2012,01,257);
REVN_ResponseInfo(1).sensitivity = 3.141950e+08;
REVN_ResponseInfo(2).tstart = datetime(2015,01,079);
REVN_ResponseInfo(2).sensitivity = 3.141950e+08;
REVN_ResponseInfo(3).tstart = datetime(2019,01,256);
REVN_ResponseInfo(3).sensitivity = 1.256780e+09;

%%
REVS_ResponseInfo(1).tstart = datetime(2012,01,257);
REVS_ResponseInfo(1).sensitivity = 3.141950e+08;
REVS_ResponseInfo(2).tstart = datetime(2015,01,105);
REVS_ResponseInfo(2).sensitivity = 3.141950e+08;

%%
for i = 1:lS
    disp(i);
    kstnm_ = listOfSensors(i);
    
    %% loop through periods with distinct sensor sensitivities
    if strcmp("CASC",kstnm_) %CASC
        %%
        respStart = pull(CASC_ResponseInfo,'tstart');
        respSensitivity = pull(CASC_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = tMaster >= respStart(ii) & tMaster <= respStart(ii+1);
                else
                    tII = tMaster >= respStart(ii);
                end
                if sum(tII)
                    p2p_CASC(tII) = p2p_(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_CASC = p2p_CASC/respSensitivity;
        end
        
    elseif strcmp("REVN",kstnm_)  %REVN
        %%
        respStart = pull(REVN_ResponseInfo,'tstart');
        respSensitivity = pull(REVN_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = tMaster >= respStart(ii) & tMaster <= respStart(ii+1);
                else
                    tII = tMaster >= respStart(ii);
                end
                if sum(tII)
                    p2p_REVN(tII) = p2p_REVN(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_REVN = p2p_REVN/respSensitivity;
        end
    else
        %%
        respStart = pull(REVS_ResponseInfo,'tstart');
        respSensitivity = pull(REVS_ResponseInfo,'sensitivity');
        lGains = length(respStart);
        if lGains > 1
            for ii = 1:lGains
                if ii ~= lGains
                    tII = tMaster >= respStart(ii) & tMaster <= respStart(ii+1);
                else
                    tII = tMaster >= respStart(ii);
                end
                if sum(tII)
                    p2p_REVS(tII) = p2p_REVS(tII) / respSensitivity(ii);
                end
            end
        else
            p2p_REVS = p2p_REVS/respSensitivity;
        end
    end
end

%%
close all;
figure();
for i = 1:lS
    kstnm_ = listOfSensors(i);
    if strcmp("CASC",kstnm_) 
        p2p_ = p2p_CASC;
    elseif strcmp("REVN",kstnm_) 
        p2p_ = p2p_REVN;
    else
        p2p_ = p2p_REVS;
    end
    ax_(i) = subplot(lS,1,i);
    semilogy(tMaster,p2p_,'o'); grid on;
    title(listOfSensors(i));
end
zoom on;
linkaxes(ax_,'x');

%%
beta = (0.5:0.1:6.5);
Q = (5:5:250)';
f = 0.4;
errVec = NaN(size(Q));

%%
ampOrig = [];
for i = 1:lS
    kstnm_ = listOfSensors(i);
    if strcmp("CASC",kstnm_) 
        ampOrig = [ampOrig log10(p2p_CASC)];
    elseif strcmp("REVN",kstnm_) 
        ampOrig = [ampOrig log10(p2p_REVN)];
    else
        ampOrig = [ampOrig log10(p2p_REVS)];
    end
end
magEst = ampOrig;

%%
[X,Y] = meshgrid(beta,Q);
sizeX = size(X);
X = X(:);
Y = Y(:);
lX = length(X);
errVec = NaN(size(X));

for i = 1:lX
    disp(i);
    beta_ = X(i);
    Q_ = Y(i);
    B = pi*f/(Q_*beta_);
    
    for j = 1:lS
        %magEst(:,j) = magEst(:,j) + 0.5*log10(d_(j));  %surface wave model
        magEst(:,j) = magEst(:,j) + log10(d_(j));       %body wave model
        magEst(:,j) = magEst(:,j) + B*d_(j)*log10(exp(1)) + 7;
    end
    magErr = mad(magEst,1,2);
    errVec(i) = mad(magErr);
    magEst = ampOrig; %clobber the old magEst
end

%


