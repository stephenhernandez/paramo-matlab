%cinchona inversion
clear; close all; clc; 

%cd ~/costa_rica/
%load ~/costa_rica/event_templates/08January2009_templates_acc_20s_50s_4poles
load ~/research/now/costa_rica/data/event_templates/08January2009_templates_acc_20s_50s_4poles
[templates,stnms,stla,stlo,ref,fs] = distillStruct(S);
keepI = true(3,length(stla));
ns = length(stla);
plotFlag = true;
dcFlag = false;
flag3 = true;
minlon = -86.5; 
maxlon = -83.5; 
minlat = 9; 
maxlat = 11;

% USGS Source Parameters
% lat: 10.197
% lon = -84.159
% depth = 4.5 km
% 08 January 2009, 19:21:34

%%
% lats = 10.197;
% lons = -84.159;
% depth = 4.5;
lats = (10.16-0.04):0.04:(10.16+0.04);
%lons = (-84.18-0.16):0.04:(-84.18+0.16);
lons = (-84.06-0.04):0.02:(-84.06+0.04);
depth = 14.5;
[lats,lons,depth] = meshgrid(lats,lons,depth);
lats = lats(:);
lons = lons(:);
depth = depth(:);

units = 'acc';
wFlag = false; %true if you want to apply weighting function
cFlag = false; %true if you want to constrain solution
gmtflag = false; %true if you want to print gmt psmeca file
synthflag = true; % true if you wanna use nicoya synthetics
threeDFlag = true;
if threeDFlag
    %save original data
    depth_orig = depth;
    lons_orig = lons;
    lats_orig = lats;
    
    depth_ = (0.5:1:20.5)';
    depth = repmat(depth_,size(lats));
    lats = kron(lats,ones(size(depth_)));
    lons = kron(lons,ones(size(depth_)));
    clear depth_
    depth = depth(:);
end

%%
lfc = 1/50;
hfc = 1/20;
% lfc = 1/60;
% hfc = 1/15;
npoles = 4;

PARENT_DIR = datestr(now,29);
if ~exist(PARENT_DIR,'dir')
    SUCCESSID = mkdir(PARENT_DIR);
end

%%
t0 = 4;
th = 1;
secDur = 120;
[T0,TH,secDur] = meshgrid(t0,th,secDur);
T0 = T0(:);
TH = TH(:);
secDur = secDur(:);

filtflag = true;
error_tmp = NaN(length(lons),length(T0));
for i = 1:length(T0)
    disp([num2str(i),'/', num2str(length(T0)),' ' num2str(T0(i)), ' ', num2str(TH(i))])
    t0 = T0(i);
    th = TH(i);
    secDurTmp = secDur(i);
    
    disp(length(lats))
    tic;
    parfor k = 1:length(depth)
        disp(k)
        [~,m_tmp(k,:),mw,~,~,error_tmp(k,i),err_lb,err_ub,duration] = ...
            tensor_inversion_4(keepI,lats(k),lons(k),depth(k),templates,stla,stlo,stnms,fs,...
            t0,filtflag,th,lfc,hfc,flag3,synthflag,gmtflag,...
            {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},cFlag,wFlag,npoles,...
            secDurTmp,units);
        dc(k,i) = percentDC(m_tmp(k,:)');
    end
    toc;
    
    if dcFlag
        [maxErr,maxErrI] = max(error_tmp(:,i).*dc(:,i));
    else
        [maxErr,maxErrI] = max(error_tmp(:,i));
    end
    error(i) = maxErr;
    errLat(i) = lats(maxErrI);
    errLon(i) = lons(maxErrI);
    errDepth(i) = depth(maxErrI);
    
    if plotFlag
        if dcFlag
            [~,mi] = max(error_tmp(:,i).*dc(:,i));
        else
            [~,mi] = max(error_tmp(:,i));
        end
        
        [synth,m_,mw(i,1),Gbig_tmp,observables,err(i),elb(i),eub(i),dur,dist,azs,bazs] = ...
            tensor_inversion_4(keepI,lats(mi),lons(mi),depth(mi),templates,stla,...
            stlo,stnms,fs,t0,filtflag,th,lfc,hfc,flag3,true,plotFlag,...
            {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,...
            false,npoles,secDurTmp,units);
        GBig(:,:,i) = Gbig_tmp;
        synthBig(1:length(synth),i) = synth;
        obsBig(1:length(synth),i) = observables;
        mBig(:,i) = m_';
        durations(:,i) = dur;
        dists(:,i) = dist;
        [dc_(i),clvd_(i)] = percentDC(mBig(:,i));
        disp(['Best Mw: ',num2str(mw(i))]);
    end
    
    disp(['Best Lon.: ',num2str(lons(maxErrI))]);
    disp(['Best Lat.: ',num2str(lats(maxErrI))]);
    disp(['Best Depth.: ',num2str(depth(maxErrI))]);
    disp(['Best Variance Reduction: ',num2str(error_tmp(maxErrI,i))]);
    disp(['Best DC: ',num2str(dc(maxErrI,i))]);
    disp(['Best DC-Scaled Variance Reduction: ',num2str(error_tmp(maxErrI,i).*dc(maxErrI,i)/100)]);
    
    [strike,dip,rake] = mt2sdr(mBig(:,i)');
    [strike_a,dip_a,rake_a] = auxplane([strike, dip, rake]);
    disp('nodal plane 1')
    disp(round([strike dip rake]))
    disp('nodal plane 2')
    disp(round([strike_a dip_a rake_a]))
end

%%
close all; 
plot_best_solution;