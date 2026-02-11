clear;
cd ~;
%eventID = 'igepn2016hobd.txt'; %17 apr 2016 07h13 M6
%eventID = 'igepn2021xhpa.txt'; %28 nov 2021, M7.5 Peru, 10h50ish
eventID = 'igepn2023fkei.txt'; % 18 mar 2023 nazca golfo de guayaquil
%eventID = 'igepn2022oloz.txt'; %chiles/potrerillos 25 jul 2022
%eventID = 'igepn2021wyvu.txt'; %23 nov. puembo M4.5
%eventID = 'igepn2016hnmu.txt'; %16 apr 2016 pedernales

yyyy = eventID(6:9);

%
E = readSCBulletin(string(eventID));
plotFlag = true;
dcFlag = false; %double-couple weighting (not "direct current"!!)
flag3 = true;
threeDFlag = true;

%
origmag = E.mag;
origdepth = E.depth;
origlat = E.lat;
origlon = E.lon;
origt = E.t;

componentsList = "Z";
if flag3
    componentsList = ["Z";"N";"E"]; %read 3 components no matter what
end
ncomps = length(componentsList);
dayStart = dateshift(origt,'end','day');
dfrac = minutes(dayStart - origt);

%
if true%origmag >= 3.
    noise = 600;        % noise window for measuring snr
    signal = 600;       % signal window for measuring snr
    secDur = 5*60;      % window duration for waveform inversion
    lfc = 1/100;
    hfc = 1/50;
    sta = 1/lfc;
    lta = 4*sta;
    swing = 0.08;
    inc = 0.01;
    maxDist = 450;
    minDist = 0;
    t0 = 0;
    th = 0.:0.5:3;

    ampThresh = 1.;     % not more than...
    npoles = 4;
    finalFs = 100;
    units = 'acc';
    tw = 0.01;
    snrThresh = 4;
end

%
nPphases = E.nPphases;
Pphases = E.Pphases;
kstnms = pull(Pphases,'stnm');
kcmpnms = pull(Pphases,'chan');
knetwks = pull(Pphases,'ntwk');
kholes = pull(Pphases,'locid');

%
allSNCLs = [knetwks,kstnms,kholes,kcmpnms];
[allSNCLs,ia] = unique(allSNCLs,'rows');
Pphases = Pphases(ia);

%
kcmpnmstmp = char(allSNCLs(:,4));
kcmpnmstmp = string(kcmpnmstmp(:,1:2));

%
[stla,stlo] = metaDataFromStationList(allSNCLs(:,2),allSNCLs(:,1),allSNCLs(:,4),allSNCLs(:,3));
refEllipse = referenceEllipsoid('wgs84');
d_ = distance(origlat,origlon,stla,stlo,refEllipse)*1e-3;
%dI = d_ <= maxDist & d_ >= minDist & (strcmp(kcmpnmstmp,"HH") | strcmp(kcmpnmstmp,"BH") | strcmp(kcmpnmstmp,"HN"));
dI = strcmp(kcmpnmstmp,"HN");
dI = find(dI);

%
allSNCLs = allSNCLs(dI,:);
d_ = d_(dI);
Pphases = Pphases(dI);

%
[~,dI] = sort(d_);
allSNCLs = allSNCLs(dI,:);
Pphases = Pphases(dI);

%
kstnms = allSNCLs(:,2);
kcmpnms = allSNCLs(:,4);
knetwks = allSNCLs(:,1);
kholes = allSNCLs(:,3);
kcmpnmstmp = char(kcmpnms);
kcmpnmstmp = string(kcmpnmstmp(:,1:2));

%
[stla,stlo] = metaDataFromStationList(allSNCLs(:,2),allSNCLs(:,1),allSNCLs(:,4),allSNCLs(:,3));
[d_,azs] = distance(origlat,origlon,stla,stlo,refEllipse);
[~,bazs] = distance(stla,stlo,origlat,origlon,refEllipse);
d_ = d_*1e-3;

%
lSNCLs = length(kstnms);

S = populateWaveforms(ncomps*lSNCLs);
S2 = S;
n = 0;

tic;
for i = 1:lSNCLs
    for j = 1:ncomps

        S1 = loadWaveforms(dayStart-1,1,kstnms(i),strcat(kcmpnmstmp(i),componentsList(j)),knetwks(i));

        if isnat(S1.ref)
            fprintf('skipping %s\n', kstnms(i));
            continue;
        end

        n = n + 1;

        S1 = differentiateWaveforms(S1);
        %S1 = cutWaveforms(S1,origt-minutes(60),0,minutes(120));
        S1 = cutWaveforms(S1,Pphases(i).t-minutes(1),0,minutes(10));
        S1 = demeanWaveforms(S1);
        S1 = detrendWaveforms(S1);
        S1 = interpolateWaveforms(S1);

        %
        S1.evla = origlat;
        S1.evlo = origlon;
        S1.evdp = origdepth;
        S1.dist = d_(i);
        S1.az = azs(i);
        S1.baz = bazs(i);
        S1.gcarc = km2deg(d_(i));
        S1.user0 = Pphases(i).t;        % analyst's p-pick

        %
        S1 = taperWaveforms(S1,0.01);
        if sum(~isfinite(S1.d))
            n = n-1;
            continue;
        end

        %
        S1orig = S1;
        S1 = transferWaveforms(S1orig,lfc/4,-inf,npoles/2,100,units,1);
        if ~isfinite(S1.stla) || ~isfinite(S1.stlo)
            n = n-1;
            continue;
        end

        %
        tt = taupTime('iasp91',origdepth,'p,P','km',d_(i));
        if isempty(tt)
            n = n-1;
            continue;
        end

        %
        S1.user1 = tt(1).time;          % travel time using tables

        %
        S(n,1) = detrendWaveforms(intWaveforms(detrendWaveforms(S1)));
        %S1 = transfer(S1,wa_zeros,wa_poles,1e3*wa_sensitivity*wa_normalizationFactor,1/10,-inf,npoles,'disp',0);
        %S1 = transferWaveforms(S1orig,1/10,-inf,npoles,100,'disp',1,true);
        S2(n,1) = resampleWaveforms(detrendWaveforms(intWaveforms(detrendWaveforms(S1))),25);

        %
        fprintf('loaded: %s.%s.%s.%s, trace: %d\n',knetwks(i),kstnms(i),kholes(i),strcat(kcmpnmstmp(i),componentsList(j)),n);
        toc;
    end
end

%
% synch data
if ~n
    return;
end

S = S(1:n);
S2 = S2(1:n);
d_ = pull(S,'dist');
stla = pull(S,'stla');
stlo = pull(S,'stlo');

%
% if strcmp(units,'vel')
%     S = differentiateWaveforms(S);
% elseif strcmp(units,'acc')
%     S = differentiateWaveforms(S,2);
% end

% resample and filter
S = resampleWaveforms(S,finalFs);
S = cutWaveforms(S,origt-minutes(10),0,minutes(20));
S2 = cutWaveforms(S2,origt-minutes(10),0,minutes(20));

%
S = detrendWaveforms(S);
S2 = detrendWaveforms(S2);

%
Sorig = S;
S2orig = S2;
dOrig = d_;
stlaOrig = stla;
stloOrig = stlo;

%
close all; 
S = Sorig;
S2 = S2orig;
d_ = dOrig;
stla = stlaOrig;
stlo = stloOrig;

%
S = filterWaveforms(S,lfc,hfc,npoles);

%% get snr and amplitudes
snr = NaN(n,1);
p2p = snr;

for i = 1:n
    d = S2(i).d;
    p2p(i) = max(abs(d));
    
    %
    [snr_,fxx] = freqDomSNR(interpolateWaveforms(S2(i)),S2(i).user0,8*60);
    fI = fxx >= lfc & fxx <= hfc;
    snr(i) = sum(snr_(fI))/sum(fI);
end
clear snr_ fI

%
rI = p2p <= 500;
if sum(rI) < 2
    fprintf('all clipped\n');
    return;
end

b_ = robustfit(log10(d_(rI)),log10(p2p(rI)));
b_ = flipud(b_);
yq = polyval(b_,log10(d_));
normdiff = (p2p-(10.^yq))./(10.^yq);
zscoreProxy = (normdiff-median(normdiff))./mad(normdiff,1);

kcmpnmstmp = char(pull(S2,'kcmpnm'));
kcmpnmstmp = string(kcmpnmstmp(:,1:2));
badI = snr < snrThresh | d_ > maxDist | d_ <= minDist | isnat(pull(S2,'ref')) | abs(zscoreProxy) > 1.75 | (~strcmp(kcmpnmstmp,"HN") & p2p > 500);
goodI = ~badI;

%
figure(); 
loglog(d_,p2p,'.'); zoom on; grid on;
hold on;
figure(1); ll = loglog(d_,10.^yq,'linewidth',3); 
ll.Color(4) = 0.5;
loglog(d_(~goodI),p2p(~goodI),'o'); zoom on; grid on;

%
sumgood = sum(goodI);
if sumgood < 4
    disp('no good events');
    return;
end
S = S(goodI);
S2 = S2(goodI);
snr = snr(goodI);
p2p = p2p(goodI);
d_ = d_(goodI);
kcmpnmstmp = kcmpnmstmp(goodI);

b_ = robustfit(log10(d_),log10(p2p));
b_ = flipud(b_);
yq = polyval(b_,log10(d_));
normdiff = (p2p-(10.^yq))./(10.^yq);
zscoreProxy = (normdiff-median(normdiff))./mad(normdiff,1);
badI = snr < snrThresh | d_ > maxDist | d_ <= minDist | isnat(pull(S2,'ref')) | abs(zscoreProxy) > 2.25 | (~strcmp(kcmpnmstmp,"HN") & p2p > 500);
goodI = ~badI;

figure(1); ll = loglog(d_,10.^yq,'linewidth',3); 
ll.Color(4) = 0.5;
loglog(d_(~goodI),p2p(~goodI),'o'); zoom on; grid on;

%
sumgood = sum(goodI);
if sumgood < 4
    disp('no good events');
    return;
end
S = S(goodI);
snr = snr(goodI);
p2p = p2p(goodI);

%
kstnms = pull(S,'kstnm');
kcmpnms = pull(S,'kcmpnm');
knetwks = pull(S,'knetwk');
kholes = pull(S,'khole');
kcmpnmstmp = char(kcmpnms);
kcmpnmstmp = string(kcmpnmstmp(:,1:2));

%
[stla,stlo,stel] = metaDataFromStationList(kstnms,knetwks,kcmpnms,kholes);
d_ = distance(stla,stlo,origlat,origlon,refEllipse)*1e-3;

%
figure(); 
subplot(121);
plot(stlo,stla,'.'); 
axis equal; zoom on; grid on; hold on; 
text(stlo,stla,kstnms);
hold on; plot(origlon,origlat,'p','linewidth',2,'markersize',20);

load ~/igdata/ec_boundaries.mat
figure(2); hold on;
plot(lonEC,latEC,'-');
axis([min([origlon; stlo])-1 max([origlon; stlo])+1 ...
    min([origlat; stla])-1 max([origlat; stla])+1]);

subplot(122);
plot(pull(S,'dist'),pull(S,'az'),'.'); zoom on; grid on; hold on;

% cut data
S = syncWaveforms(S);
S = cutWaveforms(S,origt,0,secDur);
S = detrendWaveforms(S);

%
Sdata = pull(S);
Sdata = normalizeWaveforms(detrend(Sdata));
af = 50;
figure(); hold on; 
for i = 1:length(d_)
    ll = plot(af.*Sdata(:,i)+d_(i),'linewidth',2); 
    ll.Color(4) = 0.75; 
end
zoom on; grid on;
clear Sdata;

%
lats = (origlat-swing):inc:(origlat+swing);
lons = (origlon-swing):inc:(origlon+swing);

depth = origdepth;
depth = floor(depth)+0.5;

minDepth = depth - 20;
minDepth = max([minDepth 0.5]);
maxDepth = depth + 20;
maxDepth = min([maxDepth 49.5]);

minlon = floor(min([-81.5 min([min(stloOrig) (origlon-swing)])]));   %origlon-3.;
maxlon = ceil(max([max(stloOrig) (origlon+swing)]));    %origlon+3.;
minlat = floor(min([min(stlaOrig) (origlat-swing)]));   %origlat-3.;
maxlat = ceil(max([max(stlaOrig) (origlat+swing)]));    %origlat+3.;

[lats,lons,depth] = meshgrid(lats,lons,depth);
lats = lats(:);
lons = lons(:);
depth = depth(:);

%
if threeDFlag
    %save original data
    depth_orig = depth;
    lons_orig = lons;
    lats_orig = lats;
    
    depth_ = (minDepth:maxDepth)';
    depth = repmat(depth_,size(lats));
    lats = kron(lats,ones(size(depth_)));
    lons = kron(lons,ones(size(depth_)));
    depth = depth(:);
end
clear snr_ lonEC latEC fxx fI d %S S2 S1

%%
% [templatesOrig,stnmsOrig,stlaOrig,stloOrig,refOrig,fsOrig] = distillStruct(S);
% [templates,stnms,stla,stlo,ref,fs] = distillStruct(S);
% ns = length(stlaOrig);
%
% PARENT_DIR = datestr(now,29);
% if ~exist(PARENT_DIR,'dir')
%     SUCCESSID = mkdir(PARENT_DIR);
% end
% 
% %
% secDurOrig = secDur(1);
% [T0,TH,secDur] = meshgrid(t0,th,secDurOrig);
% T0 = T0(:);
% TH = TH(:);
% secDur = secDur(:);
% 
% %
% k = 1;
% i = 1;
% 
% filtflag = true;
% dist_ = distance(origlat,origlon,stlaOrig,stloOrig,refEllipse);
% dist_ = dist_/1000;
% 
% fI = dist_ < maxDist & dist_ >= minDist;
% fIorig = fI;
% fI = find(fI);
% fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
% fI = fI';
% fI = fI(:);
% 
% 
% templates = templatesOrig(:,fI);
% stla = stlaOrig(fIorig);
% stlo = stloOrig(fIorig);
% stnms = stnmsOrig(fIorig);
% 
% threeDFlagOrig = threeDFlag;
% threeDFlag = false;
% 
% [synth,m_,mw(i,1),Gbig_tmp,observables,error_tmp(k,i),err_lb,err_ub,dur_,dist_] = ...
%     tensor_inversion_4(keepI,origlat,origlon,origdepth,templates,stla,stlo,stnms,fs,...
%     t0(1),filtflag,th(1),lfc,hfc,flag3,true,true,...
%     {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,false,npoles,...
%     secDur(1),units);
% 
% dc(k,i) = percentDC(m_(k,:)');
% if dcFlag
%     [maxErr,maxErrI] = max(error_tmp(:,i).*dc(:,i));
% else
%     [maxErr,maxErrI] = max(error_tmp(:,i));
% end
% 
% error(i) = maxErr;
% errLat(i) = lats(maxErrI);
% errLon(i) = lons(maxErrI);
% errDepth(i) = depth(maxErrI);
% synthBig(1:length(synth),i) = synth;
% obsBig(1:length(synth),i) = observables;
% mBig(:,i) = m_';
% durations(1:length(dur_),i) = dur_;
% dists(1:length(dist_),i) = dist_;
% synth2 = reshape(synth,[length(synth)/(ns) ns]);
% obs2 = reshape(observables,[length(observables)/(ns) ns]);
% normalizedAmpDiff = (max(abs(synth2)) - max(abs(obs2)))./max(abs(obs2));
% normalizedAmpDiff2 = (max(abs(obs2)) - max(abs(synth2)))./max(abs(synth2));
% 
% fI = true(size(normalizedAmpDiff));
% %fI = normalizedAmpDiff < ampThresh & normalizedAmpDiff2 < ampThresh; %true(size(normalizedAmpDiff)); %
% 
% %%
% ns = sum(fI);
% fI = fI';
% fIorig = fI;
% fI = find(fI);
% fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
% fI = fI';
% fI = fI(:);
% templatesOrig = templatesOrig(:,fI);
% stlaOrig = stlaOrig(fIorig);
% stloOrig = stloOrig(fIorig);
% stnmsOrig = stnmsOrig(fIorig);
% keepIorig = keepI(:,fIorig);
% 
% %%
% S = S(fI);
% threeDFlag = threeDFlagOrig;
% 
% %%
% error_tmp = NaN(length(lons),length(T0));
% if flag3
%     ncomps = 3;
% else
%     ncomps = 1;
% end
% GBig = zeros(fs*ns*ncomps*secDurOrig(1),5,length(T0));
% synthBig = zeros(fs*ns*ncomps*secDurOrig(1),length(T0));
% obsBig = synthBig;
% mBig = zeros(6,length(T0));
% durations = zeros(ns,length(T0));
% dists = durations;
% 
% %%
% tic;
% for i = 1:length(T0)
%     disp([num2str(i),'/', num2str(length(T0)),' ' num2str(T0(i)), ' ', num2str(TH(i))])
%     t0 = T0(i);
%     th = TH(i);
%     secDurTmp = secDur(i);
%     
%     disp(length(lats))
%     parfor k = 1:length(depth)
%         disp(k)
%         dist_ = distance(lats(k),lons(k),stlaOrig,stloOrig,refEllipse)*1e-3;
%         fI = dist_ <= maxDist & dist_ >= minDist;
%         fIorig = fI;
%         fI = find(fI);
%         fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
%         fI = fI';
%         fI = fI(:);
%         templates = templatesOrig(:,fI);
%         stla = stlaOrig(fIorig);
%         stlo = stloOrig(fIorig);
%         stnms = stnmsOrig(fIorig);
%         keepI = keepIorig(:,fIorig);
% 
%         [~,m_tmp(k,:),mw,~,~,error_tmp(k,i),err_lb,err_ub,duration] = ...
%             tensor_inversion_4(keepI,lats(k),lons(k),depth(k),templates,stla,stlo,stnms,fs,...
%             t0,filtflag,th,lfc,hfc,flag3,true,false,...
%             {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,false,npoles,...
%             secDurTmp,units);
%         dc(k,i) = percentDC(m_tmp(k,:)');
%     end
% 
%     if dcFlag
%         [maxErr,maxErrI] = max(error_tmp(:,i).*dc(:,i));
%     else
%         [maxErr,maxErrI] = max(error_tmp(:,i));
%     end
%     error(i) = maxErr;
%     errLat(i) = lats(maxErrI);
%     errLon(i) = lons(maxErrI);
%     errDepth(i) = depth(maxErrI);
% 
%     if plotFlag
%         if dcFlag
%             [~,mi] = max(error_tmp(:,i).*dc(:,i));
%         else
%             [~,mi] = max(error_tmp(:,i));
%         end
% 
%         dist_ = distance(lats(mi),lons(mi),stlaOrig,stloOrig,refEllipse);
%         dist_ = dist_/1000;
%         fI = dist_ <= maxDist & dist_ >= minDist;
%         fIorig = fI;
%         fI = find(fI);
%         fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
%         fI = fI';
%         fI = fI(:);
%         templates = templatesOrig(:,fI);
%         stla = stlaOrig(fIorig);
%         stlo = stloOrig(fIorig);
%         stnms = stnmsOrig(fIorig);
%         keepI = keepIorig(:,fIorig);
% 
%         [synth,m_,mw(i,1),Gbig_tmp,observables,err(i),elb(i),eub(i),dur,dist,azs,bazs] = ...
%             tensor_inversion_4(keepI,lats(mi),lons(mi),depth(mi),templates,stla,...
%             stlo,stnms,fs,t0,filtflag,th,lfc,hfc,flag3,true,plotFlag,...
%             {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,...
%             false,npoles,secDurTmp,units);
%         GBig(1:length(Gbig_tmp),:,i) = Gbig_tmp;
%         synthBig(1:length(synth),i) = synth;
%         obsBig(1:length(synth),i) = observables;
%         mBig(:,i) = m_';
%         durations(1:length(dur),i) = dur;
%         dists(1:length(dist),i) = dist;
%         [dc_(i),clvd_(i)] = percentDC(mBig(:,i));
%         disp(['Best Mw: ',num2str(mw(i))]);
%     end
% 
%     disp(['Best Lon.: ',num2str(lons(maxErrI))]);
%     disp(['Best Lat.: ',num2str(lats(maxErrI))]);
%     disp(['Best Depth.: ',num2str(depth(maxErrI))]);
%     disp(['Best Variance Reduction: ',num2str(error_tmp(maxErrI,i))]);
%     disp(['Best DC: ',num2str(dc(maxErrI,i))]);
%     disp(['Best DC-Scaled Variance Reduction: ',num2str(error_tmp(maxErrI,i).*dc(maxErrI,i)/100)]);
% 
%     [strike,dip,rake] = mt2sdr(mBig(:,i)');
%     [strike_a,dip_a,rake_a] = auxplane([strike, dip, rake]);
%     disp('nodal plane 1');
%     disp(round([strike dip rake]));
%     disp('nodal plane 2');
%     disp(round([strike_a dip_a rake_a]));
% end
% toc;
% save([eventID(1:end-4),'_emts'],'-v7.3');
% 
% %%
% close all;
% plot_best_solution;
