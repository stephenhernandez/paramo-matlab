clear; close all; clc;
cd ~;
%eventID = 'igepn2017wnrz.txt';
eventID = 'igepn2022ofrj.txt'; %manta 22 jul 2022 (M4.8)
%eventID = 'igepn2022oloz.txt'; %chiles/potrerillos 25 jul 2022 (M5.7)
%eventID = 'igepn2022ptut.txt';

%eventID = 'igepn2016hwwy.txt';
%eventID = 'igepn2016hwwk.txt';
%eventID = 'igepn2019druq.txt';
%eventID = 'igepn2019drun.txt';
%eventID = 'igepn2020tkhf.txt'; %pisayambo event, CCPP, 03 oct 2020
%eventID = 'igepn2021wyvu.txt'; %puembo 23 nov. 2021
%eventID = 'igepn2019drvh.txt';

%
flag3 = true;
maxDist = 450;
minDist = 50;
if flag3
    components = {'E','N','Z'}; %read 3 components no matter what
else
    components = {'Z'};
end
ncomps = length(components);

%
yyyy = eventID(6:9);
E = readSCBulletin(string(eventID));
nPphases = E.nPphases;
Pphases = E.Pphases;

%
refEllipse = referenceEllipsoid('wgs84');
origmag = E.mag;
origdepth = E.depth;
origlat = E.lat;
origlon = E.lon;
origt = E.t;

% get station lats and lons before you even attempt to read any seismic traces
stnmtmp = pull(Pphases,'stnm');
chantmp = pull(Pphases,'chan');
ntwktmp = pull(Pphases,'ntwk');
locidtmp = pull(Pphases,'locid');
[stla,stlo] = metaDataFromStationList(stnmtmp,ntwktmp,chantmp,locidtmp);
d_ = distance(origlat,origlon,stla,stlo,refEllipse)*1e-3;

%
chantmp = char(chantmp);
chantmp = string(chantmp(:,1:2));
lia = (ismember(chantmp,"HH") | ismember(chantmp,"BH") | ismember(chantmp,"HN")) ...
    & d_ < maxDist & isfinite(d_);

%
stnmtmp = stnmtmp(lia);
chantmp = chantmp(lia);
ntwktmp = ntwktmp(lia);
locidtmp = locidtmp(lia);
Pphases = Pphases(lia);
d_ = d_(lia);

[~,lia] = sort(d_);
stnmtmp = stnmtmp(lia);
chantmp = chantmp(lia);
ntwktmp = ntwktmp(lia);
locidtmp = locidtmp(lia);
Pphases = Pphases(lia);

%
lPotential = length(lia);
S = populateWaveforms(ncomps*lPotential);
S2 = S;
n = 0;
for i = 1:lPotential
    for j = 1:ncomps
        n = n + 1;
        S1 = extractWaveforms(origt-minutes(20),minutes(40),...
            string(stnmtmp(i)),string([char(chantmp(i)),components{j}]),string(ntwktmp(i)),locidtmp(i));

        if isnat(S1.ref) || strcmp(stnmtmp(i),"TOMA") || ...
                strcmp(stnmtmp(i),"IMBA") || strcmp(stnmtmp(i),"ESM1") || ...
                strcmp(stnmtmp(i),"BTER") || strcmp(stnmtmp(i),"SUCR") || ...
                strcmp(stnmtmp(i),"PIS1")
            fprintf('Skipping: %s\n',stnmtmp(i));
            continue;
        end
        S1.user0 = Pphases(i).t;
        S(n,1) = detrendWaveforms(S1);
        disp(n);
    end
end

S = padWaveforms(S);
Sorig = S;

%%
S = Sorig;
plotFlag = true;
dcFlag = false; %double-couple weighting
threeDFlag = true;
madThresh = 2.;
secDur = 300;
lfc = 1/50;
hfc = 1/20;
sta = 1/lfc; %round(1/lfc);
lta = 4*sta;
swing = 0.02;
inc = 0.01;
t0 = -1; %(-1:1)';
th = 1; %5.5;
snrThresh = 5.;     % at least...
ampThresh = 2.;     % not more than...
npoles = 4;
finalFs = 1;
units = 'vel';
tw = 0.008;
depthSwing = 60;
depthInc = 1;

%
fprintf('-------------------------------------------------------------\n');
fprintf('<strong>filtering and deconvolving instrument response</strong>\n');
fprintf('<strong>units: %s</strong>\n',units);
fprintf('-------------------------------------------------------------\n');
scaleValue = 1e6;
S = intWaveforms(...
    detrendWaveforms(...
    scaleWaveforms(...
    transferWaveforms(...
    taperWaveforms(...
    detrendWaveforms(...
    differentiateWaveforms(S)),tw),lfc,hfc,npoles,finalFs,units),scaleValue)));

refs = pull(S,'ref');
S = reshape(S,[ncomps length(S)/ncomps]);
badI = sum(isnat(reshape(refs,[ncomps length(refs)/ncomps])),1);
S(:,badI>0) = [];

S = S(:);

%
lS = length(S);
kstnm = pull(S,'kstnm');
kcmpnm = pull(S,'kcmpnm');
knetwk = pull(S,'knetwk');
khole = pull(S,'khole');

[stla,stlo] = metaDataFromStationList(kstnm,knetwk,kcmpnm,khole);
d_ = distance(stla,stlo,origlat,origlon,refEllipse)*1e-3;

% get snr, assign lats/lons
snr = NaN(lS,1);
p2p = snr;

for i = 1:lS
    d = S(i).d;
    if ~isempty(d)
        p2p(i) = max(abs(d));
    end
end

%
for i = 1:lS
    pwave = S(i).user0 - seconds(2); 
    [snr_,fxx] = freqDomSNR(S(i),pwave,secDur(1));
    fI = fxx >= lfc & fxx <= hfc;
    snr(i) = sum(snr_(fI))/sum(fI);
    S(i).stla = stla(i);
    S(i).stlo = stlo(i);
    S(i).dist = d_(i);
    S(i).evla = origlat;
    S(i).evlo = origlon;
end

goodI = snr >= snrThresh;
goodI = reshape(goodI,[ncomps length(goodI)/ncomps]);
S = reshape(S,[ncomps lS/ncomps]);
d_ = reshape(d_,[ncomps lS/ncomps]);
snr = reshape(snr,[ncomps lS/ncomps]);
p2p = reshape(p2p,[ncomps lS/ncomps]);

%

badI = sum(goodI,1) == 0;
S(:,badI) = [];
snr(:,badI) = [];
p2p(:,badI) = [];
d_(:,badI) = [];
goodI(:,badI) = [];
keepI = goodI;

S = S(:);
snr = snr(:);
p2p = p2p(:);
d_ = d_(:);

lS = length(S);
if ~lS
    disp('no good events');
    return;
end

%
close all
if flag3
    bZ = flipud(robustfit(log10(d_(3:3:end)),log10(p2p(3:3:end))));
    yZ = polyval(bZ,log10(d_(3:3:end)));

    bN = flipud(robustfit(log10(d_(2:3:end)),log10(p2p(2:3:end))));
    yN = polyval(bN,log10(d_(2:3:end)));

    bE = flipud(robustfit(log10(d_(1:3:end)),log10(p2p(1:3:end))));
    yE = polyval(bE,log10(d_(1:3:end)));

    zDiff = (p2p(3:3:end) - 10.^yZ);
    nDiff = (p2p(2:3:end) - 10.^yN);
    eDiff = (p2p(1:3:end) - 10.^yE);

    madScore = p2p;
    madScore(1:3:end) = abs((eDiff - median(eDiff))/mad(eDiff,1));
    madScore(2:3:end) = abs((nDiff - median(nDiff))/mad(nDiff,1));
    madScore(3:3:end) = abs((zDiff - median(zDiff))/mad(zDiff,1));

    figure();
    rAX(1) = subplot(311);
    loglog(d_(3:3:end),10.^yZ,'.'); title("Z");
    hold on;
    loglog(d_(3:3:end),p2p(3:3:end),'o');

    rAX(2) = subplot(312);
    loglog(d_(2:3:end),10.^yN,'.'); title("N");
    hold on;
    loglog(d_(2:3:end),p2p(2:3:end),'o');

    rAX(3) = subplot(313);
    loglog(d_(1:3:end),10.^yE,'.'); title("E");
    hold on;
    loglog(d_(1:3:end),p2p(1:3:end),'o');
    linkaxes(rAX,'x');

    figure();
    rAX(1) = subplot(311);
    semilogx(d_(3:3:end),(zDiff - median(zDiff))/mad(zDiff,1),'.'); title("Z");

    rAX(2) = subplot(312);
    semilogx(d_(2:3:end),(nDiff - median(nDiff))/mad(nDiff,1),'.'); title("N");

    rAX(3) = subplot(313);
    semilogx(d_(1:3:end),(eDiff - median(eDiff))/mad(eDiff,1),'.'); title("E");
    zoom on;
    linkaxes(rAX,'x');
else
    bZ = flipud(robustfit(log10(d_),log10(p2p)));
    yZ = polyval(bZ,log10(d_));
    madScore = (p2p - 10.^yZ); %predicted minus observed
    madScore = abs((madScore - median(madScore))/mad(madScore,1));
end

%
goodI = snr >= snrThresh & p2p <= 50*scaleValue & d_ < maxDist & d_ >= minDist & madScore < madThresh; % & ~isnat(pull(S,'ref'));
goodI = reshape(goodI,[ncomps length(goodI)/ncomps]);
S = reshape(S,[ncomps lS/ncomps]);
d_ = reshape(d_,[ncomps lS/ncomps]);

%
badI = sum(goodI,1) == 0;
S(:,badI) = [];
d_(:,badI) = [];
goodI(:,badI) = [];
keepI = goodI;
S = S(:);
d_ = d_(:);
lS = length(S);

%
if ~lS
    disp('no good events');
    return;
end

% cut data
tref = min(t0);
S = cutWaveforms(S,origt+seconds(tref),0,secDur*2);
S = syncWaveforms(S);

%
lats = (origlat-swing):inc:(origlat+swing);
lons = (origlon-swing):inc:(origlon+swing);

depth = origdepth;
depth = floor(depth)+0.5;

minDepth = depth - depthSwing;
minDepth = max([minDepth 0.5]);
maxDepth = depth + depthSwing;
maxDepth = min([maxDepth 149.5]);

minlon = floor(min([-81.5 min([min(stlo) (origlon-swing)])]));   %origlon-3.;
maxlon = ceil(max([max(stlo) (origlon+swing)]));    %origlon+3.;
minlat = floor(min([min(stla) (origlat-swing)]));   %origlat-3.;
maxlat = ceil(max([max(stla) (origlat+swing)]));    %origlat+3.;

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

    depthTmp = (minDepth:depthInc:maxDepth)';
    sizeDepth = size(depthTmp);
    depth = repmat(depthTmp,size(lats));
    lats = kron(lats,ones(sizeDepth));
    lons = kron(lons,ones(sizeDepth));
    depth = depth(:);
end

%
PARENT_DIR = datestr(now,29);
if ~exist(PARENT_DIR,'dir')
    SUCCESSID = mkdir(PARENT_DIR);
end

%
secDurOrig = secDur(1);
[T0,TH,secDur] = meshgrid(t0,th,secDurOrig);
T0 = T0(:);
TH = TH(:);
secDur = secDur(:);

filtflag = true;
clear synth m_ mw Gbig_tmp observables error_tmp err_lb err_ub dur_ dist_ error errLon errLat errDepth

%
[templates,stnms,stla,stlo,ref,fs] = distillStruct(S,ncomps);
ns = length(stla);
if isfinite(lfc)
    templates = detrend(templates);
end

% preallocate
error_tmp = NaN(length(lons),length(T0));
GBig = zeros(fs*ns*ncomps*secDurOrig(1),5,length(T0));
synthBig = zeros(fs*ns*ncomps*secDurOrig(1),length(T0));
obsBig = synthBig;
mBig = zeros(6,length(T0));
durations = zeros(ns,length(T0));
dists = durations;

%%
pwaves = pull(S,'user0');
ttp = seconds(pwaves - origt);
if flag3
    ttp = ttp(3:3:end);
end
ttp = zeros(size(ttp)); round(ttp*fs)/fs;

%
tic;
for i = 1:length(T0)
    disp([num2str(i),'/', num2str(length(T0)),' ' num2str(T0(i)), ' ', num2str(TH(i))])
    t0 = T0(i);
    th = TH(i);
    secDurTmp = secDur(i);

    disp(length(lats))
    lenDepth = length(depth);
    parfor k = 1:lenDepth
        lat_ = lats(k);
        lon_ = lons(k);
        depth_ = depth(k);
        fprintf('%d/%d %f %f %f\n',k,lenDepth,lat_,lon_,depth_);

        [~,m_tmp(k,:),mw,~,~,error_tmp(k,i),err_lb,err_ub,duration] = ...
            tensor_inversion_4(keepI,lat_,lon_,depth_,templates,stla,stlo,stnms,fs,...
            t0,tref,ttp,filtflag,th,lfc,hfc,flag3,true,false,...
            {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,false,npoles,...
            secDurTmp,units);
        dc(k,i) = percentDC(m_tmp(k,:)');
    end

    if dcFlag
        [maxErr,maxErrI] = max(error_tmp(:,i).*dc(:,i));
    else
        [maxErr,maxErrI] = max(error_tmp(:,i));
    end
    error(i,1) = maxErr;
    errLat(i,1) = lats(maxErrI);
    errLon(i,1) = lons(maxErrI);
    errDepth(i,1) = depth(maxErrI);

    if plotFlag
        if dcFlag
            [~,mi] = max(error_tmp(:,i).*dc(:,i));
        else
            [~,mi] = max(error_tmp(:,i));
        end

        [synth,m_,mw(i,1),Gbig_tmp,observables,err(i),elb(i),eub(i),dur,dist,azs,bazs] = ...
            tensor_inversion_4(keepI,lats(mi),lons(mi),depth(mi),templates,stla,...
            stlo,stnms,fs,t0,tref,ttp,filtflag,th,lfc,hfc,flag3,true,plotFlag,...
            {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,...
            false,npoles,secDurTmp,units);
        GBig(1:length(Gbig_tmp),:,i) = Gbig_tmp;
        synthBig(1:length(synth),i) = synth;
        obsBig(1:length(synth),i) = observables;
        mBig(:,i) = m_';
        durations(1:length(dur),i) = dur;
        dists(1:length(dist),i) = dist;
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
    disp('nodal plane 1');
    disp(round([strike dip rake]));
    disp('nodal plane 2');
    disp(round([strike_a dip_a rake_a]));
end
toc;

%
save([eventID(1:end-4),'_emts'],'-v7.3');

close all;
plot_best_solution;
