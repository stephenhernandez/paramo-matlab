%function emts(eventID)
clear;
close all;
more off;
% %cd ~/regions/pichincha
cd ~/;
%regions/sierra_negra/
%eventID = 'igepn2016hnmu.txt'; %7.8 mainshock
%eventID = 'igepn2016htqh.txt'; %2016 04 20  08 33 42 (2 minutes until next big event)
%eventID = 'igepn2018efav.txt'; %
%eventID = 'igepn2016htqj.txt'; %2016 04 20  08 35 08
%eventID = 'igepn2017wnrz.txt'; %2017-11-17, guayaquil earthquake
%eventID = 'igepn2019dcwv.txt';
%eventID = 'igepn2016hwwy.txt'; %22-Apr-2016 03:20, M6
%eventID = 'igepn2018mkgw.txt'; %5.3 from sierra negra
%eventID = 'igepn2016hobd.txt'; %2016 04 17	07 13 57, M6 (edited sc3 file...)
%eventID = 'igepn2016hofk.txt'; %2016 04 17	09 23 39, M5.7 (edited sc3 file...)
%eventID = 'igepn2016hwwk.txt'; %22-Apr-2016 03:03:39 6.1603 (15 minutes until next big event)
%eventID = 'igepn2016ifpw.txt'; %26-Apr-2016 21:58
%eventID = 'igepn2016hswb.txt';
%eventID = 'igepn2016hyoq.txt';
%eventID = 'igepn2016hnmk.txt'; %5.1 foreshock, best with vertical displacement records (M5.1)
%eventID = 'igepn2015qmro.txt'; %normal earthquake on the outer trench, not enough stations for inversion
%eventID = 'igepn2015qphv.txt'; %normal earthquake on the outer trench, SV solution available (M4)
%eventID = 'igepn2016fgfa.txt'; %quito earthquake, march 15, 2016
%eventID = 'igepn2016ksml.txt'; %normal faulting aftershock, best with vertical displacement records (M5)
%eventID = 'igepn2016jtkp.txt'; % second may 18 aftershock
%eventID = 'igepn2016lcdd.txt'; %june 6 aftershock
%eventID = 'igepn2016nedp.txt'; %july 6, esmeraldas 4.9 event
%eventID = 'igepn2016nhqy.txt'; %july 8, 5.5 near trench event
%eventID = 'igepn2016nmyn.txt'; %july 10, 5.9, first of july 10 doublet
%eventID = 'igepn2016ntoa.txt'; %off-shore esmeraldas event, 4.7, 14 july
%eventID = 'igepn2016odit.txt'; %5.3 near manta, evening of 19 july
%eventID = 'igepn2016pocw.txt'; %4.6 in quito, 08 Aug2016
%eventID = 'igepn2016yvmv.txt';
%eventID = 'igepn2017cdxo.txt'; %5.3 near esmeraldas
%eventID = 'igepn2016hpdo.txt'; %17-APR manta event (with 6.4 mag that is way too high)
%eventID = 'igepn2017kdtr.txt'; %pichincha event, may 2017
%eventID = 'igepn2014psys.txt'; %calderon, 2014, 5.1
eventID = 'igepn2021wyvu.txt'; %puembo/quito 23 nov 2021

%%
yyyy = eventID(6:9);
%['~/phaseInformationSC3/',yyyy,'/',eventID]
%E = readSCBulletin(string(['~/phaseInformationSC3/',yyyy,'/',eventID]));
E = readSCBulletin(string(eventID));
plotFlag = true;
dcFlag = true; %double-couple weighting (not "direct current"!!)
flag3 = true;
threeDFlag = true;

%%
origmag = E.mag;
origdepth = E.depth;
origlat = E.lat;
origlon = E.lon;

t = E.t;
[yyyy,mm,days] = datevec(t);
npoles = 4;
finalFs = 1;
units = 'disp';

components = {'E','N','Z'}; %read 3 components no matter what
ncomps = length(components);
origtime = E.t; %OabsTime;
dfrac = seconds(origtime - datetime(yyyy,mm,days))/86400; %floor(origtime);
refEllipse = referenceEllipsoid('wgs84');

if dfrac >= 0.993 %~10 minutes til end of day, might as well get next days data...
    days = [days days+1];
end

if origmag >= 3.
    noise = 600;
    signal = 600;
    secDur = 60;
    lfc = 1/32;
    hfc = 1/8;
    sta = 1/lfc; %round(1/lfc);
    lta = 4*sta;
    swing = 0.01;
    inc = 0.01;
    maxDist = 500;
    minDist = 0;
    t0 = 0;
    th = 0.:0.5:5.5;
    snrThresh = 1.; %at least...
    ampThresh = 1.;  % not more than...
    npoles = 4;
    finalFs = 4;
    units = 'vel';
    tw = 0.01;
end
nPphases = E.nPphases;
Pphases = E.Pphases;

%% get station lats and lons before you even attempt to read any seismic traces
bbdata = load('~/igdata/ecuador_bb_meta_data');
bbkstnm = bbdata.bbkstnm;
bbstla = bbdata.bbstla;
bbstlo = bbdata.bbstlo;

masterStnms = cell(nPphases,1);
masterStla = NaN(nPphases,1);
masterStlo = masterStla;
stnmSubset = false(nPphases,1);

for i = 1:nPphases
    stnmtmp = Pphases(i).stnm;
    [lia,locb] = ismember(stnmtmp,bbkstnm);
    if lia
        stnmSubset(i) = true;
        masterStnms{i} = stnmtmp;
        masterStla(i) = bbstla(locb);
        masterStlo(i) = bbstlo(locb);
    else
        disp(['Will not process: ',stnmtmp]);
    end
end

masterStnms = masterStnms(stnmSubset);
masterStla = masterStla(stnmSubset);
masterStlo = masterStlo(stnmSubset);
Pphases = Pphases(stnmSubset);

%% get events greater than minDist and nearer than maxDist
newn = sum(stnmSubset);
dep = zeros(newn,1);
az = dep;
for i = 1:newn
    [dep(i),az(i)] = distance(origlat,origlon,masterStla(i),masterStlo(i),refEllipse);
end
dep = dep/1000;
nI = dep >= minDist & dep <= maxDist;
masterStnms = masterStnms(nI);
masterStla = masterStla(nI);
masterStlo = masterStlo(nI);
Pphases = Pphases(nI);
newn = sum(nI);

%% get (and cut) data
tic;
n = 0;
S = populateWaveforms(1);
newn = 50;
for i = 1:newn
    stnmtmp = Pphases(i).stnm;
    if strcmp(stnmtmp,'SUCR') || strcmp(stnmtmp,'SLOR') || strcmp(stnmtmp,'YAHU') ...
            || strcmp(stnmtmp,'TOMA') || strcmp(stnmtmp,'BRRN') ...
            || strcmp(stnmtmp,'CAB1') || strcmp(stnmtmp,'SRAM') ...%|| strcmp(stnmtmp,'MILO')...
            || strcmp(stnmtmp,'CABP') || strcmp(stnmtmp,'POND') || strcmp(stnmtmp,'URCU')...
            || strcmp(stnmtmp,'ISPG') || strcmp(stnmtmp,'TULM') || strcmp(stnmtmp,'BONI') || strcmp(stnmtmp,'CHMA') %|| strcmp(stnmtmp,'PTGL')
        %|| strcmp(stnmtmp,'POND') || strcmp(stnmtmp,'URCU')
        %strcmp(stnmtmp,'YAHU') || strcmp(stnmtmp,'POND') || strcmp(stnmtmp,'URCU')...
        %|| strcmp(stnmtmp,'BMOR') || strcmp(stnmtmp,'VCES') || strcmp(stnmtmp,'SLOR')...
        %|| strcmp(stnmtmp,'BBIL') || strcmp(stnmtmp,'BVC2') || strcmp(stnmtmp,'BTAM')...
        %|| strcmp(stnmtmp,'BRUN') || strcmp(stnmtmp,'BRRN') || strcmp(stnmtmp,'TOMA')...
        %|| strcmp(stnmtmp,'SRAM') || strcmp(stnmtmp,'ANTS') || strcmp(stnmtmp,'BMAS')...
        %|| strcmp(stnmtmp,'BNAS') || strcmp(stnmtmp,'CHL2') || strcmp(stnmtmp,'CUSE')...
        %|| strcmp(stnmtmp,'CUSW') || strcmp(stnmtmp,'LNGL') || strcmp(stnmtmp,'CHMA')...
        %|| strcmp(stnmtmp,'TAMH') || strcmp(stnmtmp,'ANTM') || strcmp(stnmtmp,'CHSH')...
        %|| strcmp(stnmtmp,'BULB') || strcmp(stnmtmp,'BOSC') || strcmp(stnmtmp,'ARNL')...
        %|| strcmp(stnmtmp,'CAYA') || strcmp(stnmtmp,'CAB1') || strcmp(stnmtmp,'MILO')...
        %|| strcmp(stnmtmp,'TULM') || strcmp(stnmtmp,'SUCR')
        disp(' ')
        disp(['cant find response info for ',stnmtmp,', moving on...'])
        disp(' ')
        continue;
    else
        for j = 1:ncomps
            comptmp = components{j};
            %Stmp = get_bb_data(stnmtmp,comptmp,yyyy,mm,days,units,finalFs,lfc,hfc,npoles,origtime,noise,signal);
            
            %Stmp = loadWaveforms(datetime(yyyy,mm,days),1,stnmtmp,comptmp,"EC","");
            chan = ['HH',comptmp;...
                'BH',comptmp];
            Stmp = extractWaveforms(datetime(origtime)-seconds(noise),seconds(noise+signal),...
                stnmtmp,chan,"EC","",true,false);
            Stmp = Stmp(1);
            direction = 1;
            waFlag = false;
            Stmp = differentiateWaveforms(Stmp);
            Stmp = transferWaveforms(Stmp,lfc,hfc,npoles,finalFs,units,direction,waFlag);
            if isempty(Stmp)
                continue;
            else
                n = n+1;
                S(n) = Stmp;
            end
        end
    end
end
toc;

%% synch data
if n
    S = padWaveforms(S);
end
Sorig = S;

%%
%[~,stnms,stla,stlo] = distillStruct(S);
stnms = pull(S,'kstnm');
[stla,stlo,stelev] = metaDataFromStationList(stnms);
ns = length(stnms);
fs = finalFs;

%% filter based on SNR ratio
snr_e = zeros(ns,1);
snr_n = snr_e;
snr_z = snr_e;
ns = round(length(S)/3);
trueI = true(ns,3);
[~,azs] = distance(origlat,origlon,stla,stlo,refEllipse);
hFlag = true;

%
for i = 1:ns
    disp(i);
    S_ = populateWaveforms(3);
    S_(1) = S(1+(i-1)*3);
    S_(2) = S(2+(i-1)*3);
    S_(3) = S(3+(i-1)*3);
    S_ = rotateWaveforms(S_,azs(i));
    %[locs,snr,staOlta,sosSTA] = stalta(S,sta,lta,mph,hFlag,plotFlag,envFiltFlag,hfc,verboseFlag);
    [~,snrtmpe] = stalta(S_(1),sta,lta,-100,hFlag,false,false);
    snr_e(i) = max(snrtmpe);
    [~,snrtmpn] = stalta(S_(2),sta,lta,-100,hFlag,false,false);
    snr_n(i) = max(snrtmpn);
    [~,snrtmpz] = stalta(S_(3),sta,lta,-100,hFlag,false,false);
    snr_z(i) = max(snrtmpz);
end
snr_ = [snr_e'; snr_n'; snr_z'];
keepI = snr_ >= -999; %snrThresh;
snr_ = min(snr_);
clear de dn dz azs

%% cut data
t = (S(1).ref + (0:S(1).npts-1)*S(1).delta/86400)';
tstart = find(t >= origtime,1);
tend = tstart - 1 + secDur/S(1).delta;
tdum = t(tstart:tend);
newN = length(tdum);
newref = t(tstart);

for i = 1:length(S)
    dtmp = S(i).d;
    dtmp = dtmp(tstart:tend);
    S(i).d = dtmp;
    S(i).ref = newref;
    S(i).npts = newN;
end
clear dtmp newN

%%
[templatesOrig,stnmsOrig,stlaOrig,stloOrig,refOrig,fsOrig] = distillStruct(S);
[templates,stnms,stla,stlo,ref,fs] = distillStruct(S);
ns = length(stlaOrig);

%%
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

%%
if threeDFlag
    %save original data
    depth_orig = depth;
    lons_orig = lons;
    lats_orig = lats;
    
    depth_ = (minDepth:maxDepth)';
    %depth_ = (minDepth:2:maxDepth)';
    %depth_ = (5.5:29.5)';
    depth = repmat(depth_,size(lats));
    lats = kron(lats,ones(size(depth_)));
    lons = kron(lons,ones(size(depth_)));
    depth = depth(:);
end

PARENT_DIR = datestr(now,29);
if ~exist(PARENT_DIR,'dir')
    SUCCESSID = mkdir(PARENT_DIR);
end

%%
secDurOrig = secDur;
[T0,TH,secDur] = meshgrid(t0,th,secDurOrig);
T0 = T0(:);
TH = TH(:);
secDur = secDur(:);

%%
k = 1;
i = 1;

filtflag = true;
dist_ = distance(origlat,origlon,stlaOrig,stloOrig,refEllipse);
dist_ = dist_/1000;

fI = dist_ < maxDist & dist_ >= minDist;
fIorig = fI;
fI = find(fI);
fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
fI = fI';
fI = fI(:);

%%
templates = templatesOrig(:,fI);
stla = stlaOrig(fIorig);
stlo = stloOrig(fIorig);
stnms = stnmsOrig(fIorig);

threeDFlagOrig = threeDFlag;
threeDFlag = false;

[synth,m_,mw(i,1),Gbig_tmp,observables,error_tmp(k,i),err_lb,err_ub,dur_,dist_] = ...
    tensor_inversion_4(keepI,origlat,origlon,origdepth,templates,stla,stlo,stnms,fs,...
    t0(1),filtflag,th(1),lfc,hfc,flag3,true,true,...
    {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,false,npoles,...
    secDur(1),units);

dc(k,i) = percentDC(m_(k,:)');
if dcFlag
    [maxErr,maxErrI] = max(error_tmp(:,i).*dc(:,i));
else
    [maxErr,maxErrI] = max(error_tmp(:,i));
end

error(i) = maxErr;
errLat(i) = lats(maxErrI);
errLon(i) = lons(maxErrI);
errDepth(i) = depth(maxErrI);
synthBig(1:length(synth),i) = synth;
obsBig(1:length(synth),i) = observables;
mBig(:,i) = m_';
durations(1:length(dur_),i) = dur_;
dists(1:length(dist_),i) = dist_;
synth2 = reshape(synth,[length(synth)/(ns) ns]);
obs2 = reshape(observables,[length(observables)/(ns) ns]);
normalizedAmpDiff = (max(abs(synth2)) - max(abs(obs2)))./max(abs(obs2));
normalizedAmpDiff2 = (max(abs(obs2)) - max(abs(synth2)))./max(abs(synth2));
fI = true(size(normalizedAmpDiff));
%normalizedAmpDiff < ampThresh & normalizedAmpDiff2 < ampThresh; %true(size(normalizedAmpDiff)); %

%%
ns = sum(fI);
fI = fI';
fIorig = fI;
fI = find(fI);
fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
fI = fI';
fI = fI(:);
templatesOrig = templatesOrig(:,fI);
stlaOrig = stlaOrig(fIorig);
stloOrig = stloOrig(fIorig);
stnmsOrig = stnmsOrig(fIorig);
keepIorig = keepI(:,fIorig);

%%
S = S(fI);
threeDFlag = threeDFlagOrig;

%%
error_tmp = NaN(length(lons),length(T0));
if flag3
    ncomps = 3;
else
    ncomps = 1;
end
GBig = zeros(fs*ns*ncomps*secDurOrig,5,length(T0));
synthBig = zeros(fs*ns*ncomps*secDurOrig,length(T0));
obsBig = synthBig;
mBig = zeros(6,length(T0));
durations = zeros(ns,length(T0));
dists = durations;

%%
tic;
for i = 1:length(T0)
    disp([num2str(i),'/', num2str(length(T0)),' ' num2str(T0(i)), ' ', num2str(TH(i))])
    t0 = T0(i);
    th = TH(i);
    secDurTmp = secDur(i);
    
    disp(length(lats))
    parfor k = 1:length(depth)
        disp(k)
        dist_ = distance(lats(k),lons(k),stlaOrig,stloOrig,refEllipse)*1e-3;
        fI = dist_ <= maxDist & dist_ >= minDist;
        fIorig = fI;
        fI = find(fI);
        fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
        fI = fI';
        fI = fI(:);
        templates = templatesOrig(:,fI);
        stla = stlaOrig(fIorig);
        stlo = stloOrig(fIorig);
        stnms = stnmsOrig(fIorig);
        keepI = keepIorig(:,fIorig);
        
        [~,m_tmp(k,:),mw,~,~,error_tmp(k,i),err_lb,err_ub,duration] = ...
            tensor_inversion_4(keepI,lats(k),lons(k),depth(k),templates,stla,stlo,stnms,fs,...
            t0,filtflag,th,lfc,hfc,flag3,true,false,...
            {minlon,maxlon,minlat,maxlat,'M5i','black',PARENT_DIR},false,false,npoles,...
            secDurTmp,units);
        dc(k,i) = percentDC(m_tmp(k,:)');
    end
    
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
        
        dist_ = distance(lats(mi),lons(mi),stlaOrig,stloOrig,refEllipse);
        dist_ = dist_/1000;
        fI = dist_ <= maxDist & dist_ >= minDist;
        fIorig = fI;
        fI = find(fI);
        fI = [1+(fI-1)*3 2+(fI-1)*3 3+(fI-1)*3];
        fI = fI';
        fI = fI(:);
        templates = templatesOrig(:,fI);
        stla = stlaOrig(fIorig);
        stlo = stloOrig(fIorig);
        stnms = stnmsOrig(fIorig);
        keepI = keepIorig(:,fIorig);
        
        [synth,m_,mw(i,1),Gbig_tmp,observables,err(i),elb(i),eub(i),dur,dist,azs,bazs] = ...
            tensor_inversion_4(keepI,lats(mi),lons(mi),depth(mi),templates,stla,...
            stlo,stnms,fs,t0,filtflag,th,lfc,hfc,flag3,true,plotFlag,...
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
save([eventID(1:end-4),'_emts'],'-v7.3');

%%
close all;
plot_best_solution;