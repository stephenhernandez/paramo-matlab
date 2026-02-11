function [synth,m,mw,G,observables,vr,vr_l,vr_u,duration,dist,azs,bazs] = ...
    tensor_inversion_4(keepI,testLat,testLon,testDepth,data,stla,stlo,stnms,...
    fs,t0,tref,ttp,filtflag,th,lfc,hfc,flag3,synthflag,gmtflag,gmtopts,cFlag,wFlag,...
    npoles,secDur,units,rake)

if nargin < 23
    rake = 90;
end

t0shift = t0;
tdum = tref + (0:size(data,1)-1)'/fs;
tI = tdum >= t0shift;
data = data(tI,:);
ns = length(stla);
filtFlag = false; %filter observables?
unitFlag = units;
testDepth = round(2*testDepth)/2;
stf = getSTF(th,fs);

%refVel = 4; %stupid fast velocity in km/2?
%pad = 10;

%%
refEllipse = referenceEllipsoid('wgs84');
[dist,azs] = distance(testLat,testLon,stla,stlo,refEllipse);
[~,bazs] = distance(stla,stlo,testLat,testLon,refEllipse);
dist = 1e-3*dist;

%[dist_,mdi] = max(dist);
% tnew = floor(dist_ / refVel) - pad;
% tnew = max([0 tnew]);
% tnew = tnew*fs + 1;

%convert to samples
if flag3
    flagScalar = 3;
    %tmp1 = data(tnew:end,1+(mdi-1)*3);
else
    flagScalar = 1;
    %tmp1 = data(tnew:end,mdi);
end
duration_ = secDur*fs+1; %size(data,1);

duration = ones(size(stla))*duration_;
endTry = duration_;
observables = zeros(endTry*3*ns,1);

wnew = ones(endTry*3*ns,1);
G = repmat(observables,1,5);
enz = zeros(endTry,3);

% if flag3
%     flagScalar = 3;
% else
%     flagScalar = 1;
% end
totlength = 0;

for i = 1:ns
    %dist_ = dist(i);
    %tnew = floor(dist_ / refVel) - pad; %4 km/s
    tnew = ttp(i)-4*fs;
    tnew = max([0 tnew]);
    tnew = tnew*fs+1;

    %cut the observable data here
    stnm = stnms{i};
    if flag3
        try
            enz(:,1) = data(tnew:tnew+duration_-1,1+(i-1)*3);
            enz(:,2) = data(tnew:tnew+duration_-1,2+(i-1)*3);
            enz(:,3) = data(tnew:tnew+duration_-1,3+(i-1)*3);
        catch ME
            disp(i);
            disp(tnew)
            disp(tnew+duration_)
            rethrow(ME);
        end
    else
        enz = data(tnew:tnew+duration_-1,i);
    end

    [observable,~,duration(i)] = readObservable2(stnm,enz,lfc,hfc,npoles,filtFlag,flag3,fs,0,0,secDur);
    dur2 = flagScalar*duration(i);
    keepI_ = keepI(:,i);
    if flagScalar == 1
        keepI_ = keepI_(end);
    end
    wnew_ = repmat(keepI_,1,duration(i));
    wnew_ = wnew_';
    
    observables(1+(i-1)*dur2:dur2*i) = observable;
    wnew(1+(i-1)*dur2:dur2*i) = wnew_(:);
    
    %g = createGF2(bazs(i),azs(i),dist(i),filtflag,lfc,hfc,npoles,stf,(tnew-1)/2,duration(i),testDepth,flag3,synthflag,fs,stnm,unitFlag);
    g = createGF2(bazs(i),azs(i),dist(i),filtflag,lfc,hfc,npoles,stf,tnew,duration(i),testDepth,flag3,synthflag,fs,stnm,unitFlag);
    G(1+(i-1)*dur2:dur2*i,:) = g;
    totlength = totlength + length(observable);
end

%%
observables = observables(1:totlength);
G = G(1:totlength,:);
wnew = wnew(1:totlength);
goodI = find(wnew);

%% perform inversion
if cFlag %if constraining a solution...
    strike = 307;
    if testDepth <= 10
        dip = 6;
    elseif testDepth <= 15
        dip = 10;
    elseif testDepth <= 38
        dip = 25;
    else
        dip = 43;
    end
    % dip = 20;
    sdr = [strike,dip,rake];
    mt = sdr2mt(sdr);
    synth = G*mt([1 2 4 5 6])';
    observables = observables/norm(observables);
    synth = synth/norm(synth);
    %err = 100*dot(synth,observables)/(norm(synth)*norm(observables));
    [vr,vr_l,vr_u] = varRed(observables,synth,10);
    m = mt;
else
    if flag3
        if wFlag
            Ntot = length(G);
            Nsub = Ntot/(3*ns);
            W = ones(Nsub*3,1);
            W(1:Nsub*2) = 0.2;
            W = repmat(W,ns,1);
            m = lscov(G,observables,W);
            synth = G*m;
            [vr,vr_l,vr_u] = varRed(observables,synth,10);
            m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
        else
            m = lscov(G,observables,wnew);
            synth = G*m;
            synth2 = synth(goodI);
            obs2 = observables(goodI);

            [vr,vr_l,vr_u] = varRed(obs2,synth2,20);
            m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
        end
    else
        m = ((G'*G)^-1)*G'*observables;
        synth = G*m;
        [vr,vr_l,vr_u] = varRed(observables,synth,10);
        m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
    end
end

%%
if gmtflag
    disp('writing file')
    PARENT_DIR = gmtopts{7};
    disp(['m = ', num2str(m)])
    t0str = num2str(abs(t0));
    thstr = num2str(th);
    durStr = num2str(secDur);
    gmt_title = strcat(['t0 = ',t0str,'; lat. = ',num2str(testLat),'; lon. = ',num2str(testLon)]);
    evlo = testLon;
    evla = testLat;
    evdp = testDepth;
    expon = 32;
    lonb = evlo;
    latb = evla;
    
    date_string = datestr(now,'yyyymmddTHHMMSSFFF');
    latStr = num2str(round(abs(testLat*100)));
    lonStr = num2str(round(abs(testLon*100)));
    depthStr = num2str(round(abs(testDepth*1000)));
    mecfile = ['psmeca_', date_string,'_',t0str,'T0_',thstr,'TH_',durStr,'Dur_',latStr,'LAT_',lonStr,'LON_',depthStr,'DEPTH.txt'];
    script_name = [strcat(PARENT_DIR,'/tensor_'), date_string,'_',t0str,'T0_',thstr,'TH_',durStr,'Dur_',latStr,'LAT_',lonStr,'LON_',depthStr,'DEPTH.sh'];
    WEST = gmtopts{1};
    EAST = gmtopts{2};
    SOUF = gmtopts{3};
    NORF = gmtopts{4};
    PROJECTION = gmtopts{5};
    CLSTR = gmtopts{6};
    
    
    dlmwrite(script_name,'#!/bin/bash','delimiter','');
    dlmwrite(script_name,['outfile="',mecfile,'"'],'delimiter','','-append');
    dlmwrite(script_name,['gmt psbasemap -B1g1:."',gmt_title,'": -R',num2str(WEST),'/',num2str(EAST),'/',num2str(SOUF),'/',num2str(NORF),' -J',PROJECTION ,' -K > ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(script_name,'gmt pscoast -R -J -Df -W2 -K -O -Ggray >> ${outfile}.ps','delimiter','','-append');
    dlmwrite(script_name,['gmt psmeca ',mecfile, ' -R -J -Sm0.75 -O -h1 -K -G',CLSTR,' >> ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(script_name,['gmt psmeca ',mecfile, ' -R -J -Sd0.75 -O -h1 -T0 >> ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(strcat(PARENT_DIR,'/',mecfile),'evlo evla evdp mrr mtt mpp mrt mrp mtp exponent lonb latb','delimiter','');
    dlmwrite(strcat(PARENT_DIR,'/',mecfile),[num2str(evlo),' ',num2str(evla),' ',num2str(evdp),' ',num2str(m),' ',num2str(expon),' ',num2str(lonb),' ',num2str(latb)],'delimiter','','-append');
    
end
mw = mt2mw(m);
