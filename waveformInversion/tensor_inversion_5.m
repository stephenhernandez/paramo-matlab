function [synth,m,mw,G,observables,vr,vr_l,vr_u,duration,dist,azs,bazs] = tensor_inversion_5(testLat,testLon,testDepth,data,stla,stlo,stnms,fs,t0,filtflag,th,lfc,hfc,flag3,synthflag,gmtflag,gmtopts,cFlag,wFlag,npoles,secDur,units,rake,dip)
if nargin < 23
    rake = 90;
end

%t0shift = t0;
ns = length(stla);
filtFlag = false; %filter observables?
unitFlag = units;
testDepth = round(2*testDepth)/2;
disp([th fs])
stf = getSTF(th,fs);

refEllipse = referenceEllipsoid('wgs84');
[dist,azs] = distance(testLat,testLon,stla,stlo,refEllipse);
[~,bazs] = distance(stla,stlo,testLat,testLon,refEllipse);
dist = 1e-3*dist;
disp(max(dist))
[dist_,mdi] = max(dist);
tnew = floor(dist_ / 3.5) - 10;
tnew = max([0 tnew]);
tnew = tnew*fs + 1; %convert to samples
disp(tnew)
disp(size(data))
tmp1 = data(tnew:end,1+(mdi-1)*3);
duration_ = length(tmp1);

duration = ones(size(stla))*duration_;
endTry = duration_;
observables = zeros(endTry*3*ns,1);
wnew = ones(endTry*3*ns,1);
G = repmat(observables,1,5);
enz = zeros(endTry,3);

if flag3
    flagScalar = 3;
else
    flagScalar = 1;
end
totlength = 0;

disp([duration_ duration_])
for i = 1:ns
    dist_ = dist(i);
    tnew = floor(dist_ / 3.5) - 10; %3.5 km/s
    tnew = max([0 tnew]);
    tnew = tnew*fs+1;
    
    stnm = stnms{i};
    %cut the observable data here
    enz(:,1) = data(tnew:tnew+duration_-1,1+(i-1)*3);
    enz(:,2) = data(tnew:tnew+duration_-1,2+(i-1)*3);
    enz(:,3) = data(tnew:tnew+duration_-1,3+(i-1)*3);
    
    %[observable,t0shift,duration(i)] = readObservable2(stnm,enz,lfc,hfc,npoles,filtFlag,flag3,fs,t0shift,azs(i),secDur);
    observable = readObservable2(stnm,enz,lfc,hfc,npoles,filtFlag,flag3,fs,0,0,secDur); %azs(i),secDur);
    dur2 = flagScalar*duration(i);
    observables(1+(i-1)*dur2:dur2*i) = observable;
    %    g = createGF2(azs(i),dist(i),filtflag,lfc,hfc,npoles,stf,0,duration(i),testDepth,flag3,synthflag,fs,stnm,unitFlag);
    
    % cut the synthetic data in the function
    %g = createGF2(azs(i),dist(i),filtflag,lfc,hfc,npoles,stf,(tnew-1)/2,duration(i),testDepth,flag3,synthflag,fs,stnm,unitFlag);
    g = createGF2(bazs(i),azs(i),dist(i),filtflag,lfc,hfc,npoles,stf,(tnew-1)/2,duration(i),testDepth,flag3,synthflag,fs,stnm,unitFlag);
    G(1+(i-1)*dur2:dur2*i,:) = g;
    %     gunfilt = createGF2(azs(i),dist(i),filtflag,lfc,hfc,npoles,1,0,duration(i),testDepth,flag3,synthflag,fs,stnm);
    %     Gunfilt(1+(i-1)*dur2:dur2*i,:) = gunfilt;
    totlength = totlength + length(observable);
    %     else
    %         continue
    %     end
end
observables = observables(1:totlength);
G = G(1:totlength,:);
wnew = wnew(1:totlength);
goodI = find(wnew);
%Gunfilt = Gunfilt(1:totlength,:);

% %% get weights
% W = eye(length(observables));
% w = abs(hilbert(observables));
% w = w/max(abs(w));
% %w = w.^2;
% %w = ones(size(observables));
% W(W == 1) = 1./w;
% Winv = W;
% clear W

%% perform inversion
if cFlag
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
            %dD1 = synth - observables;
            %err = (1-(sum(dD1.^2)/sum(observables.^2)))*100;
            [vr,vr_l,vr_u] = varRed(observables,synth,10);
            m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
        else
            m = lscov(G,observables,wnew);
            %m = ((G'*Winv*G)^-1)*G'*Winv*observables;
            %m = ((G'*G)^-1)*G'*observables;
            synth = G*m;
            synth2 = synth(goodI);
            obs2 = observables(goodI);
            [vr,vr_l,vr_u] = varRed(obs2,synth2,20);
            m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
        end
    else
        %m = lscov(G,observables);
        m = ((G'*G)^-1)*G'*observables;
        synth = G*m;
        [vr,vr_l,vr_u] = varRed(observables,synth,10);
        m = [m(1) m(2) -(m(1)+m(2)) m(3) m(4) m(5)];
    end
end

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
    dlmwrite(script_name,['psbasemap -B1g1:."',gmt_title,'": -R',num2str(WEST),'/',num2str(EAST),'/',num2str(SOUF),'/',num2str(NORF),' -J',PROJECTION ,' -K > ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(script_name,'pscoast -R -J -Df -W2 -K -O -Ggray >> ${outfile}.ps','delimiter','','-append');
    dlmwrite(script_name,['psmeca ',mecfile, ' -R -J -Sm0.75 -O -h1 -K -G',CLSTR,' >> ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(script_name,['psmeca ',mecfile, ' -R -J -Sd0.75 -O -h1 -T0 >> ${outfile}.ps'],'delimiter','','-append');
    dlmwrite(strcat(PARENT_DIR,'/',mecfile),'evlo evla evdp mrr mtt mpp mrt mrp mtp exponent lonb latb','delimiter','');
    dlmwrite(strcat(PARENT_DIR,'/',mecfile),[num2str(evlo),' ',num2str(evla),' ',num2str(evdp),' ',num2str(m),' ',num2str(expon),' ',num2str(lonb),' ',num2str(latb)],'delimiter','','-append');
    
end
mw = mt2mw(m);
