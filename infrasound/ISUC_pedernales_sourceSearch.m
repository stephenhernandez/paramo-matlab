clear; close all; clc;

dayStart = datetime(2016,04,17);
dayEnd = datetime(2016,04,17);

dInc = 1;
dayVec = (dayStart:dInc:dayEnd)';
lDays = length(dayVec);

maxN = 1e4;
nComps = 6;
badComponents = []; %[2;5];
shiftsMaster = NaN(nComps-length(badComponents),maxN);
meanCCsMaster = NaN(maxN,1);
medCCsMaster = meanCCsMaster;
ampMaster = meanCCsMaster;
medratioMaster = meanCCsMaster;
stdratioMaster = meanCCsMaster;

tMaster = NaT(maxN,1);

newFs = 100;
Fs2 = 1e3;
secDur = 8; % seconds
nOverlap = 7/8;
winlen = round(secDur*newFs);
stride = round(winlen*(1-nOverlap))/newFs;

thresh = 0.045;
lfc = 0.2;
hfc = 0.8;
sensitivity = 1/3600;
nStrides = 2;
ampThresh = 1e-5;

%
refEllipse = referenceEllipsoid('wgs84');

% ch1	-0.64299	-78.49678	3683
% ch2	-0.64292	-78.49698	3682
% ch3	-0.64293	-78.49722	3680
% ch4	-0.64247	-78.49693	3686
% ch5	-0.64283	-78.49684	3685
% ch6	-0.6429	-78.49692	3683
kstnm = "ISUC";
knetwk = "EC";
kcmpnm = "BDF";
khole_list = ["01";"02";"03";"04";"05";"06"];

elementLats = [-0.64299; -0.64292; -0.64293; -0.64247; -0.64283; -0.6429];
elementLons = [-78.49678; -78.49698; -78.49722; -78.49693; -78.49684; -78.49692];
elementElevs = [3683; 3682; 3680; 3686; 3685; 3683];

elementLats(badComponents) = [];
elementLons(badComponents) = [];
elementElevs(badComponents) = [];

lS = length(elementLats);
waveformMaster = NaN(winlen*lS,maxN);
dx = [];
dy = [];
dz = [];
for i = 1:lS-1
    stla = elementLats(i);
    stlo = elementLons(i);
    [d_,az] = distance(stla,stlo,elementLats(i+1:end),elementLons(i+1:end),refEllipse);
    dx = [dx; d_.*cosd(90-az)];
    dy = [dy; d_.*sind(90-az)];
    dz = [dz; elementElevs(i+1:end) - elementElevs(i)];
    %
    disp([elementLats(i),elementLons(i)])
    disp([d_ az])
end
G = [dx dy]; % dz];
[Grows,Gcolumns] = size(G);
plagsMaster = NaN(Grows,maxN);
slownessMaster = NaN(Gcolumns,maxN);
%G = abs([dx dy]);

%%
n = 0;
for jj = 1:lDays
    day_ = dayVec(jj);
    S = loadWaveforms(day_,dInc,kstnm,kcmpnm,knetwk,khole_list); %this order cannot change
    Sorig = S;
    %S = detrendWaveforms(cutWaveforms(S,dateshift(S(1).ref,'start','day')+hours(00)+minutes(05),0,minutes(15))); %%<-- unique to ISUC

    if any(isnat(pull(S,'ref')))
        fprintf(2,'no data for: %s\n',datestr(day_));
        continue;
    end

    S(badComponents) = [];

    %
    fprintf(1,'processing: %s\n',datestr(day_));
    Scut = intWaveforms(filterWaveforms(detrendWaveforms(differentiateWaveforms(S)),lfc,hfc)); fprintf('done filtering...\n');
    %Scut = cutWaveforms(Scut,datetime(2021,12,25,10,00,00),0,hours(2));
    Scut = syncWaveforms(resampleWaveforms(Scut,newFs));
    dstack = [];

    lS = length(S);
    Fs = round(1./median(pull(Scut,'delta')));
    if Fs ~= newFs
        fprintf('something went wrong with resampling step\n');
        continue;
    end

    for i = 1:lS
        d = double(pull(Scut(i)));
        [dcut,~,endIndex] = cutWindows(d,winlen,nOverlap,true);
        dcut = taper(dcut,0.1);
        dstack = [dstack; dcut];
        disp(i);
    end

    t_ = getTimeVec(Scut(1));
    tdumcut = t_(endIndex);
    clear t_;
    nWindows = size(dstack,2);
    si = 1 + winlen*(0:lS-1)';
    ei = winlen*(1:lS)';

    shifts = NaN(lS,nWindows);
    meanCCs = NaN(nWindows,1);
    medCCs = meanCCs;
    medratio = meanCCs;
    stdratio = meanCCs;

    %
    dispIndex = [];
    prcWins = round(100*(1:nWindows)/nWindows);
    for i = [1 10 20 30 40 50 60 70 80 90]
        dispIndex = [dispIndex; find(prcWins>i,1)];
    end

    %
    Ndiff = 0.5*lS*(lS-1);
    plags = NaN(Ndiff,nWindows);
    slownii = NaN(size(G,2),nWindows);
    for i = 1:nWindows
        dslice = dstack(:,i);
        dslice = reshape(dslice,[winlen,lS]);
        [dslice,maxccp_,~,plags_,raw_shifts,meancc,medcc] = apply_vdcc(dslice,[],false,false,false);
        shifts(:,i) = raw_shifts;
        meanCCs(i) = meancc;
        medCCs(i) = medcc;
        plags(:,i) = plags_(1:end-1);

        % compute slowness vector
        %raw_shifts = ((G'*W*G)^(-1))*G'*W*plags_;
        plags_ = -plags_(1:end-1)/newFs;
        %m = ((G'*G).^(-1))*G'*plags_;
        m = plags_\G;
        slownii(:,i) = m;

        % super dirty, must be a way to speed up
        medampratios = NaN(1,Ndiff);
        dslice = abs(dslice);
        refBlock = dslice; %(:,2:end);

        ri = 1;
        ln = size(refBlock,2)-1;
        for j = 1:lS-1
            d1 = dslice(:,j);
            refBlock = circshift(refBlock,-1,2);
            refBlock = refBlock(:,1:end-1);
            ln = size(refBlock,2);
            medampratios(ri:ri+ln-1) = median(log10(d1./refBlock),'omitnan');
            ri = ri + ln;
        end
        medratio(i) = median(medampratios,'omitnan');
        stdratio(i) = mad(medampratios,1,2); %'omitnan');

        %
        if ismember(i,dispIndex)
            fprintf('%g\n',floor(100*i/nWindows));
        end
    end

    %figure('units','normalized','outerposition',[0 0 1 1]); plot(t(endIndex),medCCs,'.'); %,100,0.5*peak2peak(dstack)'/3600,'filled'); zoom on; SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; SS.MarkerFaceAlpha = 0.5; colorbar;
    tmpAmps = 0.5*peak2peak(dstack)';
    mI = medCCs >= thresh & tmpAmps >= ampThresh/sensitivity;
    summi = sum(mI);
    if ~summi
        continue;
    end

    %
    t = getTimeVec(Scut);
    mI = find(mI);
    if summi > 1
        ccGood = medCCs(mI);
        tGood = t(endIndex(mI));
        tdiff = seconds(diff(tGood));
        lstride = tdiff <= nStrides*stride;
        sumlstride = sum(lstride);
        while sumlstride
            lstride = find(lstride,1);
            i1 = lstride(1);
            i2 = i1 + 1;
            cc1 = ccGood(i1);
            cc2 = ccGood(i2);
            if cc1 > cc2
                %fprintf('deleting: (%s,%g), keeping (%s,%g)\n',datestr(tGood(i2)),cc2,datestr(tGood(i1)),cc1);
                mI(i2) = [];
            else
                %fprintf('deleting: (%s,%g), keeping (%s,%g)\n',datestr(tGood(i1)),cc1,datestr(tGood(i2)),cc2);
                mI(i1) = [];
            end
            ccGood = medCCs(mI);
            tGood = t(endIndex(mI));
            tdiff = seconds(diff(tGood));
            lstride = tdiff <= nStrides*stride;
            sumlstride = sum(lstride);
        end
        summi = length(mI);
    end

    %
    n = n + 1;
    dstack = dstack(:,mI);
    shifts = shifts(:,mI);
    meanCCs = meanCCs(mI);
    medCCs = medCCs(mI);
    plags = plags(:,mI);
    waveformMaster(:,n:n+summi-1) = dstack;

    fprintf('resampling %d promising candidates...\n',summi);
    for jk = 1:summi
        dslice = dstack(:,jk);
        dslice = resample(dslice,Fs2,newFs);
        dslice = reshape(dslice,[Fs2*winlen/newFs,lS]);

        [~,~,~,plags_,raw_shifts,meancc,medcc] = apply_vdcc(dslice,[],false,false,false);
        shifts(:,jk) = raw_shifts;
        meanCCs(jk) = meancc;
        medCCs(jk) = medcc;
        plags(:,jk) = plags_(1:end-1);
    end

    shiftsMaster(:,n:n+summi-1) = shifts;
    ampMaster(n:n+summi-1) = 0.5*peak2peak(dstack)';
    meanCCsMaster(n:n+summi-1) = meanCCs;
    medCCsMaster(n:n+summi-1) = medCCs;
    tMaster(n:n+summi-1) = t(endIndex(mI));
    medratioMaster(n:n+summi-1) = medratio(mI);
    stdratioMaster(n:n+summi-1) = stdratio(mI);
    slownessMaster(:,n:n+summi-1) = slownii(:,mI);
    plagsMaster(:,n:n+summi-1) = plags;
    n = n+summi-1;
end

%%
shiftsMaster(:,n+1:end) = [];
ampMaster(n+1:end) = [];
meanCCsMaster(n+1:end) = [];
medCCsMaster(n+1:end) = [];
medratioMaster(n+1:end) = [];
stdratioMaster(n+1:end) = [];
slownessMaster(:,n+1:end) = [];
waveformMaster(:,n+1:end) = [];
plagsMaster(:,n+1:end) = [];
tMaster(n+1:end) = [];

%%
close all;
lS = nComps-length(badComponents);
n = size(tMaster,1);
for i = 1:lS
    figure(1);
    hold on; ll = plot(sort(shiftsMaster(i,:))/Fs2,(0:n-1)'/n,'linewidth',2);
    ll.Color(4) = 0.8; zoom on;
    fprintf('Sensor %d: Averge Shift: %g\n',i,median(shiftsMaster(i,:))/Fs2);
end

%
figure('units','normalized','outerposition',[0 0 1 1]);
SS = scatter(tMaster,ampMaster*sensitivity,100,medCCsMaster,'filled');
zoom on; grid on; ax = gca; ax.YScale = 'log';
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; SS.MarkerFaceAlpha = 0.5; colorbar; caxis([0.5 1])

%
shiftsMaster = shiftsMaster - shiftsMaster(nComps-length(badComponents),:);
plagsMaster = getDD(shiftsMaster);

vel = NaN(size(medCCsMaster));
baz = vel;
bel = vel;

for i = 1:size(plagsMaster,2)
    %
    pm = -plagsMaster(:,i)/Fs2;
    %slow_ = lscov(G,pm);
    %slow_ = pinv(G)*pm;
    slow_ = pm\G;
    %slow_ = mean(G./pm);
    slownessMaster(:,i) = slow_';
    %

    if size(G,2) > 2
        [azimuth,elevation,r] = cart2sph(slow_(1),slow_(2),slow_(3));
        vel(i,1) = r; %rssq(slow_);
        baz(i,1) = 90 - rad2deg(azimuth);
        bel(i,1) = rad2deg(elevation);
        if bel(i) > 90
            bel(i) = bel(i) - 90;
        end
        fprintf('%g %g %g\n',vel(i),baz(i),bel(i));
    else
        baz(i,1) = atan2d(slow_(1),slow_(2))+360;
        if baz(i) >360
            baz(i) = baz(i)-360;
        end
        vel(i,1) = rssq(slow_);
    end

    if baz(i) > 180
        baz(i) = baz(i) - 360;
    end

end
strictI = medCCsMaster >= 0.55;
% & vel >= 230 & vel <= 430 & ampMaster*sensitivity >= 2*ampThresh & baz >= -65 & baz <= -5; % & bel <= 4 & baz >= -65 & baz <= -35;

figure('units','normalized','outerposition',[0 0 1/2 1]);

for i = 1:lS
    ax(i) = subplot(lS,1,i);
    plot(ax(i),tMaster(strictI),shiftsMaster(i,strictI)/Fs2,'.');
    zoom on; grid on; ylabel('lag [sec.]'); title(['Node: ',num2str(i)]);
end
linkaxes(ax,'xy'); zoom on;

%
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tMaster(strictI),baz(strictI),'o'); zoom on; grid on; ylim([-180 180]);
hold on;
plot(tMaster(~strictI),baz(~strictI),'.'); zoom on; grid on; ylim([-180 180]);

%
figure('units','normalized','outerposition',[0 0 1 1]);
ax9 = subplot(211);
SS = scatter(tMaster(strictI),ampMaster(strictI)*sensitivity,100,baz(strictI),'filled');
zoom on; grid on; ax = gca; ax.YScale = 'log';
SS.MarkerEdgeColor = 'k'; SS.MarkerEdgeAlpha = 0.5; SS.MarkerFaceAlpha = 0.5; colorbar; caxis([-180 180]);
ylabel('Presion [Pa.]');

ax9(2) = subplot(212);
histogram(tMaster(strictI),(dateshift(min(tMaster(strictI)),'start','day'):hours(1):dateshift(max(tMaster(strictI)),'end','day'))'); zoom on; grid on;
ylabel('Numero cada Hora');
cbar_dummy = colorbar;
cbar_dummy.Visible = 'off';
linkaxes(ax9,'x');

if size(G,2) > 2
    figure();
    plot(baz(strictI),bel(strictI),'o'); zoom on;
    hold on; plot(baz(~strictI),bel(~strictI),'.'); zoom on; grid on;
end

figure();
plot(tMaster(strictI),(vel(strictI))','o'); zoom on;
hold on;
plot(tMaster(~strictI),(vel(~strictI))','.'); zoom on;

tGood = tMaster(strictI);
aGood = ampMaster(strictI)*sensitivity;

figure();
clear sax;
sax(1) = subplot(211);
plot(tGood,(0:sum(strictI)-1)'/sum(strictI),'.'); zoom on; grid on;
sax(2) = subplot(212); 
histogram(tGood,dateshift(min(tGood),'start','day'):days(1):dateshift(max(tGood),'end','day')+1); zoom on; grid on; linkaxes(sax,'x');

Ncut = 101;
[~,startIndex,endIndex,~,~] = cutWindows(datenum(tGood),Ncut,Ncut-1,false);
aCut = cutWindows(aGood,Ncut,Ncut-1,false);

figure();
plot(tGood(endIndex),mean(aCut)'.*3600*Ncut./seconds(tGood(endIndex)-tGood(startIndex)),'.'); zoom on; grid on;
%close(2);